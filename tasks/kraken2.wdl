version 1.0

# kraken2 over a batch of samples on one VM. The database is extracted and read
# once, then shared by several concurrent kraken2 processes through the page cache
# (--memory-mapping), rather than being copied into RAM once per sample.
#
# viral-classify ships kraken2, samtools and python3 in one image; the kraken2
# images from other sources lack samtools.

task kraken2_batch {

    input {
        Array[String]   names
        Array[File]     bams
        File            db_tgz
        Array[Int]      filter_taxids       = []
        Boolean         remove_matching     = false
        Boolean         keep_unclassified   = false
        Float           confidence          = 0.0
        String          extra_args          = ""
        Int             cpu                 = 32
        Int             mem_gb              = 128
        Int             concurrent          = 8
        Int?            disk_gb
        Int             preemptible         = 0
        String          docker              = "quay.io/broadinstitute/viral-classify:2.5.21.0"
    }

    parameter_meta {
        names:             "Sample names, positionally matched to bams. Letters, digits, dot, underscore and hyphen only."
        bams:              "Unaligned BAMs. Paired-end is detected from the read flags; otherwise reads are treated as single-end."
        db_tgz:            "Kraken2 database tarball containing hash.k2d, opts.k2d and taxo.k2d"
        filter_taxids:     "Taxids whose clades (the taxon and everything beneath it) are kept or removed. Empty means classify only."
        remove_matching:   "Drop the reads within filter_taxids instead of keeping them. Unclassified reads are then retained."
        keep_unclassified: "When keeping, also retain unclassified reads. Has no effect with remove_matching."
        confidence:        "kraken2 --confidence threshold"
        extra_args:        "Additional kraken2 arguments, e.g. '--minimum-hit-groups 3'"
        cpu:               "CPUs on the VM, shared across the concurrent kraken2 processes"
        mem_gb:            "Memory in GB. Must exceed the uncompressed hash.k2d so the shared cache stays resident."
        concurrent:        "Samples classified at once; each gets cpu / concurrent threads"
        disk_gb:           "Disk in GB. Default is sized from the database tarball and the input BAMs."
        preemptible:       "Preemptible attempts. A preemption reruns the whole batch, so the default is 0."
        docker:            "Container image with kraken2, samtools and python3"
    }

    meta {
        description: "Classify a batch of unaligned BAMs with kraken2 against one shared copy of the database, and optionally keep or remove the reads that fall within given clades. Emits a report, a summary row and, when filtering, a filtered BAM per sample, all in input order."
    }

    Int threads = if (cpu / concurrent) < 1 then 1 else (cpu / concurrent)
    Float bam_gb = size(bams, "GB")
    # Database tarball plus its extraction, every BAM localized up front, the filtered
    # BAMs, and the FASTQ and per-read output of the samples in flight.
    Int disk = select_first([disk_gb, ceil(size(db_tgz, "GB") * 3 + bam_gb * 2 + 5 * concurrent * bam_gb / length(bams)) + 50])
    Boolean filtering = length(filter_taxids) > 0

    command <<<
        set -euo pipefail
        mkdir -p db out tmp

        tar -xf "~{db_tgz}" -C db
        HASH="$(find db -name hash.k2d -print -quit)"
        test -n "${HASH}"
        DB_DIR="$(dirname "${HASH}")"
        export DB_DIR

        # A database larger than memory still runs, but every read pages from disk.
        DB_KB=$(( $(stat -c %s "${HASH}") / 1024 ))
        MEM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo)
        if [ "${DB_KB}" -gt "${MEM_KB}" ]; then
            echo "WARNING: hash.k2d (${DB_KB} kB) exceeds memory (${MEM_KB} kB); expect heavy paging." >&2
        fi
        cat "${DB_DIR}/hash.k2d" "${DB_DIR}/taxo.k2d" > /dev/null

        paste ~{write_lines(names)} ~{write_lines(bams)} > manifest.tsv
        awk '{print "out/" $1 ".kraken2.report.txt"}' manifest.tsv > reports.list
        awk '{print "out/" $1 ".summary.tsv"}' manifest.tsv > summaries.list
        awk '{print "out/" $1 ".filtered.bam"}' manifest.tsv > filtered.list

        # Applies the taxon filter to one sample's per-read assignments and writes
        # the QNAMEs to keep and the sample's summary row.
        cat > summarize.py <<'PY'
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("report")
ap.add_argument("assignments")
ap.add_argument("names_out")
ap.add_argument("summary_out")
ap.add_argument("--taxids", default="")
ap.add_argument("--mode", default="keep")
ap.add_argument("--keep-unclassified", action="store_true")
a = ap.parse_args()

want = {int(t) for t in a.taxids.split(",") if t}

# Report rows: pct, clade reads, direct reads, [extra columns], rank, taxid, name.
# A taxon's descendants are the rows after it with deeper indentation.
rows = []
for ln in open(a.report):
    f = ln.rstrip("\n").split("\t")
    rows.append((int(f[-2]), len(f[-1]) - len(f[-1].lstrip(" ")), int(f[1]), f[-3], f[-1].strip()))

selected = set()
for i, (tid, depth, _, _, _) in enumerate(rows):
    if tid in want:
        selected.add(tid)
        j = i + 1
        while j < len(rows) and rows[j][1] > depth:
            selected.add(rows[j][0])
            j += 1
missing = sorted(want - selected)
if missing:
    print(f"{a.name}: no reads at or under taxid(s) {missing}")

top = max((r for r in rows if r[3] == "S"), key=lambda r: r[2], default=None)

total = unclassified = kept = 0
with open(a.assignments) as fh, open(a.names_out, "w") as out:
    for ln in fh:
        c = ln.split("\t", 3)
        total += 1
        if c[0] == "U":
            unclassified += 1
        if a.mode == "keep":
            keep = int(c[2]) in selected or (c[0] == "U" and a.keep_unclassified)
        else:
            keep = int(c[2]) not in selected
        if keep:
            kept += 1
            out.write(c[1] + "\n")

filtering = bool(want)
pct = lambda n: f"{100 * n / total:.3f}" if total else "NA"
row = [a.name, total, total - unclassified, unclassified,
       top[0] if top else "NA", top[4] if top else "NA", top[2] if top else "NA",
       pct(top[2]) if top else "NA",
       kept if filtering else "NA", pct(kept) if filtering else "NA"]
with open(a.summary_out, "w") as out:
    out.write("\t".join(str(x) for x in row) + "\n")
PY

        cat > process_one.sh <<'SH'
#!/bin/bash
set -euo pipefail
name="$1"; bam="$2"
w="tmp/${name}"; mkdir -p "${w}"

if [ "$(samtools view -c -f 1 "${bam}")" -gt 0 ]; then
    samtools collate -u -O "${bam}" "${w}/collate" \
        | samtools fastq -n -1 "${w}/r1.fq" -2 "${w}/r2.fq" -0 /dev/null -s /dev/null -
    reads=(--paired "${w}/r1.fq" "${w}/r2.fq")
else
    samtools fastq -n "${bam}" > "${w}/r.fq"
    reads=("${w}/r.fq")
fi

# kraken2 writes no output file for an empty BAM
: > "${w}/assignments.tsv"
kraken2 --db "${DB_DIR}" --memory-mapping --threads "${THREADS}" --confidence "${CONF}" \
    ${EXTRA} --report "out/${name}.kraken2.report.txt" --output "${w}/assignments.tsv" \
    "${reads[@]}" 2> "${w}/kraken2.log"
rm -f "${w}"/*.fq

python3 summarize.py "${name}" "out/${name}.kraken2.report.txt" "${w}/assignments.tsv" \
    "${w}/keep.txt" "out/${name}.summary.tsv" --taxids "${TAXIDS}" --mode "${MODE}" ${KEEPUNC}
rm -f "${w}/assignments.tsv"

if [ -n "${TAXIDS}" ]; then
    samtools view -b -@ "${THREADS}" -N "${w}/keep.txt" -o "out/${name}.filtered.bam" "${bam}"
fi
rm -rf "${w}"
SH
        chmod +x process_one.sh

        export THREADS=~{threads} CONF=~{confidence} EXTRA="~{extra_args}"
        export TAXIDS="~{sep=',' filter_taxids}" MODE="~{if remove_matching then 'remove' else 'keep'}"
        export KEEPUNC="~{if keep_unclassified then '--keep-unclassified' else ''}"

        # A failed sample stops the batch at once (exit 255 halts xargs) instead of
        # letting the remaining samples run on a batch that will be retried anyway.
        run_one() { ./process_one.sh "$@" || { echo "FAILED: $1" >&2; exit 255; }; }
        export -f run_one
        xargs -P ~{concurrent} -L 1 bash -c 'run_one "$@"' _ < manifest.tsv

        {
            printf 'sample\ttotal_reads\tclassified_reads\tunclassified_reads\ttop_species_taxid\ttop_species\ttop_species_reads\ttop_species_pct\tkept_reads\tkept_pct\n'
            while read -r f; do cat "${f}"; done < summaries.list
        } > batch_summary.tsv
    >>>

    output {
        File          batch_summary = "batch_summary.tsv"
        Array[File]   reports       = read_lines("reports.list")
        Array[File]   filtered_bams = if filtering then read_lines("filtered.list") else []
    }

    runtime {
        docker:         docker
        cpu:            cpu
        memory:         "~{mem_gb} GB"
        disks:          "local-disk ~{disk} SSD"
        preemptible:    preemptible
        maxRetries:     1
    }
}
