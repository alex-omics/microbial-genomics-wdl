version 1.0

# kraken2 over a batch of samples on one VM. The database is extracted and read
# once, then shared by several concurrent kraken2 processes through the page cache
# (--memory-mapping), rather than being copied into RAM once per sample.
#
# viral-classify ships kraken2, samtools and python3 in one image; the kraken2
# images from other sources lack samtools.

task check_read_inputs {

    input {
        Int     n_bams
        Int     n_fastq_r1
        Int     n_fastq_r2
        String  docker = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        n_bams:     "Length of reads_bams (0 if not supplied)"
        n_fastq_r1: "Length of reads_fastq_r1 (0 if not supplied)"
        n_fastq_r2: "Length of reads_fastq_r2 (0 if not supplied)"
        docker:     "Container image"
    }

    meta {
        description: "Fail within seconds if the read inputs are contradictory, before any database is downloaded. Emits the number of samples."
    }

    command <<<
        set -euo pipefail
        if [ ~{n_bams} -gt 0 ] && [ ~{n_fastq_r1} -gt 0 ]; then
            echo "ERROR: supply reads_bams or reads_fastq_r1, not both." >&2; exit 1
        fi
        if [ ~{n_bams} -eq 0 ] && [ ~{n_fastq_r1} -eq 0 ]; then
            echo "ERROR: supply reads_bams, or reads_fastq_r1 (with reads_fastq_r2 if paired)." >&2; exit 1
        fi
        if [ ~{n_bams} -gt 0 ] && [ ~{n_fastq_r2} -gt 0 ]; then
            echo "ERROR: reads_fastq_r2 goes with reads_fastq_r1, not reads_bams." >&2; exit 1
        fi
        if [ ~{n_fastq_r2} -gt 0 ] && [ ~{n_fastq_r2} -ne ~{n_fastq_r1} ]; then
            echo "ERROR: reads_fastq_r1 has ~{n_fastq_r1} files but reads_fastq_r2 has ~{n_fastq_r2}." >&2; exit 1
        fi
        echo $(( ~{n_bams} + ~{n_fastq_r1} )) > n_samples
    >>>

    output {
        Int n_samples = read_int("n_samples")
    }

    runtime {
        docker:         docker
        memory:         "1 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     2
    }
}


task kraken2_batch {

    input {
        Array[String]   names
        Array[String]   bam_filenames
        Array[String]   report_filenames
        Array[File]     reads
        Array[File]     reads_r2            = []
        Boolean         bam_input           = false
        File            db_tgz
        Array[Int]      filter_taxids       = []
        Array[String]   filter_taxon_names  = []
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
        names:               "Sample names, positionally matched to reads. Letters, digits, dot, underscore and hyphen only."
        bam_filenames:       "Output BAM filename per sample (e.g. name.bam), positionally matched to names. Given as an input so the output paths are known before the task runs."
        report_filenames:    "Output kraken2 report filename per sample, positionally matched to names"
        reads:               "Unaligned BAMs (bam_input) or FASTQ read 1 files (plain or gzipped)"
        reads_r2:            "FASTQ read 2 files, positionally matched to reads. Empty for BAMs and single-end FASTQ."
        bam_input:           "reads are BAMs. Paired-end is then detected from the read flags; otherwise reads are treated as single-end."
        db_tgz:              "Kraken2 database tarball containing hash.k2d, opts.k2d and taxo.k2d"
        filter_taxids:       "Taxids whose clades (the taxon and everything beneath it) are kept or removed"
        filter_taxon_names:  "Taxon names, matched exactly against the names in each sample's kraken2 report, with the same effect as filter_taxids"
        remove_matching:     "Drop the reads within the listed taxa instead of keeping them. Unclassified reads are then retained."
        keep_unclassified:   "When keeping, also retain unclassified reads. Has no effect with remove_matching."
        confidence:          "kraken2 --confidence threshold"
        extra_args:          "Additional kraken2 arguments, e.g. '--minimum-hit-groups 3'"
        cpu:                 "CPUs on the VM, shared across the concurrent kraken2 processes"
        mem_gb:              "Memory in GB. Must exceed the uncompressed hash.k2d so the shared cache stays resident."
        concurrent:          "Samples classified at once; each gets cpu / concurrent threads"
        disk_gb:             "Disk in GB. Default is sized from the database tarball and the input reads."
        preemptible:         "Preemptible attempts. A preemption reruns the whole batch, so the default is 0."
        docker:              "Container image with kraken2, samtools and python3"
    }

    meta {
        description: "Classify a batch of samples with kraken2 against one shared copy of the database, and write one unaligned BAM per sample holding either all of its reads or, if taxa are given, only those kept or removed by clade. Emits a report and a summary row per sample, all in input order."
    }

    Int threads = if (cpu / concurrent) < 1 then 1 else (cpu / concurrent)
    Boolean paired_fastq = length(reads_r2) > 0
    Float input_gb = size(reads, "GB") + size(reads_r2, "GB")
    # Database tarball plus its extraction, every input localized up front, the output
    # BAMs, and the FASTQ and per-read output of the samples in flight.
    Int disk = select_first([disk_gb, ceil(size(db_tgz, "GB") * 3 + input_gb * 3 + 5 * concurrent * input_gb / length(reads)) + 50])
    Boolean filtering = length(filter_taxids) + length(filter_taxon_names) > 0

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

        # name, read 1 or BAM, read 2 ("-" when there is none), output BAM, output report
        if [ "~{paired_fastq}" = "true" ]; then
            paste ~{write_lines(names)} ~{write_lines(reads)} ~{write_lines(reads_r2)} > manifest.tsv
        else
            paste ~{write_lines(names)} ~{write_lines(reads)} | awk '{print $0 "\t-"}' > manifest.tsv
        fi
        paste manifest.tsv ~{write_lines(bam_filenames)} ~{write_lines(report_filenames)} > manifest.full.tsv
        mv manifest.full.tsv manifest.tsv
        awk '{print "out/" $1 ".summary.tsv"}' manifest.tsv > summaries.list

        # Applies the taxon filter to one sample's per-read assignments and writes
        # the read names to keep and the sample's summary row.
        cat > summarize.py <<'PY'
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("report")
ap.add_argument("assignments")
ap.add_argument("names_out")
ap.add_argument("summary_out")
ap.add_argument("--taxids", default="")
ap.add_argument("--taxon-names", default="")
ap.add_argument("--mode", default="keep")
ap.add_argument("--keep-unclassified", action="store_true")
a = ap.parse_args()

want = {int(t) for t in a.taxids.split(",") if t}
want_names = {n.strip() for n in open(a.taxon_names)} - {""} if a.taxon_names else set()

# Report rows: pct, clade reads, direct reads, [extra columns], rank, taxid, name.
# A taxon's descendants are the rows after it with deeper indentation.
rows = []
for ln in open(a.report):
    f = ln.rstrip("\n").split("\t")
    rows.append((int(f[-2]), len(f[-1]) - len(f[-1].lstrip(" ")), int(f[1]), f[-3], f[-1].strip()))

by_name = {}
for r in rows:
    by_name.setdefault(r[4], []).append(r[0])
for n in sorted(want_names - set(by_name)):
    print(f"{a.name}: no taxon named '{n}' in the report")
want |= {t for n in want_names for t in by_name.get(n, [])}

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

filtering = bool(want or want_names)
pct = lambda n: f"{100 * n / total:.3f}" if total else "NA"
row = [a.name, total, total - unclassified, unclassified,
       top[0] if top else "NA", top[4] if top else "NA", top[2] if top else "NA",
       pct(top[2]) if top else "NA",
       kept if filtering else "NA", pct(kept) if filtering else "NA",
       total - kept if filtering else "NA"]
with open(a.summary_out, "w") as out:
    out.write("\t".join(str(x) for x in row) + "\n")
PY

        # Writes the records of a FASTQ whose read name is in the keep list.
        cat > filter_fastq.py <<'PY'
import gzip
import sys

keep_file, fastq, out_fastq = sys.argv[1:4]


def norm(name):
    name = name.split()[0]
    return name[:-2] if name.endswith(("/1", "/2")) else name


keep = {norm(ln) for ln in open(keep_file) if ln.strip()}
with open(fastq, "rb") as fh:
    opener = gzip.open if fh.read(2) == b"\x1f\x8b" else open
with opener(fastq, "rt") as fi, open(out_fastq, "w") as fo:
    while True:
        head = fi.readline()
        if not head:
            break
        rec = head + fi.readline() + fi.readline() + fi.readline()
        if norm(head[1:]) in keep:
            fo.write(rec)
PY

        cat > process_one.sh <<'SH'
#!/bin/bash
set -euo pipefail
name="$1"; in1="$2"; in2="$3"; bam_out="out/$4"; report_out="out/$5"
w="tmp/${name}"; mkdir -p "${w}"

# Read files for kraken2: converted from the BAM, or the FASTQ as given.
if [ "${BAM}" = "true" ]; then
    if [ "$(samtools view -c -f 1 "${in1}")" -gt 0 ]; then
        samtools collate -u -O "${in1}" "${w}/collate" \
            | samtools fastq -n -1 "${w}/r1.fq" -2 "${w}/r2.fq" -0 /dev/null -s /dev/null -
        r1="${w}/r1.fq"; r2="${w}/r2.fq"
    else
        samtools fastq -n "${in1}" > "${w}/r.fq"
        r1="${w}/r.fq"; r2="-"
    fi
else
    r1="${in1}"; r2="${in2}"
fi
if [ "${r2}" = "-" ]; then reads=("${r1}"); else reads=(--paired "${r1}" "${r2}"); fi
gz=""
if [ "$(head -c2 "${r1}" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then gz="--gzip-compressed"; fi

# kraken2 writes no output file for an empty input
: > "${w}/assignments.tsv"
kraken2 --db "${DB_DIR}" --memory-mapping --threads "${THREADS}" --confidence "${CONF}" \
    ${gz} ${EXTRA} --report "${report_out}" --output "${w}/assignments.tsv" \
    "${reads[@]}" 2> "${w}/kraken2.log"
if [ "${BAM}" = "true" ]; then rm -f "${w}"/*.fq; fi

python3 summarize.py "${name}" "${report_out}" "${w}/assignments.tsv" \
    "${w}/keep.txt" "out/${name}.summary.tsv" --taxids "${TAXIDS}" --mode "${MODE}" \
    ${NAMES_ARG} ${KEEPUNC}
rm -f "${w}/assignments.tsv"

# Output BAM: every read, or only those kept by the taxon filter.
if [ "${BAM}" = "true" ]; then
    if [ "${FILTER}" = "true" ]; then
        samtools view -b -@ "${THREADS}" -N "${w}/keep.txt" -o "${bam_out}" "${in1}"
    else
        cp "${in1}" "${bam_out}"
    fi
else
    if [ "${FILTER}" = "true" ]; then
        python3 filter_fastq.py "${w}/keep.txt" "${in1}" "${w}/k1.fq"
        r1="${w}/k1.fq"
        if [ "${r2}" != "-" ]; then
            python3 filter_fastq.py "${w}/keep.txt" "${in2}" "${w}/k2.fq"
            r2="${w}/k2.fq"
        fi
    fi
    # samtools import rejects an empty FASTQ, so no reads gives a header-only BAM
    if [ "${FILTER}" = "true" ]; then n_out="$(wc -l < "${w}/keep.txt")"; else n_out="$(cut -f2 "out/${name}.summary.tsv")"; fi
    if [ "${n_out}" -eq 0 ]; then
        printf '@HD\tVN:1.6\tSO:unsorted\n@RG\tID:%s\tSM:%s\n' "${name}" "${name}" | samtools view -b -o "${bam_out}" -
    else
        if [ "${r2}" = "-" ]; then src=(-0 "${r1}"); else src=(-1 "${r1}" -2 "${r2}"); fi
        samtools import -r "ID:${name}" -r "SM:${name}" "${src[@]}" -o "${bam_out}"
    fi
fi
rm -rf "${w}"
SH
        chmod +x process_one.sh

        export THREADS=~{threads} CONF=~{confidence} EXTRA="~{extra_args}" BAM=~{bam_input}
        export TAXIDS="~{sep=',' filter_taxids}" MODE="~{if remove_matching then 'remove' else 'keep'}"
        export KEEPUNC="~{if keep_unclassified then '--keep-unclassified' else ''}"
        export FILTER=~{filtering}
        export NAMES_ARG="~{if length(filter_taxon_names) > 0 then '--taxon-names ' + write_lines(filter_taxon_names) else ''}"

        # A failed sample stops the batch at once (exit 255 halts xargs) instead of
        # letting the remaining samples run on a batch that will be retried anyway.
        run_one() { ./process_one.sh "$@" || { echo "FAILED: $1" >&2; exit 255; }; }
        export -f run_one
        xargs -P ~{concurrent} -L 1 bash -c 'run_one "$@"' _ < manifest.tsv

        {
            printf 'sample\ttotal_reads\tclassified_reads\tunclassified_reads\ttop_species_taxid\ttop_species\ttop_species_reads\ttop_species_pct\tkept_reads\tkept_pct\tremoved_reads\n'
            while read -r f; do cat "${f}"; done < summaries.list
        } > batch_summary.tsv
    >>>

    output {
        File          batch_summary = "batch_summary.tsv"
        # Declared from the inputs, not read from a file: the backend must know the
        # output paths before the task runs to copy them out of the VM.
        Array[File]   reports       = prefix("out/", report_filenames)
        Array[File]   bams          = prefix("out/", bam_filenames)
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


task build_manifest {

    input {
        Array[String]  sample_ids
        Array[String]  bam_paths
        Array[String]  report_paths
        File           summary_tsv
        String         table_name = "sample"
        String         basename   = "kraken2_manifest"
        String         docker     = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        sample_ids:   "Sample names, in input order"
        bam_paths:    "Cloud paths of the output BAMs, positionally matched to sample_ids"
        report_paths: "Cloud paths of the kraken2 reports, positionally matched to sample_ids"
        summary_tsv:  "Merged per-sample summary table"
        table_name:   "Terra data table the manifest is for. The first column is entity:<table_name>_id."
        basename:     "Basename for the emitted TSV, without extension"
        docker:       "Container image"
    }

    meta {
        description: "Join each sample's output BAM and report paths with its summary statistics into one TSV that can be uploaded to a Terra data table, attaching the outputs to the per-sample rows. Every column but the ID is prefixed with the workflow name, so its origin is clear in the table."
    }

    command <<<
        set -euo pipefail
        paste ~{write_lines(sample_ids)} ~{write_lines(bam_paths)} ~{write_lines(report_paths)} > paths.tsv

        awk -F'\t' -v OFS='\t' -v t="~{table_name}" '
            NR == FNR {
                rest = $0; sub(/^[^\t]*\t/, "", rest)
                if (FNR == 1) hdr = rest; else stats[$1] = rest
                next
            }
            FNR == 1 {
                n = split(hdr, h, "\t")
                for (i = 1; i <= n; i++) h[i] = "classify_kraken2_" h[i]
                cols = h[1]; for (i = 2; i <= n; i++) cols = cols OFS h[i]
                print "entity:" t "_id", "classify_kraken2_bam", "classify_kraken2_report", cols
            }
            { print $1, $2, $3, stats[$1] }
        ' ~{summary_tsv} paths.tsv > "~{basename}.tsv"
    >>>

    output {
        File manifest = "~{basename}.tsv"
    }

    runtime {
        docker:         docker
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     2
    }
}
