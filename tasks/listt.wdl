version 1.0

# Tasks for call_ospA_type: Pfizer's LISTT (Lyme In Silico Typing Tool, Lee et
# al. 2024, Microb Genom 10:001252) wrapped without modifying its code.
#
# LISTT's own code is run as shipped (docker/listt). What this file adds:
#   - staging inputs under a name LISTT can parse (it derives the sample name by
#     splitting the filename at the first '_' or '.f')
#   - enforcing the per-IST identity thresholds in variants/thresholds.csv.
#     The paper assigns an IST to a non-exact match only if the identity meets
#     the threshold and otherwise calls "unknown IST"; the released code never
#     reads that file, so every sequence is forced into its nearest IST
#   - a trust cut-off (max_trusted_ist) for ISTs that are not in the paper
#   - tolerating LISTT crashing on a variant that is missing from its own table
#   - reading ospA locus-level quality, so bad data can be told from a novel allele
#   - reconciling an assembly call with a reads call into one answer and a flag

task validate_ospA_inputs {

    input {
        String  sample_name
        Boolean has_read1
        Boolean has_read2
        Int     max_trusted_ist
        String  docker = "aarvani1/listt:1.2.0@sha256:c016417e25c3e4b96d4fc98bedbf4d77e15608674e5e5320cdaff1df6a1a38b2"
    }

    parameter_meta {
        sample_name:     "Sample label, used only in outputs"
        has_read1:       "Whether read1 was supplied"
        has_read2:       "Whether read2 was supplied"
        max_trusted_ist: "Highest IST number reported as a call"
        docker:          "Container image"
    }

    meta {
        description: "Fail fast on a mis-specified sample before any alignment runs: reads must be supplied as a pair or not at all, because a lone read1 would otherwise be ignored without a word and the sample reported as assembly-only."
    }

    command <<<
        set -euo pipefail

        if [ "~{has_read1}" != "~{has_read2}" ]; then
            echo "ERROR: read1 and read2 must be supplied together, or both left empty." >&2
            echo "       LISTT's reads mode is built for Illumina paired-end data." >&2
            exit 1
        fi

        # Passed through a file so that no character in the name can alter the shell.
        NAME="$(cat ~{write_lines([sample_name])})"
        if [ -z "${NAME}" ]; then
            echo "ERROR: sample_name is empty." >&2
            exit 1
        fi
        if printf '%s' "${NAME}" | grep -qP '[\t\r\n]'; then
            echo "ERROR: sample_name contains a tab or newline, which would corrupt the results table." >&2
            exit 1
        fi

        if [ "~{max_trusted_ist}" -lt 0 ]; then
            echo "ERROR: max_trusted_ist must be zero or more." >&2
            exit 1
        fi

        echo "ok" > ok.txt
    >>>

    output {
        String validated = read_string("ok.txt")
    }

    runtime {
        docker:      docker
        cpu:         1
        memory:      "1 GB"
        disks:       "local-disk 1 SSD"
        preemptible: 1
    }
}

task bbsl_ospA_listt {

    input {
        String  sample_name
        String  mode                 # "assembly" or "reads"
        File?   assembly_fasta
        File?   read1
        File?   read2
        Int     max_trusted_ist      = 17
        String  validated            = "ok"   # wired to validate_ospA_inputs so it runs first

        Int     cpu                  = 2
        Int     mem_gb               = 4
        Int?    disk_gb
        Int     preemptible          = 1
        String  docker               = "aarvani1/listt:1.2.0@sha256:c016417e25c3e4b96d4fc98bedbf4d77e15608674e5e5320cdaff1df6a1a38b2"
    }

    parameter_meta {
        sample_name:     "Sample label, used only in outputs"
        mode:            "'assembly' (BLAST of a genome FASTA against the ospA references) or 'reads' (bwa alignment of Illumina paired-end reads)"
        assembly_fasta:  "Genome FASTA, required for mode='assembly'. Plain or gzip-compressed"
        read1:           "R1 FASTQ, required for mode='reads'. Plain or gzip-compressed"
        read2:           "R2 FASTQ, required for mode='reads'. LISTT's interleaved-FASTQ mode is not exposed"
        max_trusted_ist: "Highest IST number reported as a call (default = 17, the ISTs described in the paper). A closest match to a higher IST, or to a variant LISTT has not assigned to an IST, is reported as provisional and not as a call"
        validated:       "Ordering hook for validate_ospA_inputs; leave at its default"
        cpu:             "CPUs (default = 2: bwa and the samtools pipe run concurrently)"
        mem_gb:          "Memory in GB (default = 4)"
        disk_gb:         "Disk in GB. Default is 10 plus twice the size of the inputs"
        preemptible:     "Preemptible attempts (default = 1)"
        docker:          "Container image"
    }

    meta {
        description: "Run one mode of LISTT on one sample and classify its ospA consensus: closest reference variant, identity, the identity threshold for that variant's IST, and ospA locus quality. Emits a one-row TSV with a status of CALL, BELOW (identity under the IST threshold), PROVISIONAL (closest IST is above max_trusted_ist or not assigned) or QF (no usable ospA)."
    }

    Int disk = select_first([disk_gb, 10 + 2 * ceil(size(assembly_fasta, "GB") + size(read1, "GB") + size(read2, "GB"))])

    command <<<
        set -euo pipefail

        # validated="~{validated}": only here so this task waits for validate_ospA_inputs
        WD="${PWD}"
        STAGE="${WD}/stage"
        mkdir -p "${STAGE}"

        # LISTT takes the sample name from the filename: everything before the
        # first '_' for reads and before the first '.f' for a genome. A name
        # such as B31_S1.fasta would silently become "B31", so inputs are
        # staged under an ID that contains neither.
        SAMPLE="$(cat ~{write_lines([sample_name])})"
        SID="$(printf '%s' "${SAMPLE}" | tr -c 'A-Za-z0-9-' '-')"

        if [ "~{mode}" = "assembly" ]; then
            if [ -z "~{default='' assembly_fasta}" ]; then
                echo "ERROR: mode=assembly needs assembly_fasta." >&2; exit 1
            fi
            case "~{default='' assembly_fasta}" in
                *.gz) gzip -dc "~{default='' assembly_fasta}" > "${STAGE}/${SID}.fasta" ;;
                *)    cp  "~{default='' assembly_fasta}" "${STAGE}/${SID}.fasta" ;;
            esac
            LISTT_ARGS="-g ${SID}.fasta -m a -p ${STAGE}"
        elif [ "~{mode}" = "reads" ]; then
            if [ -z "~{default='' read1}" ] || [ -z "~{default='' read2}" ]; then
                echo "ERROR: mode=reads needs read1 and read2." >&2; exit 1
            fi
            ln -s "~{default='' read1}" "${STAGE}/${SID}_R1.fastq.gz"
            ln -s "~{default='' read2}" "${STAGE}/${SID}_R2.fastq.gz"
            LISTT_ARGS="-r1 ${SID}_R1.fastq.gz -r2 ${SID}_R2.fastq.gz -m r -p ${STAGE}"
        else
            echo "ERROR: mode must be 'assembly' or 'reads', not '~{mode}'." >&2; exit 1
        fi

        # LISTT reads ref/, variants/ and alleles/ relative to the working
        # directory, so it must run from its install directory. It also exits 0
        # when it fails and crashes outright on variant 68, which is missing
        # from its own serotypes.csv, so the exit status is recorded and not
        # trusted: the result is judged from the files it leaves behind.
        cd /opt/LISTT
        set +e
        python ospa_typing.py ${LISTT_ARGS} > "${WD}/listt.stdout" 2> "${WD}/listt.stderr"
        LISTT_EXIT=$?
        set -e
        cd "${WD}"

        cat > classify.py <<'PY'
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET

VARIANTS = "/opt/LISTT/variants"
(sample, sid, mode, stage, listt_exit, max_trusted,
 out_tsv, out_fasta, stdout_path, stderr_path) = sys.argv[1:11]
max_trusted = int(max_trusted)
aln = f"{stage}/aln/{sid}"
EPS = 1e-9
NA = "NA"


def read_fasta(path):
    records, name, seq = [], None, []
    if not os.path.isfile(path):
        return records
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


# ---- LISTT's own reference tables --------------------------------------
# serotypes.csv carries a BOM and CRLF line endings, repeats variant 107
# three times, and has no row for variant 68. Neither the repeat nor the gap
# is an error here: the first row for a variant wins, and a missing row means
# the variant is not assigned to any IST.
variant_info = {}
with open(f"{VARIANTS}/serotypes.csv", newline="", encoding="utf-8-sig") as fh:
    for row in csv.reader(fh):
        if len(row) < 3 or not row[0].strip().isdigit():
            continue
        variant_info.setdefault(int(row[0]), (row[1].strip(), row[2].strip()))

thresholds = {}
with open(f"{VARIANTS}/thresholds.csv", newline="", encoding="utf-8-sig") as fh:
    for row in csv.reader(fh):
        if len(row) >= 2 and row[0].strip().isdigit():
            thresholds[int(row[0])] = float(row[1])

ref_protein = {}
for name, seq in read_fasta(f"{VARIANTS}/ospA_protein.fas"):
    ref_protein[int(name.split("_")[-1])] = seq.upper().rstrip("*")

min_cov = {}
with open(f"{VARIANTS}/min_cov_metrics.csv", newline="", encoding="utf-8-sig") as fh:
    for row in csv.reader(fh):
        if len(row) >= 2:
            min_cov[row[0].strip()] = float(row[1])


def last_line(path):
    try:
        lines = [l.strip() for l in open(path, errors="replace") if l.strip()]
        return lines[-1][:200] if lines else ""
    except OSError:
        return ""


# ---- ospA locus quality ------------------------------------------------
def assembly_locus():
    """The BLAST hit LISTT itself used: the highest raw score, first one wins."""
    out = {"nt_len": 0, "nt_identity": NA, "contig": NA, "at_contig_edge": NA,
           "n_hit_contigs": 0, "ref": NA}
    path = f"{aln}/{sid}_blast.xml"
    if not os.path.isfile(path):
        return out
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return out
    best, contigs = None, set()
    for it in root.iter("Iteration"):
        contig = (it.findtext("Iteration_query-def") or "").split()[0]
        qlen = int(it.findtext("Iteration_query-len") or 0)
        for hit in it.iter("Hit"):
            ref = (hit.findtext("Hit_def") or hit.findtext("Hit_id") or "").split()[0]
            for hsp in hit.iter("Hsp"):
                score = float(hsp.findtext("Hsp_score"))
                if int(hsp.findtext("Hsp_align-len")) >= 400:
                    contigs.add(contig)
                if best is None or score > best[0]:
                    best = (score, contig, qlen, hsp, ref)
    if best is None:
        return out
    _, contig, qlen, hsp, ref = best
    qfrom, qto = int(hsp.findtext("Hsp_query-from")), int(hsp.findtext("Hsp_query-to"))
    out.update(
        nt_len=len((hsp.findtext("Hsp_qseq") or "").replace("-", "")),
        nt_identity=f"{100.0 * int(hsp.findtext('Hsp_identity')) / int(hsp.findtext('Hsp_align-len')):.2f}",
        contig=contig,
        at_contig_edge="yes" if min(qfrom, qto) <= 1 or max(qfrom, qto) >= qlen else "no",
        n_hit_contigs=len(contigs),
        ref=ref,
    )
    return out


def reads_locus():
    out = {"best_ref": NA, "aligned_reads": NA, "median_depth": NA, "breadth_5x": NA}
    try:
        text = open(stdout_path, errors="replace").read()
    except OSError:
        return out

    def grab(key):
        found = re.findall(r"metrics %s: ([0-9.eE+-]+) >=" % key, text)
        return float(found[-1]) if found else None

    reads, median, breadth = grab("min_reads"), grab("min_med"), grab("min_breadth")
    if reads is not None:
        out["aligned_reads"], out["median_depth"], out["breadth_5x"] = int(reads), f"{median:.1f}", f"{breadth:.3f}"
    thr = f"{aln}/thresholds.csv"
    if os.path.isfile(thr):
        with open(thr, newline="") as fh:
            rows = list(csv.reader(fh))
        if len(rows) > 1 and reads is not None:
            out["best_ref"] = rows[1][0]
    return out


if mode == "assembly":
    locus = assembly_locus()
else:
    locus = reads_locus()
locus_text = ";".join(f"{k}={v}" for k, v in locus.items())

# ---- the consensus protein ---------------------------------------------
recs = read_fasta(f"{aln}/{sid}_prot.fasta")
protein = recs[0][1].upper().rstrip("*") if recs else ""
listt_error = last_line(stderr_path) if listt_exit != "0" else ""

row = {
    "sample": sample, "mode": mode, "status": "QF", "call": "", "closest_ist": NA,
    "species": NA, "closest_variant": NA, "identity": NA, "threshold": NA,
    "exact_match": NA, "runner_up_ist": NA, "runner_up_identity": NA,
    "aa_length": len(protein), "reason": "", "locus_qc": locus_text,
    "listt_exit": listt_exit, "listt_error": listt_error or NA,
}


def quality_failure():
    """Why there is no usable ospA sequence, from what LISTT left behind."""
    if mode == "assembly":
        if locus["nt_len"] == 0:
            return "no_ospa_hit"
        if locus["nt_len"] < 813:
            return f"ospa_incomplete_{locus['nt_len']}nt"
        return "listt_failed"
    if locus["aligned_reads"] == NA:
        return "no_ospa_alignment"
    fails = []
    if locus["aligned_reads"] < min_cov["min_reads"]:
        fails.append(f"aligned_reads={locus['aligned_reads']}<{min_cov['min_reads']:g}")
    if float(locus["median_depth"]) < min_cov["min_med"]:
        fails.append(f"median_depth={locus['median_depth']}<{min_cov['min_med']:g}")
    if float(locus["breadth_5x"]) < min_cov["min_breadth"]:
        fails.append(f"breadth_5x={locus['breadth_5x']}<{min_cov['min_breadth']:g}")
    return "reads_qc_failed:" + ",".join(fails) if fails else "listt_failed"


identities = {}
if not protein:
    row["reason"] = quality_failure()
elif "*" in protein:
    row["reason"] = "premature_stop_codon"
else:
    exact = [v for v, s in ref_protein.items() if s == protein]
    if exact:
        closest, identity = exact[0], 100.0
        identities = {closest: 100.0}
    else:
        # LISTT writes every variant's identity (a fraction) before it looks
        # the winner up in its serotype table, so this survives its crash.
        path = f"{aln}/{sid}_alignments.csv"
        if os.path.isfile(path):
            with open(path, newline="") as fh:
                for r in csv.DictReader(fh):
                    identities[int(r["variant_id"])] = 100.0 * float(r["identity"])
        if not identities:
            row["reason"] = "no_identity_table"
        else:
            closest = None
            for v, pct in identities.items():  # LISTT's tie rule: first maximum wins
                if closest is None or pct > identities[closest]:
                    closest = v
            identity = identities[closest]

    if identities:
        info = variant_info.get(closest)
        ist = info[0] if info else NA
        row["closest_variant"] = closest
        row["identity"] = f"{identity:.2f}"
        row["exact_match"] = "yes" if exact else "no"
        row["species"] = info[1] if info else NA
        row["closest_ist"] = ist
        if len(identities) > 1:
            others = [(pct, v) for v, pct in identities.items()
                      if (variant_info.get(v) or (NA,))[0] != ist]
            if others:
                pct, v = max(others, key=lambda t: t[0])
                row["runner_up_ist"] = (variant_info.get(v) or (NA,))[0]
                row["runner_up_identity"] = f"{pct:.2f}"
        thr = thresholds.get(int(ist)) if ist.isdigit() else None
        if thr is not None:
            row["threshold"] = f"{thr:g}"

        if info is None:
            row["status"], row["reason"] = "PROVISIONAL", "variant_not_in_listt_table"
        elif not ist.isdigit():
            row["status"], row["reason"] = "PROVISIONAL", "variant_has_no_ist"
        elif int(ist) > max_trusted:
            row["status"], row["reason"] = "PROVISIONAL", f"ist_above_max_trusted_{max_trusted}"
        elif exact:
            row["status"], row["call"] = "CALL", f"IST{ist}"
        elif thr is None:
            row["status"], row["reason"] = "PROVISIONAL", "no_threshold_for_ist"
        elif identity + EPS >= thr:
            row["status"], row["call"] = "CALL", f"IST{ist}"
        else:
            row["status"], row["reason"] = "BELOW", "identity_below_ist_threshold"

if row["reason"] == "":
    row["reason"] = NA

# ---- cross-check against LISTT's own report ----------------------------
report = f"{stage}/aln/{sid}_report.csv"
if os.path.isfile(report) and row["closest_variant"] != NA:
    with open(report, newline="") as fh:
        rep = {r[0]: r[1] for r in csv.reader(fh) if len(r) >= 2}
    if rep.get("variant", "").isdigit() and int(rep["variant"]) != row["closest_variant"]:
        row["reason"] += f";differs_from_listt_report_variant_{rep['variant']}"

COLUMNS = list(row)
with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerow([row[c] for c in COLUMNS])
with open(out_fasta, "w") as fh:
    if protein:
        fh.write(f">{sample}|{mode}|variant{row['closest_variant']}|{row['status']}\n{protein}\n")
PY

        python3 classify.py "${SAMPLE}" "${SID}" '~{mode}' "${STAGE}" "${LISTT_EXIT}" '~{max_trusted_ist}' \
            listt_mode.tsv listt_protein.fasta "${WD}/listt.stdout" "${WD}/listt.stderr"

        cp "${WD}/listt.stdout" "${WD}/listt.stderr" "${STAGE}/aln/" 2>/dev/null || true
        tar -czf listt_raw.tar.gz -C "${STAGE}" aln
        echo "LISTT ${LISTT_VERSION:-unknown} ($(cut -c1-7 /opt/LISTT/COMMIT))" > listt_version.txt
    >>>

    output {
        File   mode_tsv       = "listt_mode.tsv"
        File   protein_fasta  = "listt_protein.fasta"
        File   raw_tar        = "listt_raw.tar.gz"
        String listt_version  = read_string("listt_version.txt")
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{disk} SSD"
        preemptible: preemptible
    }
}

task reconcile_ospA_calls {

    input {
        String  sample_name
        File    assembly_tsv
        File?   reads_tsv
        File    assembly_protein
        File?   reads_protein
        File    assembly_raw
        File?   reads_raw
        String  listt_version
        Int     max_trusted_ist
        String  docker = "aarvani1/listt:1.2.0@sha256:c016417e25c3e4b96d4fc98bedbf4d77e15608674e5e5320cdaff1df6a1a38b2"
    }

    parameter_meta {
        sample_name:      "Sample label"
        assembly_tsv:     "mode_tsv from the assembly run"
        reads_tsv:        "mode_tsv from the reads run, if one was made"
        assembly_protein: "protein_fasta from the assembly run"
        reads_protein:    "protein_fasta from the reads run"
        assembly_raw:     "raw_tar from the assembly run"
        reads_raw:        "raw_tar from the reads run"
        listt_version:    "LISTT version string"
        max_trusted_ist:  "Highest IST number reported as a call; recorded in the details"
        docker:           "Container image"
    }

    meta {
        description: "Combine the assembly and (optional) reads results into one call and one flag. Flags: OK = every supplied mode agrees on an IST; SM = one mode supports the call and no other mode contradicts it; UK = no mode gives a confident IST, so the closest is reported but not called; DS = modes disagree; QF = no usable ospA in any supplied mode. The call is IST<n>, Unknown (UK, QF) or Ambig (DS)."
    }

    command <<<
        set -euo pipefail

        cat > reconcile.py <<'PY'
import csv
import sys

sample, asm_path, reads_path, max_trusted, version = sys.argv[1:6]
NA = "NA"


def load(path):
    if not path:
        return None
    with open(path, newline="") as fh:
        return next(csv.DictReader(fh, delimiter="\t"))


def protein(path):
    if not path:
        return ""
    lines = [l.strip() for l in open(path) if l.strip() and not l.startswith(">")]
    return lines[0] if lines else ""


asm, reads = load(asm_path), load(reads_path or "")
modes = [(n, r) for n, r in (("assembly", asm), ("reads", reads)) if r is not None]
calls = [(n, r) for n, r in modes if r["status"] == "CALL"]
lows = [(n, r) for n, r in modes if r["status"] in ("BELOW", "PROVISIONAL")]
fails = [(n, r) for n, r in modes if r["status"] == "QF"]

asm_aa, reads_aa = protein(sys.argv[6]), protein(sys.argv[7] if len(sys.argv) > 7 else "")
same_protein = None
if asm_aa and reads_aa:
    same_protein = asm_aa == reads_aa

closest = {r["closest_ist"] for _, r in modes if r["closest_ist"] != NA}
notes = []

if modes and len(calls) == len(modes):
    if len({r["call"] for _, r in calls}) == 1:
        flag, call = ("OK" if len(modes) > 1 else "SM"), calls[0][1]["call"]
    else:
        flag, call = "DS", "Ambig"
elif calls:
    if all(r["closest_ist"] == calls[0][1]["closest_ist"] for _, r in lows):
        if lows:
            flag, call = "UK", "Unknown"
            notes.append("one mode is below threshold or provisional for the same IST")
        else:
            flag, call = "SM", calls[0][1]["call"]
            notes.append("the other mode failed QC: " + ",".join(r["reason"] for _, r in fails))
    else:
        flag, call = "DS", "Ambig"
elif lows:
    if len(closest) > 1:
        flag, call = "DS", "Ambig"
    else:
        flag, call = "UK", "Unknown"
        if same_protein:
            notes.append("both modes recovered the identical protein: a strong candidate for a new ospA allele"
                         " if the identity is good and the QC clean")
else:
    flag, call = "QF", "Unknown"

if same_protein is False and flag in ("OK", "UK"):
    notes.append("the two consensus proteins differ")

# species: from the modes behind the answer; "a|b" if they differ
src = calls if (flag in ("OK", "SM") and calls) else (lows or modes)
species = sorted({r["species"] for _, r in src if r["species"] != NA})
species_text = "|".join(species) if species else NA
if flag == "QF":
    species_text = NA

identity_text = ";".join(f"{n}={r['identity']}" for n, r in modes)


def ist_text(x):
    return f"IST{x}" if x.isdigit() else x


def describe(n, r):
    s = f"{n}={r['status']}"
    if r["status"] == "CALL":
        s += f" {r['call']} id={r['identity']} thr={r['threshold']}"
    elif r["status"] in ("BELOW", "PROVISIONAL"):
        s += f" closest={ist_text(r['closest_ist'])} variant={r['closest_variant']} id={r['identity']} thr={r['threshold']}"
        if r["runner_up_ist"] != NA:
            s += f" runner_up={ist_text(r['runner_up_ist'])}@{r['runner_up_identity']}"
        s += f" reason={r['reason']}"
    else:
        s += f" reason={r['reason']}"
    if r["listt_exit"] != "0":
        s += f" listt_exit={r['listt_exit']}({r['listt_error']})"
    return s


details = " | ".join(describe(n, r) for n, r in modes)
if notes:
    details += " | " + "; ".join(notes)
details += f" | max_trusted_ist={max_trusted} | {version}"

columns = list(asm)
final = {c: NA for c in columns}
final.update(sample=sample, mode="final", status=flag, call=call, species=species_text,
             identity=identity_text, closest_ist="|".join(sorted(closest)) or NA,
             reason=NA)
with open("ospA_results.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(columns)
    for _, r in modes:
        w.writerow([r[c] for c in columns])
    w.writerow([final[c] for c in columns])

for name, text in (("call", call), ("flag", flag), ("species", species_text),
                   ("identity", identity_text), ("details", details)):
    with open(f"{name}.txt", "w") as fh:
        fh.write(text + "\n")
PY

        python3 reconcile.py "$(cat ~{write_lines([sample_name])})" '~{assembly_tsv}' '~{default="" reads_tsv}' '~{max_trusted_ist}' '~{listt_version}' \
            '~{assembly_protein}' '~{default="" reads_protein}'

        cat '~{assembly_protein}' ~{default="" reads_protein} > ospA_proteins.fasta

        mkdir -p raw/assembly
        tar -xzf '~{assembly_raw}' -C raw/assembly
        if [ -n '~{default="" reads_raw}' ]; then
            mkdir -p raw/reads
            tar -xzf '~{default="" reads_raw}' -C raw/reads
        fi
        tar -czf ospA_raw.tar.gz -C raw .
    >>>

    output {
        String final_call     = read_string("call.txt")
        String final_flag     = read_string("flag.txt")
        String final_species  = read_string("species.txt")
        String final_identity = read_string("identity.txt")
        String final_details  = read_string("details.txt")
        File   results_tsv   = "ospA_results.tsv"
        File   protein_fasta = "ospA_proteins.fasta"
        File   raw_tar       = "ospA_raw.tar.gz"
    }

    runtime {
        docker:      docker
        cpu:         1
        memory:      "1 GB"
        disks:       "local-disk 5 SSD"
        preemptible: 1
    }
}
