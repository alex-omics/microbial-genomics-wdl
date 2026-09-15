#!/usr/bin/env bash
#
# Fixture tests for this repo's workflows.
#
# These run the REAL WDL tasks via miniwdl rather than a copy of their logic,
# so the tests cannot drift away from what the workflows actually execute.
#
# Methylation tasks: both bugs these cover were caught by hand and produced
# plausible-looking wrong numbers rather than errors:
#   - a bedtools column off-by-one that wrote the feature end coordinate into
#     the feature_strand column
#   - a missing percent-modified filter that let every evaluated base through,
#     turning the ortholog matrix into a proxy for gene length
#
# pyseer_gwas / fetch_reads_from_sra: pyseer_annotate_results is the one part
# of pyseer_gwas that isn't just a thin wrapper around pyseer itself, so it's
# the part most likely to produce a plausible-looking wrong answer rather than
# an error: a silently dropped unmatched variant, a p-value that rounds to
# 0.000000, or a comma inside a free-text annotation shifting every column
# after it. The pyseer_gwas section runs the whole workflow end to end against
# the real pyseer image, which is what actually catches wiring bugs between
# tasks - see the comment there. The fetch_reads_from_sra section runs
# fasterq_dump against real, tiny public accessions to prove paired-end and
# single-end layout detection both work against the actual tool, not a mocked
# one - see the comment there for why single-end coverage specifically
# matters.
#
# Requires: miniwdl, docker, python3. The fetch_reads_from_sra section also
# needs outbound network access (NCBI SRA + ENA).
# Usage: tests/run_unit_tests.sh

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIX="${REPO}/tests/fixtures/methylation"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

PASS=0
FAIL=0

check() {
    local label="$1" expected="$2" actual="$3"
    if [ "${expected}" = "${actual}" ]; then
        printf '  ok    %-52s %s\n' "${label}" "${actual}"
        PASS=$((PASS + 1))
    else
        printf '  FAIL  %-52s expected [%s] got [%s]\n' "${label}" "${expected}" "${actual}"
        FAIL=$((FAIL + 1))
    fi
}

command -v miniwdl >/dev/null || { echo "miniwdl not found"; exit 1; }
docker info >/dev/null 2>&1  || { echo "docker not running";  exit 1; }

# --------------------------------------------------------------------------
echo "annotate_methylation"
# Fixture exercises every filter: one clean site, one below the percent floor,
# one below the coverage floor, one at 15% with strong read support (below the
# default percent floor but recoverable by lowering it), one 4mC site, and one
# site in an upstream window.
miniwdl run "${REPO}/tasks/annotate_methylation.wdl" \
    bedmethyl="${FIX}/pileup.bedmethyl.bed" \
    sample_name=test \
    reference_gff="${FIX}/reference.gff3" \
    reference_fai="${FIX}/reference.fa.fai" \
    --dir "${WORK}/annotate" --verbose > "${WORK}/annotate.log" 2>&1 || {
        echo "  task failed; see ${WORK}/annotate.log"; sed -n '$p' "${WORK}/annotate.log"; exit 1; }

OUT="$(find "${WORK}/annotate" -name 'test_methylation_annotated.tsv' | head -1)"

check "rows kept after all three floors" "3" "$(awk 'NR>1' "${OUT}" | wc -l | tr -d ' ')"
check "genic sites"                      "2" "$(awk -F'\t' 'NR>1 && $14=="genic"    {n++} END{print n+0}' "${OUT}")"
check "upstream sites"                   "1" "$(awk -F'\t' 'NR>1 && $14=="upstream" {n++} END{print n+0}' "${OUT}")"
# The off-by-one regression: these two must not hold coordinates.
check "feature_strand column is a strand" "-" "$(awk -F'\t' 'NR==2{print $11}' "${OUT}")"
check "region_length is gene length"    "500" "$(awk -F'\t' 'NR==2{print $12}' "${OUT}")"
check "product survives its comma"  "L-aspartate oxidase, FAD-binding" \
      "$(awk -F'\t' 'NR>1 && $9=="TEST_00002" && $14=="genic" {print $13; exit}' "${OUT}")"
check "6mA code preserved"               "a" "$(awk -F'\t' 'NR==2{print $6}' "${OUT}")"

# Lowering the percent floor must recover the well-supported 15% site.
miniwdl run "${REPO}/tasks/annotate_methylation.wdl" \
    bedmethyl="${FIX}/pileup.bedmethyl.bed" \
    sample_name=lowthresh \
    reference_gff="${FIX}/reference.gff3" \
    reference_fai="${FIX}/reference.fa.fai" \
    min_percent=10 \
    --dir "${WORK}/lowthresh" > "${WORK}/lowthresh.log" 2>&1 || {
        echo "  task failed; see ${WORK}/lowthresh.log"; exit 1; }
OUT2="$(find "${WORK}/lowthresh" -name 'lowthresh_methylation_annotated.tsv' | head -1)"
check "min_percent=10 recovers partial site" "4" "$(awk 'NR>1' "${OUT2}" | wc -l | tr -d ' ')"

# --------------------------------------------------------------------------
echo "methylation_orthologs"
miniwdl run "${REPO}/tasks/methylation_orthologs.wdl" \
    gene_presence_absence="${FIX}/gene_presence_absence.csv" \
    annotated_tables="${FIX}/isoA_methylation_annotated.tsv" \
    annotated_tables="${FIX}/isoB_methylation_annotated.tsv" \
    sample_names=isoA sample_names=isoB \
    --dir "${WORK}/ortho" > "${WORK}/ortho.log" 2>&1 || {
        echo "  task failed; see ${WORK}/ortho.log"; sed -n '$p' "${WORK}/ortho.log"; exit 1; }

MAT="$(find "${WORK}/ortho" -name 'methylation_by_ortholog_matrix.tsv' | head -1)"
LONG="$(find "${WORK}/ortho" -name 'methylation_by_ortholog_long.tsv'   | head -1)"

get() { awk -F'\t' -v g="$1" -v c="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $1==g{print $(h[c])}' "$3"; }

# The distinction the whole matrix rests on: absent gene vs present-unmethylated.
check "absent gene is NA, not 0"    "NA" "$(get group_B isoA "${MAT}")"
check "present but unmethylated is 0" "0" "$(get group_A isoB "${MAT}")"
check "methylated gene counts sites" "2"  "$(get group_A isoA "${MAT}")"
check "quoted comma kept in annotation" "efflux pump, RND family" \
      "$(get group_A annotation "${MAT}")"
check "paralogues collapse into one group" "1" "$(get group_C isoA "${MAT}")"

DENS="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
        $1=="group_A" && $5=="isoA" {print $(h["genic_sites_per_kb"])}' "${LONG}")"
check "density normalises by length" "2.000" "${DENS}"

# --------------------------------------------------------------------------
echo "rebase_blastp (real BLASTP + regression guard on the evalue bug)"
# evalue is a String, not a Float, on purpose: WDL renders a Float this small
# (1e-25) as the fixed-decimal text "0.000000", which blastp then rejects as
# a non-positive e-value cutoff. This ran cleanly against a live default
# before that fix was in place -- miniwdl check cannot catch it, since it is
# a runtime string-interpolation behaviour, not a syntax error. Only a real
# execution with the actual default value catches a regression here.
REBASE_FIX="${REPO}/tests/fixtures/rebase"
miniwdl run "${REPO}/tasks/rebase_mtase_search.wdl" --task rebase_blastp \
    faa="${REBASE_FIX}/query.faa" \
    sample_name=test \
    rebase_goldset_fasta="${REBASE_FIX}/goldset.faa" \
    --dir "${WORK}/blastp" --verbose > "${WORK}/blastp.log" 2>&1 || {
        echo "  task failed; see ${WORK}/blastp.log"; sed -n '$p' "${WORK}/blastp.log"; exit 1; }
HITS="$(find "${WORK}/blastp" -path '*/out/*' -name '*_blast_hits.tsv' | head -1)"
check "blastp finds the planted 100% identity hit" "ISO1_00010	M.AacDam" \
    "$(cut -f1,2 "${HITS}")"

# --------------------------------------------------------------------------
echo "rebase_mtase_search (join logic)"
# BLASTP hit with a characterised REBASE motif must carry it through; a hit
# whose REBASE entry has no known recognition sequence (Motif="?") must be
# kept as NA rather than silently dropped, since "homology confirmed, motif
# unknown" is still real information about the isolate.
MOTIF_FIX="${REPO}/tests/fixtures/motif_landscape"
# rebase_blastp itself was already exercised above through real docker. This
# section is join.py specifically -- extracted and run directly, the way the
# ortholog join is tested, since the branching (characterised vs. unknown
# motif) is the part actually worth a regression check.
sed -n '/^        cat > join.py/,/^        PY$/p' "${REPO}/tasks/rebase_mtase_search.wdl" | sed '1d;$d' | sed 's/^        //' > "${WORK}/rebase_join.py"
python3 "${WORK}/rebase_join.py" \
    "${MOTIF_FIX}/blast_hits.tsv" "${MOTIF_FIX}/rebase_motif_data.tsv" iso1 "${WORK}/rebase_out.tsv"
check "characterised motif carried through" "GATC" \
    "$(awk -F'\t' '$2=="ISO1_00010"{print $7}' "${WORK}/rebase_out.tsv")"
check "unknown-motif hit kept as NA, not dropped" "NA" \
    "$(awk -F'\t' '$2=="ISO1_00099"{print $7}' "${WORK}/rebase_out.tsv")"
# REBASE entries recognising more than one degenerate sequence for the same
# enzyme (e.g. M.BceSV: "GGCC,GCNGC,CCGG,GGNCC") pack them into one
# comma-joined Motif field. Passing that string through unsplit produced a
# non-IUPAC "motif" that crashed motif_landscape's regex translation on the
# literal comma -- caught only by running the real 14-isolate panel.
check "comma-joined REBASE motif splits into separate rows" "2" \
    "$(awk -F'\t' '$2=="ISO1_00040"' "${WORK}/rebase_out.tsv" | wc -l | tr -d ' ')"
check "each split motif is valid IUPAC, no literal comma" "GCNGC
GGCC" \
    "$(awk -F'\t' '$2=="ISO1_00040"{print $7}' "${WORK}/rebase_out.tsv" | sort)"

# --------------------------------------------------------------------------
echo "build_motif_list"
miniwdl run "${REPO}/tasks/motif_landscape.wdl" --task build_motif_list \
    find_motifs_tsv="${MOTIF_FIX}/find_motifs.tsv" \
    rebase_mtases_tsv="${MOTIF_FIX}/rebase_mtases.tsv" \
    --dir "${WORK}/bml" --verbose > "${WORK}/bml.log" 2>&1 || {
        echo "  task failed; see ${WORK}/bml.log"; exit 1; }
LIST="$(find "${WORK}/bml" -path '*/out/*' -name 'motif_list.tsv' | head -1)"
check "de novo and REBASE agreement dedupes to one row" "1" \
    "$(awk -F'\t' 'NR>1 && $1=="GATC" && $2=="a"' "${LIST}" | wc -l | tr -d ' ')"
check "REBASE-only motif also present" "1" \
    "$(awk -F'\t' 'NR>1 && $1=="CCWGG" && $2=="m"' "${LIST}" | wc -l | tr -d ' ')"
check "row count matches expected dedup" "2" "$(awk 'NR>1' "${LIST}" | wc -l | tr -d ' ')"

# --------------------------------------------------------------------------
echo "motif_landscape (tier 1 + tier 2, against ground-truth planted signal)"
# iso1.bedmethyl.bed plants 6mA at every genomic GATC copy but makes a third
# of those copies unmethylated (phase-variation-like); iso2 plants the same
# motif uniformly high everywhere. Expected numbers were hand-verified
# against the planted values before being fixed here as regression checks.
miniwdl run "${REPO}/tasks/motif_landscape.wdl" --task motif_landscape \
    bedmethyl="${MOTIF_FIX}/iso1.bedmethyl.bed" \
    motif_list="${MOTIF_FIX}/motif_list.tsv" \
    sample_name=iso1 \
    reference_fasta="${MOTIF_FIX}/reference.fa" \
    --dir "${WORK}/ml1" > "${WORK}/ml1.log" 2>&1 || {
        echo "  task failed; see ${WORK}/ml1.log"; exit 1; }
ML1="$(find "${WORK}/ml1" -path '*/out/*' -name '*_motif_landscape.tsv' | head -1)"
check "tier1: heterogeneous isolate enrichment" "65.25" "$(awk -F'\t' 'NR==2{print $6}' "${ML1}")"
check "tier1: background rate near planted 3%"   "3.00"  "$(awk -F'\t' 'NR==2{print $7}' "${ML1}")"
check "tier2: heterogeneous isolate CV is high"  "0.665" "$(awk -F'\t' 'NR==2{print $10}' "${ML1}")"
check "tier2: recovers exact count of planted low copies" "153" "$(awk -F'\t' 'NR==2{print $12}' "${ML1}")"

miniwdl run "${REPO}/tasks/motif_landscape.wdl" --task motif_landscape \
    bedmethyl="${MOTIF_FIX}/iso2.bedmethyl.bed" \
    motif_list="${MOTIF_FIX}/motif_list.tsv" \
    sample_name=iso2 \
    reference_fasta="${MOTIF_FIX}/reference.fa" \
    --dir "${WORK}/ml2" > "${WORK}/ml2.log" 2>&1 || {
        echo "  task failed; see ${WORK}/ml2.log"; exit 1; }
ML2="$(find "${WORK}/ml2" -path '*/out/*' -name '*_motif_landscape.tsv' | head -1)"
check "tier1: uniform isolate enrichment"        "97.00" "$(awk -F'\t' 'NR==2{print $6}' "${ML2}")"
check "tier2: uniform isolate CV near zero (housekeeping null)" "0.006" "$(awk -F'\t' 'NR==2{print $10}' "${ML2}")"
check "tier2: uniform isolate has no low-methylation copies"    "0"     "$(awk -F'\t' 'NR==2{print $12}' "${ML2}")"

# --------------------------------------------------------------------------
echo "motif_landscape_summary (tier 3, run directly -- see note below)"
# Run via the extracted script rather than through miniwdl/docker: chaining
# one miniwdl run's file outputs into a second run's inputs hits a sandbox-
# specific virtiofs mount restriction in this environment that does not
# reflect anything about the task itself (confirmed by successfully running
# motif_landscape standalone above, against plain fixture files, in the same
# environment). The checked-in landscape fixtures are themselves the
# validated output of the miniwdl run above, so this still exercises real,
# previously-miniwdl-verified numbers through the actual summarise.py logic.
sed -n '/^        cat > summarise.py/,/^        PY$/p' "${REPO}/tasks/motif_landscape_summary.wdl" | sed '1d;$d' | sed 's/^        //' > "${WORK}/summarise.py"
python3 "${WORK}/summarise.py" \
    "${MOTIF_FIX}/iso1_landscape.tsv,${MOTIF_FIX}/iso2_landscape.tsv" \
    "${MOTIF_FIX}/ortholog_long.tsv" \
    "mex,opr" \
    "${WORK}/by_motif.tsv" "${WORK}/by_gene.tsv"
check "motif-level cross-isolate mean" "81.12" "$(awk -F'\t' 'NR==2{print $5}' "${WORK}/by_motif.tsv")"
check "motif-level cross-isolate CV"   "0.196" "$(awk -F'\t' 'NR==2{print $8}' "${WORK}/by_motif.tsv")"
check "gene-level: highlighted gene flagged" "yes" \
    "$(awk -F'\t' '$1=="group_A"{print $8}' "${WORK}/by_gene.tsv")"
check "gene-level: non-AMR gene not flagged" "" \
    "$(awk -F'\t' '$1=="group_B"{print $8}' "${WORK}/by_gene.tsv")"

python3 "${WORK}/summarise.py" \
    "${MOTIF_FIX}/iso1_landscape.tsv,${MOTIF_FIX}/iso2_landscape.tsv" \
    "${MOTIF_FIX}/ortholog_long_conflict.tsv" \
    "mex,opr" \
    /dev/null "${WORK}/by_gene_conflict.tsv"
check "highlighted gene sorts first even with lower CV" "group_A" \
    "$(awk -F'\t' 'NR==2{print $1}' "${WORK}/by_gene_conflict.tsv")"

# --------------------------------------------------------------------------
echo "align_modbam length filter"
# The filter that would have caught the false "different isolate" alarm this
# session actually hit: a raw modBAM's short-fragment tail (~50-90bp junk,
# empirically) is mechanically incapable of a confident minimap2 placement,
# so leaving it in the denominator made percent_mapped measure read-length
# composition, not whether the modbam and assembly correspond. Extracted
# directly from the task rather than copied, so it can't drift from what
# actually runs. Fixture covers the boundary: 499bp dropped, exactly 500bp
# (the default) kept, plus a trivially short and a comfortably long read.
ALIGN_FIX="${REPO}/tests/fixtures/align_modbam"
sed -n '/^        awk -v minlen=~{min_read_length}/,/^        .*'"'"' reads.fastq/p' "${REPO}/tasks/align_modbam.wdl" \
    | sed '1s/.*/awk -v minlen=500 \x27/; $s/.*/\x27/' > "${WORK}/lenfilter.sh"
bash "${WORK}/lenfilter.sh" < "${ALIGN_FIX}/minireads.fastq" > "${WORK}/filtered.fastq"

check "boundary: exactly-500bp read is kept" "1" \
    "$(grep -c '^@exactly500$' "${WORK}/filtered.fastq")"
check "boundary: 499bp read is dropped" "0" \
    "$(grep -c '^@short2$' "${WORK}/filtered.fastq")"
check "trivially short read is dropped" "0" \
    "$(grep -c '^@short1$' "${WORK}/filtered.fastq")"
check "long read is kept" "1" \
    "$(grep -c '^@long1$' "${WORK}/filtered.fastq")"
check "output FASTQ structure is not corrupted (2 records = 8 lines)" "8" \
    "$(wc -l < "${WORK}/filtered.fastq" | tr -d ' ')"

# --------------------------------------------------------------------------
echo "pyseer_annotate_results"
# Fixture exercises: a variant with a comma inside its annotation text (must
# survive real CSV parsing, not comma-splitting), a variant present in the
# annotation table with an empty gene-name field (must come back as "",
# distinct from "no match"), a variant absent from the annotation table
# entirely (must come back labelled "no match", not dropped), a p-value small
# enough to require scientific notation rather than rounding to 0.000000, and
# a non-numeric notes/qc_flags value that must pass through unchanged.
PYSEER_FIX="${REPO}/tests/fixtures/pyseer_gwas"
miniwdl run "${REPO}/tasks/pyseer.wdl" --task pyseer_annotate_results \
    pyseer_results="${PYSEER_FIX}/pyseer_gene_results.tsv" \
    annotation_table="${PYSEER_FIX}/gene_presence_absence.csv" \
    id_column="Gene" \
    name_column="Non-unique Gene name" \
    annotation_column="Annotation" \
    basename=test \
    docker="python:3.11-slim" \
    --dir "${WORK}/annotate_pyseer" --verbose > "${WORK}/annotate_pyseer.log" 2>&1 || {
        echo "  task failed; see ${WORK}/annotate_pyseer.log"; sed -n '$p' "${WORK}/annotate_pyseer.log"; exit 1; }

ANNOTATED="$(find "${WORK}/annotate_pyseer" -name 'test_annotated.tsv' | head -1)"
READABLE="$(find "${WORK}/annotate_pyseer" -name 'test_readable.tsv' | head -1)"

check "annotated: pyseer's own columns untouched" \
    "variant	gene_name	annotation	af	filter-pvalue	lrt-pvalue	beta	beta-std-err	variant_h2	notes" \
    "$(head -1 "${ANNOTATED}")"
check "annotated: comma-in-annotation survives real CSV parsing" \
    "efflux pump, RND family" \
    "$(awk -F'\t' '$1=="group_0001" {print $3}' "${ANNOTATED}")"
check "annotated: matched-but-empty gene name is '', not 'no match'" \
    "" \
    "$(awk -F'\t' '$1=="group_0002" {print $2}' "${ANNOTATED}")"
check "annotated: variant absent from annotation table is labelled" \
    "no match	no match" \
    "$(awk -F'\t' '$1=="group_0003" {print $2"\t"$3}' "${ANNOTATED}")"

check "readable: renamed header" \
    "gene_family_id	gene_name	annotation	allele_frequency	prefilter_pvalue	lrt_pvalue	effect_size	effect_size_stderr	variance_explained	qc_flags" \
    "$(head -1 "${READABLE}")"
check "readable: p-value below 1e-4 stays scientific, not 0.000000" \
    "1.200e-15" \
    "$(awk -F'\t' '$1=="group_0001" {print $6}' "${READABLE}")"
check "readable: p-value above 1e-4 is fixed-decimal" \
    "0.050000" \
    "$(awk -F'\t' '$1=="group_0002" {print $6}' "${READABLE}")"
check "readable: af is fixed at 6 decimal places" \
    "0.500000" \
    "$(awk -F'\t' '$1=="group_0001" {print $4}' "${READABLE}")"
check "readable: notes passes through as qc_flags unchanged" \
    "bad-chisq" \
    "$(awk -F'\t' '$1=="group_0002" {print $10}' "${READABLE}")"

# --------------------------------------------------------------------------
echo "pyseer_annotate_results (no annotation_table)"
# Standalone path: skip the join cleanly, still produce the renamed/
# reformatted readable file, and do not produce an annotated file at all.
miniwdl run "${REPO}/tasks/pyseer.wdl" --task pyseer_annotate_results \
    --none annotation_table \
    pyseer_results="${PYSEER_FIX}/pyseer_gene_results.tsv" \
    basename=test_noannot \
    docker="python:3.11-slim" \
    --dir "${WORK}/noannot" --verbose > "${WORK}/noannot.log" 2>&1 || {
        echo "  task failed; see ${WORK}/noannot.log"; sed -n '$p' "${WORK}/noannot.log"; exit 1; }

check "no annotation_table: annotated output is not produced" \
    "0" \
    "$(find "${WORK}/noannot" -name 'test_noannot_annotated.tsv' | wc -l | tr -d ' ')"
NOANNOT_READABLE="$(find "${WORK}/noannot" -name 'test_noannot_readable.tsv' | head -1)"
check "no annotation_table: readable header has no gene_name/annotation columns" \
    "gene_family_id	allele_frequency	prefilter_pvalue	lrt_pvalue	effect_size	effect_size_stderr	variance_explained	qc_flags" \
    "$(head -1 "${NOANNOT_READABLE}")"

# --------------------------------------------------------------------------
echo "pyseer_gwas (end-to-end smoke test, real pyseer)"
# Runs the actual workflow against a tiny synthetic 8-sample panel through the
# real pyseer image, exercising kinship, association, lineage effects,
# per-covariate scanning, and the annotation join together. This is the check
# that catches wiring bugs a unit test on a single task cannot: it was this
# exact test that caught pyseer's --lineage silently requiring a distance
# matrix (not just lineage_clusters) - the workflow ran clean up to this
# point and failed only inside the container, with no WDL-level type error to
# flag it. Not asserting on pyseer's own statistics (those are its job, not
# ours) - only on the shape of what came back.
#
# miniwdl's -i does not tolerate the "__comment*" keys used elsewhere in this
# repo's input JSONs for human-readable notes, so strip them first.
python3 -c "
import json
d = json.load(open('${REPO}/tests/pyseer_gwas_smoke_inputs.json'))
json.dump({k: v for k, v in d.items() if not k.startswith('__')}, open('${WORK}/smoke_inputs.json', 'w'))
"
miniwdl run "${REPO}/workflows/pyseer_gwas/pyseer_gwas.wdl" \
    -i "${WORK}/smoke_inputs.json" \
    --dir "${WORK}/smoke" --verbose > "${WORK}/smoke.log" 2>&1 || {
        echo "  workflow failed; see ${WORK}/smoke.log"; sed -n '$p' "${WORK}/smoke.log"; exit 1; }

SMOKE_GENE_RESULTS="$(find "${WORK}/smoke" -path '*call-pyseer_association*' -name 'pyseer_gene_results.tsv' | head -1)"
SMOKE_LINEAGE="$(find "${WORK}/smoke" -path '*call-pyseer_lineage_effects*' -name 'lineage_effects.txt' | head -1)"
SMOKE_READABLE="$(find "${WORK}/smoke" -path '*call-annotate_all*' -name 'pyseer_gene_results_readable.tsv' | head -1)"
SMOKE_COVSCAN_ANNOTATED="$(find "${WORK}/smoke" -path '*call-annotate_covariate_scan*' -name 'pyseer_covariate_scan_combined_annotated.tsv' | head -1)"

check "gene_results: invariant core gene filtered by max_af" \
    "0" \
    "$(awk -F'\t' '$1=="group_0004"' "${SMOKE_GENE_RESULTS}" | wc -l | tr -d ' ')"
check "lineage_effects: one row per BAPS cluster (2)" \
    "2" \
    "$(($(wc -l < "${SMOKE_LINEAGE}" | tr -d ' ') - 1))"
check "covariate_scan_results: one file per combination (RST, BAPS, RST+OspC)" \
    "3" \
    "$(find "${WORK}/smoke" -path '*call-pyseer_association*out/covariate_scan_results*' -name '*.tsv' | wc -l | tr -d ' ')"
check "covariate_scan_results: multi-column combination ran as one joint pyseer call" \
    "1" \
    "$(find "${WORK}/smoke" -path '*call-pyseer_association*out/covariate_scan_results*' -name 'RST+OspC.tsv' | wc -l | tr -d ' ')"
check "readable: variant absent from smoke_annotation.csv is 'no match'" \
    "no match" \
    "$(awk -F'\t' '$1=="group_0005" {print $2}' "${SMOKE_READABLE}")"
check "covariate_scan_combined_annotated: gene_name/annotation joined without disturbing the covariates column" \
    "variant	gene_name	annotation	covariates	af	filter-pvalue	lrt-pvalue	beta	beta-std-err	variant_h2	notes" \
    "$(head -1 "${SMOKE_COVSCAN_ANNOTATED}")"

# --------------------------------------------------------------------------
echo "fetch_reads_from_sra (fasterq_dump, real sra-tools, paired-end)"
# DRR727328 is a real, tiny (~25 KB gzipped) paired-end Illumina WGS run for
# Lacticaseibacillus rhamnosus GG. Confirms the common case still works
# exactly as before: both _1/_2 produced, layout detected from the actual
# files (not from SRA/ENA metadata, which can be stale), platform looked up.
miniwdl run "${REPO}/tasks/sra_tools.wdl" --task fasterq_dump \
    accession=DRR727328 \
    disk_gb=20 \
    --dir "${WORK}/sra_pe" --verbose > "${WORK}/sra_pe.log" 2>&1 || {
        echo "  task failed; see ${WORK}/sra_pe.log"; sed -n '$p' "${WORK}/sra_pe.log"; exit 1; }

check "paired: layout" "paired" \
    "$(find "${WORK}/sra_pe" -name LAYOUT -exec cat {} \;)"
check "paired: read1 produced and non-empty" "1" \
    "$(find "${WORK}/sra_pe" -path '*out/read1*' -name '*.fastq.gz' -size +0 | wc -l | tr -d ' ')"
check "paired: read2 produced and non-empty" "1" \
    "$(find "${WORK}/sra_pe" -path '*out/read2*' -name '*.fastq.gz' -size +0 | wc -l | tr -d ' ')"
check "paired: platform" "ILLUMINA" \
    "$(find "${WORK}/sra_pe" -name PLATFORM -exec cat {} \;)"

# --------------------------------------------------------------------------
echo "fetch_reads_from_sra (fasterq_dump, real sra-tools, single-end)"
# DRR572312 is a real, tiny (~19 KB) single-end ONT run for E. coli. This is
# the regression test for the bug that motivated porting this workflow in:
# the original script hardcoded gzip on _1/_2, which does not exist for a
# single-end run, and would have failed outright here. read2 must NOT be
# produced (not even empty) - it's an unset optional output, not a file.
miniwdl run "${REPO}/tasks/sra_tools.wdl" --task fasterq_dump \
    accession=DRR572312 \
    disk_gb=20 \
    --dir "${WORK}/sra_se" --verbose > "${WORK}/sra_se.log" 2>&1 || {
        echo "  task failed; see ${WORK}/sra_se.log"; sed -n '$p' "${WORK}/sra_se.log"; exit 1; }

check "single: layout" "single" \
    "$(find "${WORK}/sra_se" -name LAYOUT -exec cat {} \;)"
check "single: read1 produced and non-empty" "1" \
    "$(find "${WORK}/sra_se" -path '*out/read1*' -name '*.fastq.gz' -size +0 | wc -l | tr -d ' ')"
check "single: read2 is not produced" "0" \
    "$(find "${WORK}/sra_se" -path '*out/read2*' -name '*.fastq.gz' | wc -l | tr -d ' ')"
check "single: platform" "OXFORD_NANOPORE" \
    "$(find "${WORK}/sra_se" -name PLATFORM -exec cat {} \;)"

# --------------------------------------------------------------------------
echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
