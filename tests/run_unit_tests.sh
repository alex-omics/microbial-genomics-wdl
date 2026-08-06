#!/usr/bin/env bash
#
# Fixture tests for the methylation tasks.
#
# These run the REAL WDL tasks via miniwdl rather than a copy of their logic,
# so the test cannot drift away from what the workflow actually executes.
#
# Both bugs these cover were caught by hand and produced plausible-looking
# wrong numbers rather than errors:
#   - a bedtools column off-by-one that wrote the feature end coordinate into
#     the feature_strand column
#   - a missing percent-modified filter that let every evaluated base through,
#     turning the ortholog matrix into a proxy for gene length
#
# Requires: miniwdl, docker.
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
echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
