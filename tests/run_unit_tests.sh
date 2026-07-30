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
echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
