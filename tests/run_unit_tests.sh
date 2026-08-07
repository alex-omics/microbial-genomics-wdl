#!/usr/bin/env bash
#
# Fixture tests for the pyseer_gwas tasks.
#
# These run the REAL WDL task via miniwdl rather than a copy of its logic, so
# the test cannot drift away from what the workflow actually executes.
#
# pyseer_annotate_results is the one part of pyseer_gwas that isn't just a
# thin wrapper around pyseer itself, so it's the part most likely to produce
# a plausible-looking wrong answer rather than an error: a silently dropped
# unmatched variant, a p-value that rounds to 0.000000, or a comma inside a
# free-text annotation shifting every column after it. The final section runs
# the whole workflow end to end against the real pyseer image, which is what
# actually catches wiring bugs between tasks - see the comment there.
#
# Requires: miniwdl, docker, python3.
# Usage: tests/run_unit_tests.sh

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIX="${REPO}/tests/fixtures/pyseer_gwas"
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
echo "pyseer_annotate_results"
# Fixture exercises: a variant with a comma inside its annotation text (must
# survive real CSV parsing, not comma-splitting), a variant present in the
# annotation table with an empty gene-name field (must come back as "",
# distinct from "no match"), a variant absent from the annotation table
# entirely (must come back labelled "no match", not dropped), a p-value small
# enough to require scientific notation rather than rounding to 0.000000, and
# a non-numeric notes/qc_flags value that must pass through unchanged.
miniwdl run "${REPO}/tasks/pyseer.wdl" --task pyseer_annotate_results \
    pyseer_results="${FIX}/pyseer_gene_results.tsv" \
    annotation_table="${FIX}/gene_presence_absence.csv" \
    id_column="Gene" \
    name_column="Non-unique Gene name" \
    annotation_column="Annotation" \
    basename=test \
    docker="python:3.11-slim" \
    --dir "${WORK}/annotate" --verbose > "${WORK}/annotate.log" 2>&1 || {
        echo "  task failed; see ${WORK}/annotate.log"; sed -n '$p' "${WORK}/annotate.log"; exit 1; }

ANNOTATED="$(find "${WORK}/annotate" -name 'test_annotated.tsv' | head -1)"
READABLE="$(find "${WORK}/annotate" -name 'test_readable.tsv' | head -1)"

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
    pyseer_results="${FIX}/pyseer_gene_results.tsv" \
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
echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
