# pyseer_gwas

Standalone bacterial GWAS with [pyseer](https://github.com/mgalardini/pyseer), independently
callable in Terra. Runs pyseer's LMM (likelihood-ratio test) over a presence/absence Rtab,
correcting for population structure with a phylogeny-derived kinship matrix.

`presence_absence_rtab` is deliberately generic: pyseer's own docs describe the Rtab format
it accepts via `--pres` as usable "flexibly to represent variants from other sources" — a
`block_id x sample` 0/1 matrix, nothing gene-specific about it. A Panaroo
`gene_presence_absence.Rtab` and a pangenome-network module Rtab (module × isolate) both
plug in unchanged.

Deliberately separate from `multimodal_gwas` (`feature/multimodal-gwas`), which
re-implements the much larger microGWAS Snakemake DAG (unitigs, structural variants,
panfeed, whole-genome elastic net) as one component of a bigger pipeline with its own
upstream dependencies (mash distances, MLST lineages, a matched sample intersection).
This workflow runs pyseer on its own, with no such prerequisites. The two pipelines pin the
same pyseer version but use separate Docker images (`aarvani1/pyseer` vs.
`aarvani1/microgwas-pyseer`) and are expected to drift independently.

## What it does

| Step | Tool | Output |
| ---- | ---- | ------ |
| Kinship + distance matrix | `phylogeny_distance.py`, with and without `--lmm` | `kinship_matrix`, `distance_matrix` |
| Association | `pyseer --lmm` | `gene_results`, `gene_significant` (pattern-counted Bonferroni via `count_patterns.py`) |
| Lineage effects *(optional)* | `pyseer --lineage` | `lineage_effects` |
| Covariate combination scan *(optional)* | `pyseer --use-covariates`, once per combination | `covariate_scan_results`, `covariate_scan_combined` |
| SNP association *(optional)* | `pyseer --vcf` | `snp_results` |
| Name/annotation join *(optional)* | pandas left join, never touches pyseer's own files | `gene_results_annotated`, `gene_results_readable` (+ the same pair for `gene_significant`) |

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `phylogeny_newick` | `File` | **Required.** Midpoint-rooted core-genome phylogeny |
| `presence_absence_rtab` | `File` | **Required.** Gene or module presence/absence matrix |
| `phenotype_tsv` | `File` | **Required.** `sample_id\tphenotype_value` |
| `variant_vcf` | `File?` | Optional additional SNP-based association |
| `covariates_file` | `File?` | `sample_id` + one named column per covariate (e.g. RST, OspC, MLST, BAPS) |
| `use_covariates` | `String?` | pyseer's raw `--use-covariates` value, applied jointly to the main association |
| `covariate_combinations` | `Array[String]` | Groups of `covariates_file` column names to test together, one pyseer run per entry — see below |
| `run_lineage_effects` | `Boolean` | Adds `--lineage` to the main association (default `false`). Needs a distance matrix, which the workflow derives from `phylogeny_newick` automatically — nothing extra to pass in. pyseer refuses `--lineage` without one even when `lineage_clusters` is set; this is not just the fallback for MDS-derived lineages |
| `lineage_clusters` | `File?` | `sample_id\tcluster_id`, e.g. a BAPS export, for `--lineage-clusters` |
| `annotation_table` | `File?` | Panaroo/Roary `gene_presence_absence.csv`, or an equivalent module-annotation table |

### Covariates: what actually happens, verified against the real image

This section exists because every claim in it was checked by actually running `pyseer`
against `tests/fixtures/pyseer_gwas/smoke_covariates.tsv` (`sample_id, RST, OspC, MLST,
BAPS`) through the built image, not assumed from the docs. Three things are easy to get
wrong here, and getting them wrong doesn't error — it silently gives you a result that
looks fine.

**1. Supplying `covariates_file` does nothing by itself.** pyseer's own `--help` says it
plainly: `--use-covariates` defaults to *"load covariates but don't use them."* Confirmed
directly — `--covariates smoke_covariates.tsv` with no `--use-covariates` produces
output byte-for-byte identical to no covariates at all. You must also set either
`use_covariates` (raw pyseer syntax, applied jointly to the main run) or
`covariate_combinations` (below) — pointing at the file is not enough on its own.

**2. RST, OspC, MLST, and BAPS are categorical, not numeric, even though some are stored
as integers.** BAPS and MLST look like numbers in the TSV, but they're group labels —
BAPS cluster 3 isn't "more" than cluster 1, and there's no meaningful arithmetic between
ST types. Never append pyseer's `q` (quantitative) suffix to them. This only visibly
matters once a variable has 3+ levels: testing BAPS (2 levels in the fixture) as
categorical (`--use-covariates 5`) vs. quantitative (`--use-covariates 5q`) gave
*identical* results, purely because a 2-level categorical variable and a continuous one
are mathematically the same regressor. That equivalence breaks the moment a variable has
3+ levels — real BAPS clusters, and almost certainly MLST, both will. `covariate_combinations`
below always treats every column as categorical; there is no `q` option for it, and that's
deliberate. `use_covariates` is the only place a `q` could even be typed, and RST/OspC/
MLST/BAPS should never get one there either.

**3. Loading several categorical covariates jointly is a real degrees-of-freedom risk, not
a style preference.** Each categorical covariate is one-hot encoded — a k-level column
costs roughly k−1 fitted parameters, on top of the intercept and the variant term, all
drawn from your sample count. Measured on the 8-sample smoke fixture:

| Covariates in one joint run | Real variants that survived (of 4) |
| --- | --- |
| none | 4 |
| BAPS alone | 4 |
| OspC alone | 4 |
| RST alone | 4 |
| MLST alone (4 levels / 8 samples) | 2 — already halved |
| RST + BAPS | 4 |
| RST + OspC + BAPS (MLST excluded) | 4 |
| RST + MLST | 2 — MLST costs power even paired with just one other column |
| RST + OspC + MLST + BAPS (all four) | **0** — every real variant failed `lrt-filtering-failed` |

MLST is the one to watch: it's the finest-grained of the four (most distinct levels), so
it is the first to exhaust degrees of freedom, alone or in combination. This is a function
of level-count vs. sample-count, not a fixed rule — a much larger panel tolerates more —
but the mechanism doesn't go away, and it fails silently (a suspiciously empty results
table, not an error) unless you're watching for it.

### `covariate_combinations`: deliberate groups, not an automatic search

`covariate_combinations` runs one pyseer call per entry, each using exactly that group of
`covariates_file` columns as the sole covariates for that run — a bare name is a
single-covariate run (`"BAPS"`), and `+` joins several into one joint run
(`"RST+OspC+BAPS"`). This is what answers "does this hit survive once I specifically
remove the OspC effect" and, separately, "does it survive once I remove OspC *and* BAPS
together" — different, complementary questions, each an explicit, auditable pyseer
invocation you chose.

This is **not** an all-combinations powerset search, and that's intentional rather than a
missing feature. Compute isn't the constraint (each run here takes seconds); the
constraint is table 3 above — most combinations involving MLST are going to fail outright
or silently lose power, so enumerating all 15 non-empty subsets of four columns would
mostly generate combinations already known to be unusable, cluttering the output with
noise instead of the specific comparisons that answer a real question. Pick the groups
that matter (individual markers, plus perhaps the coarser ones jointly) rather than
generating everything.

Results land in `covariate_scan_results` (one TSV per entry, named `<entry>.tsv` — e.g.
`RST+OspC.tsv` — pyseer's native columns unchanged) and also in `covariate_scan_combined`
(all of them concatenated long-format with a `covariates` column, for a quick
side-by-side).

Why this matters for the dissemination/virulence question specifically: the kinship
matrix already corrects for genome-wide relatedness *continuously* — it doesn't tell you
which *named* lineage marker is responsible for an association. The combination scan
answers that sharper question directly. OspC major group in particular is a
well-established dissemination-associated marker in *B. burgdorferi* biology, so a module
that loses significance specifically in the OspC-adjusted run is more likely "this module
tags OspC type" than a genuinely new finding — worth checking first, not last.

### The annotation join

`annotation_table` needs an `id_column` (default `"Gene"`) whose values match pyseer's
`variant` column, plus optionally `name_column`/`annotation_column` (defaults match
Panaroo's own header names). The join is a **left** join — a `variant` with no match is
kept and labelled `"no match"`, never silently dropped, since the Rtab and the annotation
table could in principle come from different pangenome builds.

Two files come back, both **in addition to** pyseer's own output, which is never modified:

- `*_annotated.tsv` — pyseer's own columns, unchanged, with `gene_name`/`annotation` inserted
- `*_readable.tsv` — renamed headers (`gene_family_id`, `allele_frequency`,
  `prefilter_pvalue`, `lrt_pvalue`, `effect_size`, `effect_size_stderr`,
  `variance_explained`, `qc_flags`) and fixed-decimal formatting: `af`/`beta`/
  `beta-std-err`/`variant_h2` at 6 decimal places, and the two p-value columns at 6 decimal
  places above `0.0001`, scientific notation below it (so a p-value like `1e-15` doesn't
  round to `0.000000`)

Both are produced for `gene_results` and `gene_significant`. Leave `annotation_table` unset
and the workflow still runs — `*_annotated` is simply not produced, and `*_readable` still
comes back with renamed/reformatted headers but no gene name or annotation columns.

## Outputs

`gene_results` / `gene_significant` are the primary deliverables; `significance_threshold`
records the pattern-counted Bonferroni cutoff applied to `gene_significant`. Everything
else is conditional on the corresponding input being supplied.

## Running

```bash
miniwdl run workflows/pyseer_gwas/pyseer_gwas.wdl -i tests/pyseer_gwas_inputs.json
```

For a quick end-to-end check against a tiny synthetic 8-sample panel (no real data needed),
see `tests/pyseer_gwas_smoke_inputs.json` and `tests/run_unit_tests.sh`.
