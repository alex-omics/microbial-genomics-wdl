# pyseer_gwas

Bacterial GWAS with [pyseer](https://github.com/mgalardini/pyseer), independently callable
in Terra. Runs pyseer's LMM (likelihood-ratio test) over a presence/absence Rtab, correcting
for population structure with a phylogeny-derived kinship matrix.

`presence_absence_rtab` is a generic `block_id x sample` 0/1 matrix, so a Panaroo
`gene_presence_absence.Rtab` and a pangenome-network module Rtab both work unchanged.

## What it does

| Step | Tool | Output |
| ---- | ---- | ------ |
| Kinship + distance matrix | `phylogeny_distance.py`, with and without `--lmm` | `kinship_matrix`, `distance_matrix` |
| Association | `pyseer --lmm` | `gene_results`, `gene_significant` (pattern-counted Bonferroni via `count_patterns.py`) |
| Lineage effects *(optional)* | `pyseer --lineage` | `lineage_effects` |
| Covariate combination scan *(optional)* | `pyseer --use-covariates`, once per combination | `covariate_scan_results`, `covariate_scan_combined` |
| SNP association *(optional)* | `pyseer --vcf` | `snp_results` |
| Name/annotation join *(optional)* | pandas left join; pyseer's own files are never modified | `gene_results_annotated`, `gene_results_readable`, and the same pair for `covariate_scan_combined` |

`gene_significant` is a subset of `gene_results`, so it is not annotated separately; filter
`gene_results_annotated` instead.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `phylogeny_newick` | `File` | **Required.** Midpoint-rooted core-genome phylogeny |
| `presence_absence_rtab` | `File` | **Required.** Gene or module presence/absence matrix |
| `phenotype_tsv` | `File` | **Required.** `sample_id\tphenotype_value` |
| `variant_vcf` | `File?` | Adds a SNP-based association |
| `covariates_file` | `File?` | `sample_id` plus one named column per covariate (e.g. MLST, BAPS) |
| `use_covariates` | `String?` | pyseer's raw `--use-covariates` value, applied jointly to the main association |
| `covariate_combinations` | `Array[String]` | Groups of `covariates_file` column names to test together, one pyseer run per entry (see below) |
| `run_lineage_effects` | `Boolean` | Reports per-lineage effects in a separate task that runs `--lineage` on only a small slice of `presence_absence_rtab`, enough for pyseer to initialise. Adding `--lineage` to the full `--lmm` call can exhaust memory at scale. Default `false`. The required distance matrix is derived from `phylogeny_newick` automatically |
| `lineage_clusters` | `File?` | `sample_id\tcluster_id` (e.g. a BAPS export), for `--lineage-clusters` |
| `annotation_table` | `File?` | Panaroo/Roary `gene_presence_absence.csv`, or an equivalent module-annotation table |

## Covariates

Three points that fail silently rather than with an error:

1. **Supplying `covariates_file` does nothing by itself.** pyseer loads covariates but does
   not use them unless told to. Set `use_covariates` or `covariate_combinations`.
2. **Lineage and sequence-type columns (MLST, BAPS, ...) are categorical**, even when stored
   as integers: they are group labels, so never add pyseer's `q` (quantitative) suffix.
   `covariate_combinations` always treats columns as categorical.
3. **Joint categorical covariates use up degrees of freedom.** A k-level column costs about
   k-1 fitted parameters, drawn from your sample count. With few samples or many levels, a
   joint model can leave no variant able to fit, which shows up as an empty results table.
   The finest-grained column (most levels, typically MLST) is the first to cause this.

### `covariate_combinations`

Each entry runs one pyseer call using exactly that group of columns as the covariates: a bare
name is a single-covariate run (`"BAPS"`), and `+` joins several into one joint run
(`"MLST+BAPS"`). This shows whether an association survives adjustment for a specific lineage
marker, or for a specific group of them.

It is deliberately not an all-combinations search. Many combinations would exhaust degrees of
freedom (point 3) and only add noise; choose the groups that answer a real question.

Results land in `covariate_scan_results` (one TSV per entry, named `<entry>.tsv`, pyseer's
columns unchanged) and `covariate_scan_combined` (all of them in long format with a
`covariates` column).

### The annotation join

`annotation_table` needs an `id_column` (default `"Gene"`) whose values match pyseer's
`variant` column, and optionally `name_column` / `annotation_column` (defaults match
Panaroo's headers). It is a **left** join: a variant with no match is kept and labelled
`"no match"`, never dropped.

Two files come back in addition to pyseer's own output, which is never modified:

- `*_annotated.tsv` — pyseer's columns unchanged, with `gene_name` / `annotation` inserted
- `*_readable.tsv` — renamed headers (`gene_family_id`, `allele_frequency`,
  `prefilter_pvalue`, `lrt_pvalue`, `effect_size`, `effect_size_stderr`, `variance_explained`,
  `qc_flags`) and fixed-decimal formatting. p-values use scientific notation below `0.0001`,
  so a value like `1e-15` does not round to `0.000000`

Without `annotation_table`, `*_annotated` is not produced and `*_readable` has no gene name
or annotation columns.

## Outputs

`gene_results` and `gene_significant` are the primary outputs; `significance_threshold` is
the pattern-counted Bonferroni cutoff applied to `gene_significant`. Everything else depends
on the corresponding input being supplied.

## Running

```bash
miniwdl run workflows/pyseer_gwas/pyseer_gwas.wdl -i tests/pyseer_gwas_inputs.json
```

For a quick end-to-end check on a small synthetic panel, see
`tests/pyseer_gwas_smoke_inputs.json` and `tests/run_unit_tests.sh`.
