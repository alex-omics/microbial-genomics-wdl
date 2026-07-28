# multimodal_gwas

**Status:** 🚧 In development (Milestone 1 — translation complete, validation pending)

Bacterial genome-wide association study across multiple classes of genetic
variation, run natively on Terra.

This workflow is a WDL translation of
[microGWAS](https://github.com/microbial-pangenomes-lab/microGWAS), the
Snakemake pipeline described in Burgaya et al. 2025, *microGWAS: a
computational pipeline to perform large-scale bacterial genome-wide association
studies*, Microbial Genomics 11:001349. It is a re-implementation of the
orchestration, not of the science: every underlying tool (pyseer, Panaroo,
panfeed, unitig-counter, snippy, mash, eggNOG-mapper) is called with upstream's
own arguments.

## What it tests

| Variant class | Tool | Association |
|---|---|---|
| Unitigs | unitig-counter | pyseer LMM |
| Gene presence/absence | Panaroo | pyseer LMM |
| Structural variants | Panaroo | pyseer LMM |
| Gene-cluster-specific k-mers | panfeed | pyseer LMM |
| Common SNPs vs. reference | snippy + bcftools | pyseer LMM |
| Rare variants (gene burden) | snippy + bcftools + Sequence UNET | pyseer LMM — *off by default, see below* |
| Whole genome | — | pyseer elastic net (ridge α=0.01, lasso α=1) |

Population structure is corrected via a core-genome-SNP kinship matrix
(pyseer's LMM), with mash distances and MLST lineages as additional
covariates. The significance threshold is set by counting unique
presence/absence patterns rather than by naive Bonferroni.

Downstream: bwa/bedtools mapping-back to genes, eggNOG-mapper annotation,
COG/GO/KEGG enrichment (Fisher's exact + FDR), Manhattan and QQ plots,
per-gene hit summaries, and narrow-sense heritability.

## Inputs

Minimum viable run:

| Input | Description |
|---|---|
| `samples_tsv` | Strain ID in column 1, plus one column per phenotype in `targets` |
| `sample_fastas` | Assemblies named `SAMPLE.fasta`, matching the strain IDs |
| `targets` | Phenotype column names to run associations for |
| `genus`, `species`, `reference_strain` | For reference genome bootstrapping |
| `assembly_accessions` | NCBI RefSeq accessions for the reference genomes |
| `summary_references`, `annotation_references`, `enrichment_reference` | Reference strain flag strings, as in upstream's config |
| `mlst_scheme` | e.g. `ecoli` |

Naming matters: `SAMPLE.fasta` / `SAMPLE.gff` basenames are how sample IDs are
recovered throughout, because that is how upstream identifies strains. A
mismatch between a FASTA basename and its row in `samples_tsv` will silently
drop that sample at the `prepare_pyseer` intersection step (the task fails
loudly if the intersection empties entirely).

Worth setting:

- **`eggnog_db_tarball`** — run `download_eggnog` once, park the tarball in a
  GCS bucket, and pass it here. Leaving it unset re-downloads tens of GB every
  run.
- **`sample_gffs`** — supply pre-computed annotations to skip ggCaller.
- **`lineages_file`** — bring your own lineage assignments instead of MLST.
- **`length` / `min_hits` / `max_genes`** — spurious-hit filters. Upstream's
  config defaults (10 / 1 / 25) are deliberately loose; the paper's own
  *E. coli* analysis used 30 / 9 / 10.

## Repository layout

```
workflows/multimodal_gwas/wf_multimodal_gwas.wdl   # orchestration + call graph
tasks/gwas/task_setup.wdl                          # reference bootstrap, regions, GO DAG
tasks/gwas/task_variants.wdl                       # unitigs, Panaroo, panfeed, ggCaller
tasks/gwas/task_population_structure.wdl           # mash, MLST, FastTree, snp-sites, kinship
tasks/gwas/task_rare_variants.wdl                  # snippy, bcftools rare/common split
tasks/gwas/task_pyseer.wdl                         # the 15 pyseer-environment rules
tasks/gwas/task_mapback.wdl                        # map-back, Manhattan plot
tasks/gwas/task_heritability.wdl                   # limix + ALBI
tasks/gwas/task_annotation.wdl                     # eggNOG-mapper
tasks/gwas/task_summaries.wdl                      # per-gene summaries, annotation joins
tasks/gwas/task_enrichment.wdl                     # COG/GO/KEGG enrichment + plots
tasks/gwas/task_amr.wdl                            # optional abritamr side branch
docker/microgwas-*/Dockerfile                      # one image per upstream conda env
tests/multimodal_gwas/*.inputs.json                # Tier 1 and Tier 2 input templates
```

## Translation notes

Things a reader comparing this against the Snakefile should know.

**Manifests are rebuilt per task.** Upstream writes path manifests once
(`out/mash_input.txt`, `out/unitigs_input.tsv`, …) and every rule reads them,
which works because all rules share a working directory. Under Cromwell each
task gets its own container and its own localized paths, so a manifest written
by one task is meaningless to the next. Each task therefore reconstructs what
it needs from its `Array[File]` inputs. This is safe because every tool
involved derives sample IDs from file basenames — verified for `square_mash`
(`os.path.split(...)[-1].split('.')[0]`), `sanitize_mlst.py`, and
`map_back.py`.

**The 12 collector rules are gone, by design.** `manhattan_plots`, `wg`,
`qq_plots`, `map_back`, `map_summary`, `annotate_summary`, `enrichment`,
`enrichment_plots`, `pyseer`, `pyseer_vcf`, `heritability` and `panfeed` have
no `shell:` body; they exist only to force Snakemake to materialise scattered
outputs via `expand()`. A WDL task consuming an `Array[File]` from a scatter
already means "wait for all of these".

**Some upstream rule inputs were dropped.** `run_panfeed` declares `distances`
and `lineages`, `run_pyseer_vcf` declares `lineages`, and `run_wg` declares
`similarity` — none of which their command bodies read. In Snakemake these
exist to order the DAG. Cromwell orders on the explicit call graph, so keeping
them would mean localizing large matrices into tasks that ignore them.

**Five `run_annotate_*` rules became one task.** Their shell bodies are
identical apart from paths, so they collapse to a single parameterised
`annotate_summary`. Same for the three `run_enrich*` and three
`run_enrichment_plots*` rules. The `map_summary` variants genuinely differ in
their flags and stay separate.

**One upstream quirk preserved deliberately.** `run_pyseer` filters the gene
presence/absence and structural variant results against the *unitig* pattern
threshold, not their own — even though it writes `gpa_patterns.txt` and
`struct_patterns.txt`. This looks like a bug, but changing it would shift which
hits pass and break comparability with the published results, so it is
reproduced verbatim and flagged in the task file.

**Shell strictness.** Upstream's `shell:` blocks run without `pipefail`, and a
couple of them depend on that: `grep -v 'cv.glmnet'` legitimately filters every
line when glmnet is quiet, and `head -1` on a large file SIGPIPEs its producer.
Those pipelines are decomposed into separate steps so `set -euxo pipefail` can
stay on without changing behaviour.

## Docker images

One image per upstream `workflow/envs/*.yaml`. See [`docker/README.md`](../../docker/README.md)
for the full table and build commands. Two things worth calling out:

- **`aarvani1/pyseer:1.4.0` is superseded** by `aarvani1/microgwas-pyseer:1.4.2`.
  Upstream pins `pyseer>=1.4.0`, so this is a forward bump within the pinned
  range. The new image also vendors the microGWAS helper scripts and adds
  bwa/bedtools, which `map_back.py` shells out to.
- **Several upstream "single-tool" environments are not single-tool.**
  `mash.yaml` also requires pyseer (for the `square_mash` helper the `distance`
  rule pipes into); `snp-sites.yaml` also requires bcftools; `mlst.yaml` also
  requires GNU `parallel`. Pulling the stock StaPH-B image alone would leave a
  pipeline stage missing, so those images layer the extra pieces on.

### Version pins that could not be matched exactly

Flagged rather than silently substituted, per the milestone brief:

| Pin | Situation | Impact |
|---|---|---|
| `abritamr=1.0.14` | No such version exists on bioconda — published builds go 1.0.9 → 1.1.0. Defaulted to 1.0.9. | **None on Tier 2.** `find_amr_vag` is a side branch whose output feeds no association task, and it is off by default. |
| `mlst=2.16` | No bioconda build (2.9 → 2.33.1). Using StaPH-B `2.16.2`, the patch release of the pinned line. | Low — same minor version and scheme database generation. |
| `limix=3.0.4` | Not published on biocontainers at all; built from the conda recipe. | Affects heritability only, not association results. |

All other pins matched exactly: `fasttree=2.1.11`, `bcftools=1.13`,
`snp-sites=2.5.1`, `snippy=4.6.0`, `mash=2.1`, `panaroo>=1.5.0`,
`unitig-counter=1.1.0`, `panfeed>=1.6.1`, `eggnog-mapper>=2.1.6`.

## Validation

Two tiers, per the Milestone 1 definition of done.

**Tier 1 — mechanical correctness.** Run against microGWAS's own bundled small
test dataset (`test/small_fastas.tgz`, `test/stripped_small_gffs.tgz`,
`test/test_data.tsv`). Confirms the DAG executes end-to-end and scatter logic
resolves. Template: [`tests/multimodal_gwas/tier1_smoke_test.inputs.json`](../../tests/multimodal_gwas/tier1_smoke_test.inputs.json).

Note that upstream's own CI only ever ran a Snakemake **dry run** — the
pipeline has never been verified end-to-end by its own test suite. Tier 1 here
needs to be a genuine execution.

**Tier 2 — statistical fidelity.** Reproduce the paper's *E. coli* virulence
validation on the 370-strain public dataset. Template:
[`tests/multimodal_gwas/tier2_ecoli_validation.inputs.json`](../../tests/multimodal_gwas/tier2_ecoli_validation.inputs.json).

Pass criteria: recover the three iron-uptake-system associations (HPI,
aerobactin, *sitABCD* operon), the core-gene hits (*zinT*, *mtfA*, *shiA*), and
a lasso whole-genome R² in the neighbourhood of 0.48. Exact numbers will drift
— tool versions have moved since publication. **Getting the same loci is the
bar, and that judgment requires a human reading the Manhattan plots and hit
tables.** It is not something this workflow can self-certify.

Read: `annotated_summary_unitigs`, `annotated_summary_gpa`, `manhattan_png`,
and `wg_metrics`.

## Known limitations carried forward from upstream

- **Heritability assumes normally distributed errors.** The microGWAS
  maintainers flag this as potentially inappropriate for binary phenotypes —
  which is what most clinical-manifestation phenotypes are. Treat these
  estimates as indicative.
- **The Manhattan plot script does not correctly handle multi-chromosome or
  plasmid-bearing references.** This is an open upstream limitation, and a
  directly relevant one for *Borrelia*.
- **`square_mash` truncates strain names at the first dot.** A reference named
  `K-12_substr._MG1655.fasta` becomes `K-12_substr` in the distance matrix.
  Harmless for references (they carry no phenotype and drop out at the
  intersection step), but avoid dots in study sample names.
- **Covariate handling looks broken upstream and is untested here.**
  `config["covariates"]` documents 1-based column indices into the phenotype
  file (e.g. `--use-covariates 6q 7`), but `prepare_pyseer.py` reduces that
  file to just the index plus the single target column before pyseer sees it —
  so those columns no longer exist by the time the flag is applied. Upstream's
  CI never caught this because it only ever ran a Snakemake dry run. The
  `covariates` input is wired through faithfully, but treat it as unverified.
  This does not affect Tier 2: the *E. coli* validation dataset has no
  covariate columns.

## Out of scope for Milestone 1

Deferred deliberately, not overlooked: Sequence UNET deleteriousness scoring
and therefore the rare-variant gene burden test; `bac-nortia` as an alternative
whole-genome model; PopPUNK lineages; hogwash/treeWAS phylogenetic-convergence
cross-checks; the phenotype manifest/long-format system; new variant classes
(copy-number, mobile elements, plasmid architecture).

## Attribution

This workflow is a translation of
[microGWAS](https://github.com/microbial-pangenomes-lab/microGWAS), MIT-licensed,
Copyright (c) 2022, Marco Galardini. Helper scripts from
`workflow/scripts/` are vendored unmodified into the Docker images at pinned
commit [`1307250`](https://github.com/microbial-pangenomes-lab/microGWAS/tree/130725016f023d78c13f08fdb961779afc2190ce),
with the upstream LICENSE alongside them at `/opt/microgwas/LICENSE`.

If you use this workflow, cite the original:

> Burgaya J, Damaris BF, Fiebig J, Galardini M. microGWAS: a computational
> pipeline to perform large-scale bacterial genome-wide association studies.
> *Microbial Genomics*. 2025;11(2):001349.

The *E. coli* validation phenotypes come from:

> Galardini M, Clermont O, Baron A, et al. Major role of iron uptake systems in
> the intrinsic extra-intestinal virulence of the genus *Escherichia* revealed
> by genome-wide association studies. *PLoS Genetics*. 2020;16(10):e1009065.
