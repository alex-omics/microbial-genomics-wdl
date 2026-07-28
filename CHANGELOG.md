# Changelog

All notable changes to this repository are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`multimodal_gwas` workflow (Milestone 1: Snakemake → WDL translation).**
  Bacterial GWAS across unitigs, gene presence/absence, structural variants,
  gene-cluster-specific k-mers, common SNPs and a whole-genome elastic net,
  with population-structure correction via pyseer's LMM. Translated from
  [microGWAS](https://github.com/microbial-pangenomes-lab/microGWAS) at pinned
  commit `1307250`.
  - `workflows/multimodal_gwas/wf_multimodal_gwas.wdl` — orchestration, with
    the Snakemake DAG reconstructed as an explicit WDL call graph.
  - `tasks/gwas/` — 11 task files covering 40+ upstream rules, grouped by
    upstream conda environment. Shell commands and script invocations preserved
    verbatim except where documented in the task files.
  - `docker/microgwas-*` — 11 Dockerfiles, one per upstream environment.
  - `docs/workflows/multimodal_gwas.md` — usage, translation notes, validation
    plan, and carried-forward upstream limitations.
  - `tests/multimodal_gwas/` — Tier 1 (smoke) and Tier 2 (E. coli statistical
    validation) input templates.

### Changed

- **`aarvani1/pyseer:1.4.0` superseded by `aarvani1/microgwas-pyseer:1.4.2`.**
  Forward bump within upstream's `pyseer>=1.4.0` pin. The new image also
  vendors the microGWAS helper scripts and adds bwa/bedtools, which
  `map_back.py` shells out to.

### Notes

- Three upstream version pins could not be matched exactly and are flagged
  rather than silently substituted: `abritamr=1.0.14` (does not exist on
  bioconda; no Tier 2 impact — the rule feeds no association task),
  `mlst=2.16` (using StaPH-B `2.16.2`), and `limix=3.0.4` (unpublished on
  biocontainers; built from the conda recipe). See
  `docs/workflows/multimodal_gwas.md`.
- Sequence UNET, and therefore the rare-variant gene burden test, is
  deliberately out of scope for Milestone 1 — the microGWAS paper reports the
  burden test as the one modality that did not reproduce the earlier
  validation results.
- Neither validation tier has been executed yet. The WDL type-checks under
  `miniwdl check`; that is not the same as running.
