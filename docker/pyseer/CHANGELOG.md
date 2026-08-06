# Changelog

## 1.4.2 — 5 August 2026
- Initial build of the standalone image, split off from the earlier
  `aarvani1/pyseer:1.4.0` build used by the original `task_pyseer.wdl` /
  `wf_pyseer_gwas.wdl` sketch.
- Bumped pyseer 1.4.0 → 1.4.2 (released upstream since the previous build).
- Reuses the dependency pins already validated for `microgwas-pyseer:1.4.2`
  (python 3.9, `numpy<1.23.0`, `scipy==1.7.0` — required for the glmnet
  C-ABI compatibility; scipy 1.7.0 has no cp310 wheel, hence python 3.9 and
  not 3.10) rather than re-deriving them, but trims what's specific to the
  microGWAS DAG only: no eggNOG-mapper deps, no vendored microGWAS scripts,
  no seaborn/mash/bwa/bedops. `matplotlib-base`, `pybedtools` and `bedtools`
  are kept — the first build attempt dropped them as apparently
  microGWAS-only and failed, since they're actually pyseer's own
  `install_requires`, not extras. Also kept: `pandas`, for the
  annotation-join step in `tasks/pyseer.wdl`'s `pyseer_annotate_results`.
- Base: `mambaorg/micromamba:1.5.10-bookworm-slim`
- Built and pushed to Docker Hub on 6 August 2026, after the previous attempt
  was blocked by `conda.anaconda.org` being unreachable on the MGH network
  (see `docker/build_and_push.sh`'s `NEEDS_CONDA` note) — built from a network
  that could reach it instead.
- Digest: sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c
  - Verified against the registry, not just locally. Pin with:
    `aarvani1/pyseer@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c`
