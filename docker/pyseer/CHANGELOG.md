# Changelog

## 1.4.2 — 5 August 2026
- Initial build of the standalone pyseer image.
- pyseer 1.4.2. Dependency pins: python 3.9, `numpy<1.23.0`, `scipy==1.7.0` (required for glmnet
  C-ABI compatibility; scipy 1.7.0 has no cp310 wheel, hence python 3.9).
- Also included: `matplotlib-base`, `pybedtools` and `bedtools` (pyseer's own `install_requires`),
  and `pandas` for the annotation join in `tasks/pyseer.wdl`.
- Base: `mambaorg/micromamba:1.5.10-bookworm-slim`
- Digest: sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c
  - Pin with `aarvani1/pyseer@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c`
