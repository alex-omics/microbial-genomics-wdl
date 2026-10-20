# Changelog

## 1.2.0 — 20 October 2026
- Initial build of the LISTT image, tracking LISTT release 1.2.0 (commit `b9c6569`, scheme-update-1.2).
- LISTT is cloned at that commit to `/opt/LISTT` (its `.git` removed, `COMMIT` records the hash); the build fails if
  the checkout differs. `LISTT_VERSION=1.2.0` is set in the environment.
- Tool versions follow LISTT's `environment.yml`: bwa 0.7.17, samtools/bcftools/htslib 1.18, blast 2.15.0,
  mafft 7.520, python 3.11, biopython 1.81, numpy 1.26, pandas 2.1, screed 1.1.2, pyvcf 0.6.8. The paper used
  bcftools 1.17 and mafft 7.480.
- PyVCF 0.6.8 only builds with `setuptools<58` and `--no-build-isolation` (its `setup.py` uses `use_2to3`).
- `audit_listt_tables.py` checks LISTT's reference tables and re-derives the per-IST thresholds.
- Base: `mambaorg/micromamba:1.5.10-bookworm-slim`, built for `linux/amd64`.
- Digest: sha256:c016417e25c3e4b96d4fc98bedbf4d77e15608674e5e5320cdaff1df6a1a38b2
  - Pin with `aarvani1/listt@sha256:c016417e25c3e4b96d4fc98bedbf4d77e15608674e5e5320cdaff1df6a1a38b2`
