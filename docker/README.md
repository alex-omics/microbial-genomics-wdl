# Docker images

Custom images for `microbial-genomics-wdl` workflows. All images are built
`--platform linux/amd64` for Terra/GCP compatibility and pushed to the
[`aarvani1`](https://hub.docker.com/u/aarvani1) Docker Hub namespace.

```bash
docker buildx build --platform linux/amd64 -t aarvani1/<image>:<tag> docker/<image>/
docker push aarvani1/<image>:<tag>
```

## multimodal_gwas image set

Each image corresponds to one `workflow/envs/*.yaml` conda environment in
[microGWAS](https://github.com/microbial-pangenomes-lab/microGWAS), the
Snakemake pipeline these WDL tasks are translated from. Where an upstream
environment pulls in more than one tool (e.g. `mash.yaml` needs both `mash`
*and* pyseer's `square_mash` helper), the image layers the extra pieces onto a
verified single-tool base rather than substituting a different version.

Every image that runs a `workflow/scripts/*.py` helper vendors the microGWAS
source at pinned commit
[`1307250`](https://github.com/microbial-pangenomes-lab/microGWAS/tree/130725016f023d78c13f08fdb961779afc2190ce)
(VERSION `0.9.1-dev`) into `/opt/microgwas`, with the scripts directory on
`PATH`. microGWAS is MIT-licensed, Copyright (c) 2022, Marco Galardini; the
license file travels with the vendored source at
`/opt/microgwas/LICENSE` in every image.

| Image | Upstream env (pinned version) | Base | Notes |
|---|---|---|---|
| `microgwas-base` | `environment.yml` | `mambaorg/micromamba` | Shared python/pandas/biopython/gffutils image for the ~24 rules with no dedicated env, plus `ncbi-genome-download` for reference bootstrapping |
| `microgwas-pyseer` | `pyseer.yaml` (`pyseer>=1.4.0`) | `mambaorg/micromamba` | pyseer 1.4.2 + bwa/bedtools for `map_back.py`. Supersedes `aarvani1/pyseer:1.4.0` |
| `microgwas-mash` | `mash.yaml` (`mash=2.1`) | `python:3.10-slim` + marbl mash 2.1 release binary | Adds pyseer for the `square_mash` helper used by the `distance` rule. Does *not* build on `staphb/mash:2.1`: that image is Ubuntu-16.04-era (Python 3.5, pip 8.1.1) and cannot install a current pandas, so the layering is inverted and mash comes from marbl's own release tarball at the pinned version |
| `microgwas-mlst` | `mlst.yaml` (`mlst=2.16`) | `staphb/mlst:2.16.2` | Adds GNU `parallel` and python for `sanitize_mlst.py` |
| `microgwas-snpsites` | `snp-sites.yaml` (`snp-sites=2.5.1`, `bcftools=1.13`) | `debian:bookworm-slim`, both tools from source | `aln2vcf` needs snp-sites *and* bcftools. Does not build on `staphb/snp-sites:2.5.1` — that image is EOL Ubuntu 21.04 whose apt repos 404 |
| `microgwas-bcftools` | `bcftools.yaml` (`bcftools=1.13`) | `python:3.10-slim` + bcftools 1.13 from source | Adds biopython/pandas/pysam for `vcf2deleterious.py`. Does not build on `staphb/bcftools:1.13` — Ubuntu 16.04, pip 8.1.1, no setuptools |
| `microgwas-snippy` | `nucmer.yaml` (`snippy=4.6.0`) | `staphb/snippy:4.6.0` | Composite upstream env; snippy vendors its own mummer/snpeff/bcftools |
| `microgwas-panfeed` | `panfeed.yaml` (`panfeed>=1.6.1`) | `quay.io/biocontainers/panfeed` | Helper scripts arrive via multi-stage `COPY --from`; the base has no package manager |
| `microgwas-enrich` | `enrich.yaml` | `mambaorg/micromamba` | goatools + scipy/statsmodels |
| `microgwas-limix` | `limix.yaml` (`limix=3.0.4`) | `mambaorg/micromamba` | Also vendors the `albi` submodule needed by `run_heritability` |
| `microgwas-eggnog` | `eggnog-mapper.yaml` | `quay.io/biocontainers/eggnog-mapper` | Adds the `eggnog-mapper-fixurl` pip patch |

Images pulled and used **unmodified** (no Dockerfile here):

| Image | Upstream env | Used by |
|---|---|---|
| `quay.io/biocontainers/panaroo:1.5.0--pyhdfd78af_0` | `panaroo.yaml` | `pangenome` |
| `quay.io/biocontainers/unitig-counter:1.1.0--h56fc30b_0` | `unitig-counter.yaml` | `unitigs` |
| `staphb/fasttree:2.1.11` | `fasttree.yaml` | `tree` |

### Version pins that could not be matched exactly

Per the milestone brief's version-tag fallback rule, these are flagged rather
than silently substituted — see `docs/workflows/multimodal_gwas.md` for detail.

- **`abritamr=1.0.14`** — no such version is published on bioconda
  (biocontainers jumps 1.0.9 → 1.1.0). The `find_amr_vag` rule it backs is a
  side-branch that produces no input to any association task, so no image is
  built here and the task defaults to a documented near-version.

### Build status

Verified on an amd64 build from a network that can reach quay.io, Docker Hub
and PyPI but **not** `conda.anaconda.org` (blocked on the MGH network).

| Image | Builds | Pushed | Notes |
|---|---|---|---|
| `microgwas-mash:2.1` | ✅ | ✅ | Also functionally tested: `mash sketch` → `mash dist` → `square_mash` yields a square matrix labelled by bare sample ID |
| `microgwas-mlst:2.16.2` | ✅ | ✅ | |
| `microgwas-snippy:4.6.0` | ✅ | ✅ | |
| `microgwas-eggnog:2.1.13` | ✅ | ✅ | |
| `microgwas-bcftools:1.13` | ✅ | — | |
| `microgwas-snpsites:2.5.1` | ✅ | — | |
| `microgwas-panfeed:1.6.1` | ✅ | — | |
| `microgwas-base:0.9.1` | ⛔ conda | — | needs anaconda.org |
| `microgwas-pyseer:1.4.2` | ⛔ conda | — | needs anaconda.org |
| `microgwas-enrich:0.9.1` | ⛔ conda | — | needs anaconda.org |
| `microgwas-limix:3.0.4` | ⛔ conda | — | needs anaconda.org |

### Base images that could not be used as-is

Three StaPH-B / biocontainers images turned out to be unusable as bases, each
for a different reason. All three are documented in the relevant Dockerfile
header; recorded here so the reasoning is not lost:

- **`staphb/bcftools:1.13`** — Ubuntu 16.04 (Xenial), ships pip 8.1.1 with no
  setuptools. Cannot install biopython or pandas at any version.
- **`staphb/snp-sites:2.5.1`** — Ubuntu 21.04 ("hirsute"), which is end of
  life. `archive.ubuntu.com` 404s for it, so `apt-get update` fails outright.
- **`staphb/mash:2.1`** — same Xenial-era pip problem as bcftools.

In each case the layering was inverted: a current base, with the pinned tool
built from its own upstream release (marbl for mash, samtools for bcftools,
the `v2.5.1` git tag for snp-sites). This keeps the version pins exact rather
than drifting to whatever a distro or conda channel resolves.

**`quay.io/biocontainers/panfeed:1.6.1`** has no package manager at all — no
apt, no conda, no micromamba. The microGWAS helper scripts are fetched in a
throwaway builder stage and `COPY --from`'d in, which places no requirements
on the runtime base.

## busco-prokaryota

See [`busco-prokaryota/`](busco-prokaryota/).
