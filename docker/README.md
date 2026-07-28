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
| `microgwas-snpsites` | `snp-sites.yaml` (`snp-sites=2.5.1`, `bcftools=1.13`) | `staphb/snp-sites:2.5.1` | Adds bcftools 1.13 for the `aln2vcf` rule |
| `microgwas-bcftools` | `bcftools.yaml` (`bcftools=1.13`) | `staphb/bcftools:1.13` | Adds biopython/pandas/pysam for `vcf2deleterious.py` |
| `microgwas-snippy` | `nucmer.yaml` (`snippy=4.6.0`) | `staphb/snippy:4.6.0` | Composite upstream env; snippy vendors its own mummer/snpeff/bcftools |
| `microgwas-panfeed` | `panfeed.yaml` (`panfeed>=1.6.1`) | `quay.io/biocontainers/panfeed` | |
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

| Image | Built | Functionally tested |
|---|---|---|
| `microgwas-mash:2.1` | ✅ | ✅ `mash sketch` → `mash dist` → `square_mash` produces a square matrix labelled by bare sample ID |
| `microgwas-mlst:2.16.2` | ✅ | Not yet |
| all others | Not yet | Not yet |

The conda-based images (`base`, `pyseer`, `enrich`, `limix`) could not be built
in the session that authored them, because `conda.anaconda.org` was unreachable
from that network. Build them before the first Tier 1 run.

## busco-prokaryota

See [`busco-prokaryota/`](busco-prokaryota/).
