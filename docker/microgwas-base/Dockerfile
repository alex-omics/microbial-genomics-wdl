# microgwas-base
#
# Shared image for the ~24 microGWAS rules that have no dedicated conda
# environment and inherit the pipeline's base `environment.yml`:
# lightweight python + pandas + biopython work, plus the reference-genome
# bootstrapping that upstream does in `bootstrap.sh`.
#
# Backs these WDL tasks: bootstrap_references, prepare_regions,
# sample_whole_pangenome, lineages2covariance, combine_heritability,
# gpa_summary, rare_summary, map_summary (all variants), annotate_summary
# (all variants), annotate_reference, download_obo, combine_panfeed.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="Shared base image for microGWAS-derived WDL tasks (pandas/biopython/gffutils + microGWAS helper scripts)"
LABEL version="0.9.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

# Terra runs containers as root; micromamba's default USER is not.
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    gzip \
    wget \
    curl \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# environment.yml: python>=3.6, ncbi-genome-download, biopython, gffutils.
# pandas is pulled in by nearly every workflow/scripts/*.py helper.
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.10 \
    pandas \
    numpy \
    scipy \
    biopython \
    gffutils \
    ncbi-genome-download \
    pyyaml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# Vendor the microGWAS source at a pinned commit so script behaviour is
# reproducible across image rebuilds. MIT-licensed, Copyright (c) 2022,
# Marco Galardini — LICENSE travels with the source at /opt/microgwas/LICENSE.
ARG MICROGWAS_COMMIT=130725016f023d78c13f08fdb961779afc2190ce
RUN git clone https://github.com/microbial-pangenomes-lab/microGWAS.git /opt/microgwas \
    && cd /opt/microgwas \
    && git checkout ${MICROGWAS_COMMIT} \
    && rm -rf .git \
    && chmod +x workflow/scripts/*.py

# Put the helper scripts on PATH so WDL command blocks can call them by name
# rather than by the `python workflow/scripts/X.py` path the Snakefile uses.
ENV PATH="/opt/microgwas/workflow/scripts:${PATH}"
ENV MICROGWAS_SCRIPTS="/opt/microgwas/workflow/scripts"
# map_back.py does `from bwa import bwa_index`, and several other helpers
# import their siblings, so the scripts dir must be importable, not just on PATH
ENV PYTHONPATH="/opt/microgwas/workflow/scripts:${PYTHONPATH:-}"

RUN python -c "import pandas, Bio, gffutils; print('base deps OK')" \
    && ncbi-genome-download --version \
    && test -f /opt/microgwas/workflow/scripts/prepare_pyseer.py \
    && test -f /opt/microgwas/LICENSE \
    && echo "=== All checks passed ==="

WORKDIR /data
