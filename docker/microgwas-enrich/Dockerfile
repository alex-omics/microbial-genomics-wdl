# microgwas-enrich
#
# Backs the microGWAS `run_enrich`, `run_enrich_alt` and `run_enrich_wg` rules
# (workflow/envs/enrich.yaml: scipy, pandas, goatools, statsmodels>=0.10).
#
# Not a published tool - the environment exists to run
# workflow/scripts/functional_enrichment.py, which does Fisher's exact tests
# over COG categories, GO terms and KEGG pathways with FDR correction
# (statsmodels) and GO DAG traversal (goatools).
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="goatools + scipy/statsmodels for microGWAS COG/GO/KEGG functional enrichment"
LABEL version="0.9.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    wget \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.10 \
    scipy \
    pandas \
    numpy \
    goatools \
    "statsmodels>=0.10" \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

ARG MICROGWAS_COMMIT=130725016f023d78c13f08fdb961779afc2190ce
RUN git clone https://github.com/microbial-pangenomes-lab/microGWAS.git /opt/microgwas \
    && cd /opt/microgwas \
    && git checkout ${MICROGWAS_COMMIT} \
    && rm -rf .git \
    && chmod +x workflow/scripts/*.py

ENV PATH="/opt/microgwas/workflow/scripts:${PATH}"
ENV MICROGWAS_SCRIPTS="/opt/microgwas/workflow/scripts"
# map_back.py does `from bwa import bwa_index`, and several other helpers
# import their siblings, so the scripts dir must be importable, not just on PATH
ENV PYTHONPATH="/opt/microgwas/workflow/scripts:${PYTHONPATH}"

RUN python -c "import goatools, scipy, pandas, statsmodels; print('enrich deps OK')" \
    && test -f /opt/microgwas/workflow/scripts/functional_enrichment.py \
    && echo "=== All checks passed ==="

WORKDIR /data
