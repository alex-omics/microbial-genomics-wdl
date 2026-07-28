# microgwas-limix
#
# Backs the microGWAS `run_heritability` rule (workflow/envs/limix.yaml:
# limix=3.0.4, pandas<1.1.0).
#
# Two things make this image non-trivial:
#
#  1. limix is not published on quay.io/biocontainers at all, so there is no
#     image to reuse - it has to be built from the conda recipe.
#  2. The rule also runs `python albi/albi.py`, which lives in the `albi`
#     *git submodule* of microGWAS (github.com/cozygene/albi), not in
#     workflow/scripts. Cloning microGWAS non-recursively would silently
#     produce an empty albi/ and the heritability confidence intervals would
#     come back blank (the rule swallows this with `|| true`). So this image
#     clones with --recursive and verifies albi.py is really there.
#
# limix 3.0.4 requires an older scipy/pandas stack; the pins below are what
# actually resolve rather than the newest available.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="limix 3.0.4 + albi, for microGWAS narrow-sense heritability estimation"
LABEL version="3.0.4"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    wget \
    git \
    build-essential \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# limix 3.0.4 predates python 3.9 wheels for several of its deps
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.8 \
    "pandas<1.1.0" \
    "numpy<1.21" \
    scipy \
    limix=3.0.4 \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# --recursive is load-bearing here: albi is a submodule, see header comment
ARG MICROGWAS_COMMIT=130725016f023d78c13f08fdb961779afc2190ce
RUN git clone --recursive https://github.com/microbial-pangenomes-lab/microGWAS.git /opt/microgwas \
    && cd /opt/microgwas \
    && git checkout ${MICROGWAS_COMMIT} \
    && git submodule update --init --recursive \
    && rm -rf .git \
    && chmod +x workflow/scripts/*.py

ENV PATH="/opt/microgwas/workflow/scripts:${PATH}"
ENV MICROGWAS_SCRIPTS="/opt/microgwas/workflow/scripts"
# map_back.py does `from bwa import bwa_index`, and several other helpers
# import their siblings, so the scripts dir must be importable, not just on PATH
ENV PYTHONPATH="/opt/microgwas/workflow/scripts:${PYTHONPATH}"
ENV MICROGWAS_ALBI="/opt/microgwas/albi"

RUN python -c "import limix; print('limix', limix.__version__)" \
    && test -f /opt/microgwas/albi/albi.py \
    && test -f /opt/microgwas/workflow/scripts/estimate_heritability.py \
    && test -f /opt/microgwas/workflow/scripts/prepare_fiesta.py \
    && echo "=== All checks passed ==="

WORKDIR /data
