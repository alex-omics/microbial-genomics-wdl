# microgwas-panfeed
#
# Backs the microGWAS `panfeed_kmers` and `panfeed_downstream` rules
# (workflow/envs/panfeed.yaml, which pins `panfeed>=1.6.1`).
#
# panfeed_downstream needs the `panfeed-get-clusters`, `panfeed-get-kmers` and
# `panfeed-plot` console scripts in addition to `panfeed` itself; all ship in
# the same distribution.
FROM quay.io/biocontainers/panfeed:1.6.1--pyhdfd78af_0

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="panfeed 1.6.1 + microGWAS helper scripts, for gene-cluster-specific k-mer generation"
LABEL version="1.6.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

# biocontainers images are minimal; git is needed to vendor the helper scripts
RUN if command -v apt-get > /dev/null; then \
      apt-get update && apt-get install -y --no-install-recommends git procps ca-certificates \
      && rm -rf /var/lib/apt/lists/*; \
    else \
      conda install -y -c conda-forge git procps-ng && conda clean --all --yes; \
    fi

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
ENV PYTHONPATH="/opt/microgwas/workflow/scripts:${PYTHONPATH:-}"

ENV MPLBACKEND=Agg
ENV MPLCONFIGDIR=/tmp/matplotlib

RUN panfeed --version \
    && panfeed-get-clusters --help > /dev/null \
    && panfeed-get-kmers --help > /dev/null \
    && panfeed-plot --help > /dev/null \
    && echo "=== All checks passed ==="

WORKDIR /data
