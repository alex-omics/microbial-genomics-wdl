# microgwas-panfeed
#
# Backs the microGWAS `panfeed_kmers` and `panfeed_downstream` rules
# (workflow/envs/panfeed.yaml, which pins `panfeed>=1.6.1`).
#
# panfeed_downstream needs the `panfeed-get-clusters`, `panfeed-get-kmers` and
# `panfeed-plot` console scripts in addition to `panfeed` itself; all ship in
# the same distribution.
#
# Layering note: this biocontainer has NO package manager - no apt-get, no
# conda, no micromamba. It is a stripped Debian userland with a conda env
# baked in. So the microGWAS helper scripts cannot be git-cloned in place;
# they are fetched in a throwaway builder stage and copied in. This is the
# more robust pattern generally, since it places no requirements at all on the
# runtime base.

# --- Builder: exists only to fetch the pinned microGWAS source --------------
FROM debian:bookworm-slim AS microgwas-src

RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ARG MICROGWAS_COMMIT=130725016f023d78c13f08fdb961779afc2190ce
RUN git clone https://github.com/microbial-pangenomes-lab/microGWAS.git /opt/microgwas \
    && cd /opt/microgwas \
    && git checkout ${MICROGWAS_COMMIT} \
    && rm -rf .git \
    && chmod +x workflow/scripts/*.py

# --- Runtime ----------------------------------------------------------------
FROM quay.io/biocontainers/panfeed:1.6.1--pyhdfd78af_0

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="panfeed 1.6.1 + microGWAS helper scripts, for gene-cluster-specific k-mer generation"
LABEL version="1.6.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

# MIT-licensed, Copyright (c) 2022, Marco Galardini - LICENSE travels along
COPY --from=microgwas-src /opt/microgwas /opt/microgwas

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
    && python3 /opt/microgwas/workflow/scripts/count_patterns.py --help > /dev/null \
    && test -f /opt/microgwas/LICENSE \
    && echo "=== All checks passed ==="

WORKDIR /data
