# microgwas-mlst
#
# Backs the microGWAS `lineage_st` rule (workflow/envs/mlst.yaml, which pins
# `mlst=2.16` and also requires GNU `parallel`).
#
# The rule builds a job file of per-genome `mlst` invocations, fans them out
# with `parallel -j {threads}`, concatenates the results, then pipes them
# through workflow/scripts/sanitize_mlst.py - so the image needs mlst,
# parallel, AND python, none of which come together in a stock image.
#
# Version note: bioconda has no mlst 2.16 build (it jumps 2.9 -> 2.33.1), so
# this builds on StaPH-B's 2.16.2, the patch release of the pinned 2.16 line.
FROM staphb/mlst:2.16.2

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="mlst 2.16.2 + GNU parallel + microGWAS sanitize_mlst.py, for lineage assignment"
LABEL version="2.16.2"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    parallel \
    python3 \
    python3-pip \
    procps \
    gzip \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# `parallel` prints a citation notice to stderr on first run unless silenced,
# which pollutes Cromwell task logs.
RUN mkdir -p /root/.parallel && touch /root/.parallel/will-cite

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

RUN mlst --version \
    && parallel --version | head -1 \
    && python3 --version \
    && test -f /opt/microgwas/workflow/scripts/sanitize_mlst.py \
    && echo "=== All checks passed ==="

WORKDIR /data
