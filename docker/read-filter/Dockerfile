# read-filter
#
# Image for the filter_reads workflow: samtools (filters a BAM, record for record,
# so MM/ML modification tags survive), chopper (filters and crops FASTQ on read
# length and mean quality) and seqkit (before/after read statistics) in one
# container. No single public image carries all three.
#
# Build for Terra's architecture on any host:
#   docker build --platform linux/amd64 -t aarvani1/read-filter:1.0.0 docker/read-filter
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="samtools, chopper and seqkit for filtering ONT reads by length and quality"
LABEL version="1.0.0"

# Run everything as root for Terra compatibility
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# pigz compresses the filtered FASTQ in parallel; gzip alone would make the
# compression step the slowest part of filtering a multi-Gb run.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    samtools=1.24 \
    chopper=0.14.1 \
    seqkit=2.14.0 \
    pigz \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# Fail the build if a tool is missing or is not the version asked for
RUN samtools --version | head -n1 | grep -q "samtools 1.24" \
    && chopper --version | grep -q "0.14.1" \
    && seqkit version | grep -q "2.14.0" \
    && pigz --version \
    && echo "=== All checks passed ==="

WORKDIR /data
