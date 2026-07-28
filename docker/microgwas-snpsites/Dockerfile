# microgwas-snpsites
#
# Backs the microGWAS `variant_sites` and `aln2vcf` rules
# (workflow/envs/snp-sites.yaml, which pins BOTH `snp-sites=2.5.1` and
# `bcftools=1.13`).
#
# `aln2vcf` runs snp-sites to emit a VCF, then immediately normalises and
# indexes it with bcftools, so a single-tool snp-sites image is not enough.
FROM staphb/snp-sites:2.5.1

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="snp-sites 2.5.1 + bcftools 1.13, for core-genome variant site extraction and VCF conversion"
LABEL version="2.5.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    procps \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libperl-dev \
    libgsl-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Build bcftools at the exact pinned 1.13 rather than taking whatever the
# distro ships, since the pin is load-bearing for Tier 2 fidelity.
ARG BCFTOOLS_VERSION=1.13
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd .. \
    && rm -rf bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2

# bgzip/tabix ship with htslib, built as part of bcftools above
RUN snp-sites -V \
    && bcftools --version | head -1 \
    && bgzip --version | head -1 \
    && echo "=== All checks passed ==="

WORKDIR /data
