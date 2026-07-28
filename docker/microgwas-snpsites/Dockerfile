# microgwas-snpsites
#
# Backs the microGWAS `variant_sites` and `aln2vcf` rules
# (workflow/envs/snp-sites.yaml, which pins BOTH `snp-sites=2.5.1` and
# `bcftools=1.13`).
#
# `aln2vcf` runs snp-sites to emit a VCF, then immediately normalises and
# indexes it with bcftools, so a single-tool snp-sites image is not enough.
#
# Layering note: the obvious base, staphb/snp-sites:2.5.1, is Ubuntu 21.04
# ("hirsute"), which reached end of life - its archive.ubuntu.com repos now
# 404, so `apt-get update` fails outright and nothing can be layered on.
# Rather than repoint apt at old-releases.ubuntu.com and inherit an EOL
# userland, both tools are built from source on a current base. snp-sites
# publishes no release tarball, so it is built from the v2.5.1 git tag via
# autoreconf; bcftools has a release tarball at the pinned version.
FROM debian:bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="snp-sites 2.5.1 + bcftools 1.13, for core-genome variant site extraction and VCF conversion"
LABEL version="2.5.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    check \
    git \
    wget \
    procps \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libgsl-dev \
    libperl-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# snp-sites at the exact pinned version, from its git tag
ARG SNP_SITES_VERSION=2.5.1
RUN git clone --depth 1 --branch "v${SNP_SITES_VERSION}" \
      https://github.com/sanger-pathogens/snp-sites.git /tmp/snp-sites \
    && cd /tmp/snp-sites \
    && autoreconf -i -f \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf /tmp/snp-sites \
    && ldconfig

# bcftools at the exact pinned version
ARG BCFTOOLS_VERSION=1.13
RUN wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && tar -xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && cd "bcftools-${BCFTOOLS_VERSION}" \
    && ./configure --prefix=/usr/local --enable-libgsl \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf "/bcftools-${BCFTOOLS_VERSION}" "/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

RUN snp-sites -V \
    && bcftools --version | head -1 \
    && test "$(bcftools --version | head -1 | awk '{print $2}')" = "${BCFTOOLS_VERSION}" \
    && bgzip --version | head -1 \
    && echo "=== All checks passed ==="

WORKDIR /data
