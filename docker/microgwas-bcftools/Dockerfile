# microgwas-bcftools
#
# Backs the microGWAS `prepare_vcf_variants` rule
# (workflow/envs/bcftools.yaml: bcftools=1.13, biopython=1.78, pandas=1.3.3,
# pysam=0.17.0).
#
# The rule merges per-sample snippy VCFs, splits multiallelics, partitions
# into rare (<5% AF) and common (>1% AF) sets, and then runs
# workflow/scripts/vcf2deleterious.py over the rare set - hence the python
# stack alongside bcftools itself.
#
# Layering note: the obvious base, staphb/bcftools:1.13, is Ubuntu 16.04
# (Xenial), shipping pip 8.1.1 with no setuptools - it cannot install
# biopython or pandas at any version. Same trap as microgwas-mash. So the
# layering is inverted: modern python base, with bcftools 1.13 built from
# samtools' own release tarball. That keeps the pin exact rather than taking
# whatever a distro happens to ship.
FROM python:3.10-slim-bookworm

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="bcftools 1.13 + biopython/pandas/pysam, for microGWAS rare and common variant preparation"
LABEL version="1.13"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    git \
    procps \
    gzip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libgsl-dev \
    libperl-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# bcftools at the exact pinned version. bgzip/tabix come with the bundled
# htslib, and vcf2deleterious.py pipes into bgzip.
ARG BCFTOOLS_VERSION=1.13
RUN wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && tar -xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && cd "bcftools-${BCFTOOLS_VERSION}" \
    && ./configure --prefix=/usr/local --enable-libgsl \
    && make -j"$(nproc)" \
    && make install \
    # bcftools' own `make install` does not install bgzip/tabix - those belong
    # to the bundled htslib, which lives INSIDE the bcftools source tree and
    # needs its own install step. vcf2deleterious.py pipes into bgzip, so a
    # missing bgzip would surface as a runtime failure inside the workflow
    # rather than here.
    && make -C "htslib-${BCFTOOLS_VERSION}" install \
    && cd / \
    && ldconfig \
    && rm -rf "/bcftools-${BCFTOOLS_VERSION}" "/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

# Versions pinned to match workflow/envs/bcftools.yaml exactly.
# numpy is pinned even though the env file does not name it: pandas 1.3.3 was
# built against the numpy 1.x C ABI, and pip will otherwise resolve numpy 2.x,
# giving "numpy.dtype size changed" on the first `import pandas`.
RUN pip install --no-cache-dir \
    "numpy==1.22.4" \
    "biopython==1.78" \
    "pandas==1.3.3" \
    "pysam==0.17.0"

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

RUN bcftools --version | head -1 \
    && test "$(bcftools --version | head -1 | awk '{print $2}')" = "${BCFTOOLS_VERSION}" \
    && bgzip --version | head -1 \
    && python3 -c "import Bio, pandas, pysam; print('python deps OK')" \
    && test -f /opt/microgwas/workflow/scripts/vcf2deleterious.py \
    && echo "=== All checks passed ==="

WORKDIR /data
