# listt
#
# Image for the call_ospA_type workflow: Pfizer's LISTT (Lyme In Silico Typing
# Tool, https://github.com/pfizer-opensource/LISTT) baked in at a pinned commit,
# with the tool versions from LISTT's own environment.yml.
#
# LISTT is not an installable package: ospa_typing.py reads ref/, variants/ and
# alleles/ relative to the current directory, so the repository is cloned to
# /opt/LISTT and tasks run from there, passing absolute input paths with -p.
# The upstream conda export (environment.yml) is a linux-64 lockfile with ~80
# Perl packages pulled in by entrez-direct; only the tools LISTT actually
# calls are installed here at the same versions.
#
# LISTT's own code is left unmodified. The wrapper logic (staging, thresholds,
# reconciliation) lives in tasks/listt.wdl.
#
# Build for Terra's architecture on any host:
#   docker build --platform linux/amd64 -t aarvani1/listt:1.2.0 docker/listt
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="LISTT v1.2.0 - Lyme In Silico Typing Tool (OspA IST) for Borrelia burgdorferi s.l."
LABEL version="1.2.0"
LABEL upstream="https://github.com/pfizer-opensource/LISTT"
LABEL license="Apache-2.0 (LISTT)"

# Run everything as root for Terra compatibility
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    gzip \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Versions match LISTT's environment.yml (bwa 0.7.17, samtools/bcftools/htslib
# 1.18, blast 2.15.0, mafft 7.520, python 3.11, biopython 1.81, numpy 1.26,
# pandas 2.1). The paper (Lee et al. 2024) was run with bcftools 1.17 and mafft
# 7.480; the repository's environment is the later of the two, so it is used.
# htslib provides bgzip, which build_consensus.py calls directly.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    bwa=0.7.17 \
    samtools=1.18 \
    bcftools=1.18 \
    htslib=1.18 \
    blast=2.15.0 \
    mafft=7.520 \
    biopython=1.81 \
    "numpy=1.26" \
    "pandas=2.1" \
    screed=1.1.2 \
    pip \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# pyvcf 0.6.8 as pinned by LISTT's environment.yml (imported as `vcf`). Its
# setup.py uses `use_2to3`, which setuptools removed in v58, so it only builds
# against an older setuptools and without pip's isolated build environment.
RUN pip install --no-cache-dir "setuptools<58" wheel \
    && pip install --no-cache-dir --no-build-isolation pyvcf==0.6.8

# Pin LISTT by commit, not tag: tags can be moved. b9c6569 is tag 1.2.0
# (Merge of PR #2, scheme-update-1.2). The build fails if the checkout differs.
ARG LISTT_REF=1.2.0
ARG LISTT_COMMIT=b9c6569
RUN git clone --branch ${LISTT_REF} https://github.com/pfizer-opensource/LISTT.git /opt/LISTT \
    && cd /opt/LISTT \
    && git rev-parse HEAD | grep -q "^${LISTT_COMMIT}" \
    && git rev-parse HEAD > /opt/LISTT/COMMIT \
    && rm -rf /opt/LISTT/.git

# Fail the build if anything LISTT shells out to is missing or cannot import
RUN { bwa 2>&1 || true; } | grep -q "Version: 0.7.17" \
    && samtools --version | head -n1 \
    && bcftools --version | head -n1 \
    && bgzip --version | head -n1 \
    && blastn -version | head -n1 \
    && mafft --version 2>&1 | head -n1 \
    && cd /opt/LISTT \
    && python -c "import vcf, Bio, screed, numpy, pandas; print('python imports ok')" \
    && python ospa_typing.py --help > /dev/null \
    && echo "=== All checks passed ==="

# Audit of LISTT's reference tables; see its docstring
COPY audit_listt_tables.py /opt/audit/audit_listt_tables.py

ENV LISTT_VERSION=1.2.0

WORKDIR /opt/LISTT
