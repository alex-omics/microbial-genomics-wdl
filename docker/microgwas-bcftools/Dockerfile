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
FROM staphb/bcftools:1.13

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="bcftools 1.13 + biopython/pandas/pysam, for microGWAS rare and common variant preparation"
LABEL version="1.13"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    procps \
    gzip \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Versions pinned to match workflow/envs/bcftools.yaml exactly
RUN pip3 install --no-cache-dir --break-system-packages \
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
ENV PYTHONPATH="/opt/microgwas/workflow/scripts:${PYTHONPATH}"

RUN bcftools --version | head -1 \
    && python3 -c "import Bio, pandas, pysam; print('python deps OK')" \
    && test -f /opt/microgwas/workflow/scripts/vcf2deleterious.py \
    && echo "=== All checks passed ==="

WORKDIR /data
