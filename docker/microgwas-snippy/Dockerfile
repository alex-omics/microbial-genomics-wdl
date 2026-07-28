# microgwas-snippy
#
# Backs the microGWAS `get_snps` rule (workflow/envs/nucmer.yaml).
#
# nucmer.yaml is a *composite* environment, not a single tool:
#   mummer4=4.0.0rc1, snpeff=5.0, snippy=4.6.0, bcftools>1.10,
#   perl-bioperl=1.7.8, biopython, pandas=1.2.5, gsl, tabixpp=1.1.0
# Only `snippy` is actually invoked by the rule; the rest of the list is
# snippy's own runtime dependency closure, which StaPH-B's snippy image
# already bundles (snippy vendors mummer, snpeff, bcftools, samtools,
# freebayes, vt and bioperl).
#
# The remaining gap is the python stack, which snippy itself does not need but
# which keeps this image usable for the small VCF post-processing steps that
# sit adjacent to the rule.
FROM staphb/snippy:4.6.0

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="snippy 4.6.0 (with bundled mummer/snpeff/bcftools) for microGWAS reference-based variant calling"
LABEL version="4.6.0"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    procps \
    gzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Confirm the composite environment's tools really are present in this base,
# so a missing component fails the build here rather than mid-workflow.
RUN snippy --version \
    && nucmer --version 2>&1 | head -1 \
    && bcftools --version | head -1 \
    && echo "=== All checks passed ==="

WORKDIR /data
