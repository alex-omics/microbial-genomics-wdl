# microgwas-pyseer
#
# Backs the 15 microGWAS rules using `workflow/envs/pyseer.yaml`
# (which pins `pyseer>=1.4.0`): similarity, prepare_wg, run_pyseer,
# run_pyseer_vcf, run_qq_plot, run_wg, run_wg_metrics, run_map_back,
# run_map_back_all, run_map_back_wg, run_manhattan_plot, run_panfeed,
# run_enrichment_plots{,_alt,_wg}.
#
# Bumped from the previous aarvani1/pyseer:1.4.0 to pyseer 1.4.2 (released
# upstream since that image was built). 1.4.2 still satisfies the upstream
# `>=1.4.0` pin, so this is a forward bump within the pinned range, not a
# substitution.
#
# bwa + bedtools/pybedtools are required by workflow/scripts/map_back.py,
# which shells out to `bwa fastmap` and intersects hits against GFF features.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="pyseer v1.4.2 - microbial pangenome-wide association studies (GWAS), with microGWAS helper scripts"
LABEL version="1.4.2"
LABEL upstream="https://github.com/mgalardini/pyseer"

# Run everything as root for Terra compatibility
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    gzip \
    wget \
    curl \
    git \
    gfortran \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# pyseer's well-behaved dependencies via conda
# (everything except glmnet_py, which is broken in the conda channels)
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.10 \
    numpy \
    scipy \
    pandas \
    scikit-learn \
    statsmodels \
    pysam \
    dendropy \
    matplotlib-base \
    seaborn \
    pybedtools \
    bedtools \
    bedops \
    bwa \
    mash \
    tqdm \
    pip \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# Install pyseer at the pinned tag - it pulls its own glmnet_python fork
ARG PYSEER_VERSION=1.4.2
RUN pip install --no-cache-dir \
    git+https://github.com/mgalardini/pyseer.git@${PYSEER_VERSION}

# pyseer's scripts/ helpers (phylogeny_distance.py, count_patterns.py,
# summarise_annotations.py, ...) are shipped in the repo but not installed by
# pip, so symlink them onto PATH. Note that square_mash is NOT one of these -
# it is a console_script entry point (pyseer.mash:main) already created by the
# pip install above.
RUN git clone --depth 1 --branch ${PYSEER_VERSION} \
      https://github.com/mgalardini/pyseer.git /opt/pyseer-scripts \
    && chmod +x /opt/pyseer-scripts/scripts/*.py \
    && for s in /opt/pyseer-scripts/scripts/*.py; do \
         ln -sf "$s" /opt/conda/bin/; \
       done

# Vendor microGWAS at a pinned commit for its workflow/scripts helpers
# (count_patterns.py, map_back.py, manhattan_plot.py, qq_plot.py,
# prepare_pyseer.py, prediction_metrics.py, make_skyline.py, enrich_plots.py).
# MIT-licensed, Copyright (c) 2022, Marco Galardini.
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

# Matplotlib needs a writable, headless config dir under Cromwell
ENV MPLBACKEND=Agg
ENV MPLCONFIGDIR=/tmp/matplotlib

# Verify installation - fail the build if anything is broken
RUN pyseer --version \
    && square_mash --help > /dev/null \
    && python -c "import pyseer; print('pyseer imported successfully')" \
    && python -c "import glmnet_python; print('glmnet_py imported successfully')" \
    && python -c "import pybedtools; print('pybedtools imported successfully')" \
    && bwa 2>&1 | grep -q Version \
    && test -f /opt/microgwas/workflow/scripts/count_patterns.py \
    && test -f /opt/microgwas/LICENSE \
    && echo "=== All checks passed ==="

WORKDIR /data
