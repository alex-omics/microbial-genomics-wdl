# pyseer
#
# Standalone pyseer image for the pyseer_gwas workflow: LMM association
# (--lrt / --lineage) over a Panaroo (or other) presence/absence Rtab, plus a
# phylogeny-derived kinship matrix. Deliberately lighter than
# microgwas-pyseer, which backs the much larger multimodal_gwas translation
# (unitigs, panfeed, snippy, eggNOG, whole-genome elastic net) and carries
# that pipeline's dependencies along with it. The two images are pinned to
# the same pyseer version but serve different workflows and are expected to
# drift independently as each pipeline evolves.
#
# Version pins below (python 3.9, numpy<1.23, scipy==1.7.0) are reused as-is
# from microgwas-pyseer's Dockerfile, not re-derived: pyseer 1.4.x and its
# glmnet_python fork are built against the numpy 1.x C ABI, and an unpinned
# numpy resolves to 2.x and fails at import with "numpy.dtype size changed".
# scipy 1.7.0 has no cp310 wheel, which is why python is 3.9 and not 3.10.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="pyseer v1.4.2 - standalone microbial pangenome-wide association studies (GWAS)"
LABEL version="1.4.2"
LABEL upstream="https://github.com/mgalardini/pyseer"

# Run everything as root for Terra compatibility
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    gzip \
    wget \
    git \
    gfortran \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# pyseer's well-behaved dependencies via conda (everything except
# glmnet_py, which is broken in the conda channels). pandas is required by
# the annotation-join step in tasks/pyseer.wdl's pyseer_annotate_results.
#
# matplotlib-base, pybedtools and tqdm are NOT microGWAS-script extras -
# pyseer's own setup.py lists matplotlib, pybedtools and tqdm in
# install_requires. Installing pybedtools via conda (which brings its own
# compiled C extension) is what lets the pip install below skip building it
# from source; without it, pip tries to compile pybedtools itself and fails
# on a missing g++. bedtools is pybedtools' own runtime dependency. bwa,
# mash, bedops and seaborn genuinely are microGWAS-script-only and stay out.
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.9 \
    "numpy<1.23.0" \
    "scipy==1.7.0" \
    pandas \
    scikit-learn \
    statsmodels \
    pysam \
    dendropy \
    matplotlib-base \
    pybedtools \
    bedtools \
    tqdm \
    pip \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# Install pyseer at the pinned tag - it pulls its own glmnet_python fork
ARG PYSEER_VERSION=1.4.2
RUN pip install --no-cache-dir \
    git+https://github.com/mgalardini/pyseer.git@${PYSEER_VERSION}

# pyseer's scripts/ helpers aren't installed by pip. Only phylogeny_distance.py
# (kinship matrix) and count_patterns.py (pattern-counted significance
# threshold, used in place of hand-rolled Bonferroni) are needed here.
RUN git clone --depth 1 --branch ${PYSEER_VERSION} \
      https://github.com/mgalardini/pyseer.git /opt/pyseer-scripts \
    && chmod +x /opt/pyseer-scripts/scripts/*.py \
    && ln -sf /opt/pyseer-scripts/scripts/phylogeny_distance.py /opt/conda/bin/ \
    && ln -sf /opt/pyseer-scripts/scripts/count_patterns.py /opt/conda/bin/

# Verify installation - fail the build if anything is broken. pybedtools and
# matplotlib are checked explicitly (not just implied by `import pyseer`)
# because they are pyseer's own install_requires, not microGWAS extras - the
# first cut of this image dropped them as apparently GWAS-DAG-only and only
# failed at the pip install step above, not here.
RUN pyseer --version \
    && python -c "import pyseer; print('pyseer imported successfully')" \
    && python -c "import glmnet_python; print('glmnet_py imported successfully')" \
    && python -c "import pandas; print('pandas imported successfully')" \
    && python -c "import pybedtools; print('pybedtools imported successfully')" \
    && python -c "import matplotlib; print('matplotlib imported successfully')" \
    && phylogeny_distance.py --help > /dev/null \
    && count_patterns.py --help > /dev/null \
    && echo "=== All checks passed ==="

WORKDIR /data
