# microgwas-mash
#
# Backs the microGWAS `mash_sketch` and `distance` rules
# (workflow/envs/mash.yaml, which pins `mash=2.1` AND `pyseer>=1.3.6`).
#
# The pyseer dependency is easy to miss: the `distance` rule is
#   mash dist -p {threads} {input} {input} | square_mash > {output}
# and `square_mash` is a pyseer console_script (pyseer.mash:main) that
# reshapes mash's long-format output into the square distance matrix every
# downstream task expects. StaPH-B's mash image alone would leave that
# pipeline stage missing.
#
# mash stays at the pinned 2.1 by building on StaPH-B's exact-version image
# rather than pulling mash from conda (where 2.1 is unavailable and only 2.3
# is published).
FROM staphb/mash:2.1

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="mash 2.1 + pyseer's square_mash helper, for microGWAS distance-matrix tasks"
LABEL version="2.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    procps \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# pyseer.mash (the module behind square_mash) imports only pandas, so install
# pyseer with --no-deps and supply that one dependency explicitly. This avoids
# dragging in glmnet_python/tensorflow-scale build dependencies for a helper
# that reshapes a text table.
ARG PYSEER_VERSION=1.4.2
RUN pip3 install --no-cache-dir --break-system-packages pandas \
    && pip3 install --no-cache-dir --break-system-packages --no-deps \
       git+https://github.com/mgalardini/pyseer.git@${PYSEER_VERSION}

RUN mash --version \
    && square_mash --help > /dev/null \
    && echo "=== All checks passed ==="

WORKDIR /data
