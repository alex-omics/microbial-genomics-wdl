# microgwas-mash
#
# Backs the microGWAS `mash_sketch` and `distance` rules
# (workflow/envs/mash.yaml, which pins `mash=2.1` AND `pyseer>=1.3.6`).
#
# The pyseer dependency is easy to miss: the `distance` rule is
#   mash dist -p {threads} {input} {input} | square_mash > {output}
# and `square_mash` is a pyseer console_script (pyseer.mash:main) that
# reshapes mash's long-format output into the square distance matrix every
# downstream task expects. A mash-only image would leave that stage missing.
#
# Layering note: the obvious approach - FROM staphb/mash:2.1, then pip install
# pyseer - does not work. That image is Ubuntu 16.04-era, shipping Python 3.5
# and pip 8.1.1, which cannot build any remotely current pandas. So the
# layering is inverted: a modern Python base, with mash 2.1 taken from marbl's
# own release tarball. That is the canonical upstream distribution of the
# pinned version, and it keeps mash at exactly 2.1 rather than drifting to the
# 2.3 that conda would resolve.
FROM python:3.10-slim-bookworm

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="mash 2.1 + pyseer's square_mash helper, for microGWAS distance-matrix tasks"
LABEL version="2.1"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    git \
    procps \
    ca-certificates \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# mash 2.1 from marbl's official release assets - the exact pinned version
ARG MASH_VERSION=2.1
RUN wget -q "https://github.com/marbl/Mash/releases/download/v${MASH_VERSION}/mash-Linux64-v${MASH_VERSION}.tar" \
    && tar -xf "mash-Linux64-v${MASH_VERSION}.tar" \
    && mv "mash-Linux64-v${MASH_VERSION}/mash" /usr/local/bin/mash \
    && chmod +x /usr/local/bin/mash \
    && rm -rf "mash-Linux64-v${MASH_VERSION}" "mash-Linux64-v${MASH_VERSION}.tar"

# pyseer.mash (the module behind square_mash) imports only pandas, so install
# pyseer with --no-deps and supply that one dependency explicitly. This avoids
# dragging in glmnet_python and the rest of the association stack for a helper
# that reshapes a text table.
ARG PYSEER_VERSION=1.4.2
RUN pip install --no-cache-dir pandas \
    && pip install --no-cache-dir --no-deps \
       git+https://github.com/mgalardini/pyseer.git@${PYSEER_VERSION}

RUN mash --version \
    && test "$(mash --version)" = "${MASH_VERSION}" \
    && square_mash --help > /dev/null \
    && echo "=== All checks passed ==="

WORKDIR /data
