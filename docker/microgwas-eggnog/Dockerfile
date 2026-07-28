# microgwas-eggnog
#
# Backs the microGWAS `download_eggnog` and `annotate_pangenome` rules
# (workflow/envs/eggnog-mapper.yaml: eggnog-mapper>=2.1.6, diamond>=2.0.11,
# plus a pip-installed `eggnog-mapper-fixurl` patch).
#
# That fixurl patch is not cosmetic - the `download_eggnog` rule calls
# `eggnog-mapper-fixurl` before `download_eggnog_data.py`, because upstream
# eggNOG-mapper ships a stale database download URL. Any base image, however
# recent, still needs this layer on top.
#
# Note on operating cost: the eggNOG database is tens of GB. It is downloaded
# once by the download_eggnog task and should be persisted to a GCS bucket and
# passed back in as `eggnog_db_tarball`, NOT re-downloaded per run.
FROM quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2

LABEL maintainer="Alex Arvanitis <aarvanitis1@mgh.harvard.edu>"
LABEL description="eggNOG-mapper 2.1.13 + eggnog-mapper-fixurl patch, for microGWAS pangenome functional annotation"
LABEL version="2.1.13"
LABEL upstream="https://github.com/microbial-pangenomes-lab/microGWAS"

USER root

RUN if command -v apt-get > /dev/null; then \
      apt-get update && apt-get install -y --no-install-recommends procps wget ca-certificates \
      && rm -rf /var/lib/apt/lists/*; \
    fi

# The URL-fixing patch required by the download_eggnog rule
RUN pip install --no-cache-dir eggnog-mapper-fixurl

RUN emapper.py --version \
    && download_eggnog_data.py --help > /dev/null \
    && command -v eggnog-mapper-fixurl \
    && diamond --version \
    && echo "=== All checks passed ==="

WORKDIR /data
