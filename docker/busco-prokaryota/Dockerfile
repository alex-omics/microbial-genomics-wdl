FROM staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08

LABEL maintainer="aarvani1"
LABEL description="BUSCO 5.7.1 with the complete OrthoDB v10 prokaryote lineage set (99 datasets) \
and placement files pre-loaded for fully offline Terra execution, including --auto-lineage-prok."

# Single consolidated download path. The base image ships bacteria_odb10 and
# file_versions.tsv under /busco_downloads; earlier revisions of this image
# added lineages under /data/busco_downloads instead, which left two disjoint
# trees and made auto-lineage impossible (no manifest at the path the task
# pointed at). Everything now lives under /busco_downloads.
ENV BUSCO_DOWNLOAD_PATH=/busco_downloads

COPY fetch_busco_data.sh /usr/local/bin/fetch_busco_data.sh

RUN chmod +x /usr/local/bin/fetch_busco_data.sh \
    && /usr/local/bin/fetch_busco_data.sh "${BUSCO_DOWNLOAD_PATH}" \
    && ls "${BUSCO_DOWNLOAD_PATH}/lineages" | wc -l
