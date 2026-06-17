FROM staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08

LABEL maintainer="aarvani1"
LABEL description="BUSCO 5.7.1 with all prokaryote lineages (bacteria & archaeae) \ 
pre-downloaded for offline execution in Terra."

ENV BUSCO_DOWNLOADS=/busco_downloads

RUN busco --download prokaryota \
    --download_path ${BUSCO_DOWNLOADS} \
    && echo "Prokaryote lineages downloaded successfully!" \
    && ls ${BUSCO_DOWNLOADS}/lineages/