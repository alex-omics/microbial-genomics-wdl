FROM staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08

LABEL maintainer="aarvani1"
LABEL description="BUSCO 5.7.1 with spirochaetales_odb10 and spirochaetia_odb10 pre-loaded for offline Terra execution."

COPY spirochaetales_odb10.2024-01-08.tar.gz /tmp/
COPY spirochaetia_odb10.2024-01-08.tar.gz /tmp/

RUN mkdir -p /data/busco_downloads/lineages \
    && tar -xzf /tmp/spirochaetales_odb10.2024-01-08.tar.gz \
       -C /data/busco_downloads/lineages/ \
    && tar -xzf /tmp/spirochaetia_odb10.2024-01-08.tar.gz \
       -C /data/busco_downloads/lineages/ \
    && rm /tmp/*.tar.gz \
    && echo "Installed lineages:" \
    && ls /data/busco_downloads/lineages/