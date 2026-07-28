version 1.0

# task_mapback.wdl
#
# Mapping significant variants back onto genomes and reference coordinates,
# and the Manhattan plot built from that mapping.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: run_map_back, run_map_back_all,
# run_map_back_wg, run_manhattan_plot.
#
# map_back.py accepts either a file-of-files or a directory for its genome and
# GFF arguments. These tasks stage symlink directories rather than path lists:
# under Cromwell that is the more robust of the two, since directory mode looks
# up `{strain}.gff` by name instead of string-matching basenames across two
# independently-localized path lists.

task run_map_back {

  meta {
    version: "0.1.0"
    description: "Map significant unitigs back onto the study genomes with bwa fastmap and annotate hits with the gene they fall in. Translated from the microGWAS `run_map_back` rule."
  }

  input {
    String phenotype

    File filtered_variants          # unitigs_filtered.tsv
    Array[File] sample_fastas
    Array[File] reference_fastas = []
    Array[File] sample_gffs
    Array[File] reference_gffs = []
    File gene_presence_absence_csv

    Int memory    = 32
    Int cpu       = 4
    Int disk_size = 200
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # Stage genomes and annotations as name-keyed directories
    mkdir -p genomes annotations tmp_mapback
    for f in ~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas}; do
      ln -s "$f" "genomes/$(basename "$f")"
    done
    for f in ~{sep=' ' sample_gffs} ~{sep=' ' reference_gffs}; do
      ln -s "$f" "annotations/$(basename "$f")"
    done

    # map_back.py exits 0 without output when the filtered table is empty
    # (no significant hits), so the header is written first and unconditionally
    printf 'strain\tunitig\tcontig\tstart\tend\tstrand\tupstream\tgene\tdownstream\n' > mapped.tsv

    python3 "$SCRIPTS/map_back.py" ~{filtered_variants} genomes \
      --tmp-prefix tmp_mapback \
      --gff annotations \
      --print-details \
      --pangenome ~{gene_presence_absence_csv} \
      >> mapped.tsv

    rm -rf tmp_mapback

    echo $(( $(wc -l < mapped.tsv) - 1 )) | tee N_MAPPED
  >>>

  output {
    String date  = read_string("DATE")
    Int n_mapped = read_int("N_MAPPED")
    File mapped  = "mapped.tsv"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB" # TES
    preemptible: 0
    maxRetries:  2
  }
}

task run_map_back_all {

  meta {
    version: "0.1.0"
    description: "Map ALL tested unitigs onto the reference genomes and collapse them into per-position skyline coordinates for plotting. Translated from the microGWAS `run_map_back_all` rule."
  }

  input {
    String phenotype

    File all_variants               # unitigs.tsv (unfiltered)
    Array[File] reference_fastas
    Array[File] reference_gffs

    Int memory    = 64
    Int cpu       = 8
    Int disk_size = 300
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # This one maps against references only - the Manhattan plot needs a single
    # consistent coordinate system, which the study genomes cannot provide
    mkdir -p ref_fastas ref_gffs tmp_mapback
    for f in ~{sep=' ' reference_fastas}; do ln -s "$f" "ref_fastas/$(basename "$f")"; done
    for f in ~{sep=' ' reference_gffs};   do ln -s "$f" "ref_gffs/$(basename "$f")"; done

    printf 'strain\tunitig\tcontig\tstart\tend\tstrand\tupstream\tgene\tdownstream\n' > mapped_raw.tsv

    python3 "$SCRIPTS/map_back.py" ~{all_variants} ref_fastas \
      --tmp-prefix tmp_mapback \
      --gff ref_gffs \
      --print-details \
      >> mapped_raw.tsv

    # make_skyline joins the mapping back onto the association p-values, which
    # is what the Manhattan plot actually consumes
    python3 "$SCRIPTS/make_skyline.py" ~{all_variants} mapped_raw.tsv > mapped_all.tsv

    rm -rf tmp_mapback

    echo $(( $(wc -l < mapped_all.tsv) - 1 )) | tee N_MAPPED_ALL
  >>>

  output {
    String date      = read_string("DATE")
    Int n_mapped_all = read_int("N_MAPPED_ALL")
    File mapped_all  = "mapped_all.tsv"
    File mapped_raw  = "mapped_raw.tsv"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB" # TES
    preemptible: 0
    maxRetries:  2
  }
}

task run_map_back_wg {

  meta {
    version: "0.1.0"
    description: "Map the unitigs selected by a whole-genome elastic net model back onto the study genomes. Translated from the microGWAS `run_map_back_wg` rule."
  }

  input {
    String phenotype
    String model                    # "ridge" or "lasso"

    File model_variants             # <model>.tsv
    Array[File] sample_fastas
    Array[File] reference_fastas = []
    Array[File] sample_gffs
    Array[File] reference_gffs = []
    File gene_presence_absence_csv

    Int memory    = 32
    Int cpu       = 4
    Int disk_size = 200
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    mkdir -p genomes annotations tmp_mapback
    for f in ~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas}; do
      ln -s "$f" "genomes/$(basename "$f")"
    done
    for f in ~{sep=' ' sample_gffs} ~{sep=' ' reference_gffs}; do
      ln -s "$f" "annotations/$(basename "$f")"
    done

    printf 'strain\tunitig\tcontig\tstart\tend\tstrand\tupstream\tgene\tdownstream\n' > mapped_~{model}.tsv

    python3 "$SCRIPTS/map_back.py" ~{model_variants} genomes \
      --tmp-prefix tmp_mapback \
      --gff annotations \
      --print-details \
      --pangenome ~{gene_presence_absence_csv} \
      >> mapped_~{model}.tsv

    rm -rf tmp_mapback
  >>>

  output {
    String date = read_string("DATE")
    File mapped = "mapped_~{model}.tsv"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB" # TES
    preemptible: 0
    maxRetries:  2
  }
}

task run_manhattan_plot {

  meta {
    version: "0.1.0"
    description: "Draw the Manhattan plot of unitig association p-values along the reference genome, with the pattern-counted significance line. Translated from the microGWAS `run_manhattan_plot` rule."
  }

  input {
    String phenotype

    File mapped_all
    File unitigs_patterns

    # config["enrichment_reference"] - the reference strain whose coordinate
    # system the plot is drawn in, e.g. "IAI39"
    String reference_strain

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    THRESHOLD=$(python3 "$SCRIPTS/count_patterns.py" --threshold ~{unitigs_patterns})
    echo "$THRESHOLD" | tee THRESHOLD

    # Upstream scatters over {format} to produce both raster and vector output
    for fmt in png svg; do
      python3 "$SCRIPTS/manhattan_plot.py" ~{mapped_all} ~{reference_strain} \
        "manhattan.$fmt" \
        -t "$THRESHOLD"
    done
  >>>

  output {
    String date = read_string("DATE")
    String threshold = read_string("THRESHOLD")
    File manhattan_png = "manhattan.png"
    File manhattan_svg = "manhattan.svg"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " HDD"
    disk:        disk_size + " GB" # TES
    preemptible: 1
    maxRetries:  2
  }
}
