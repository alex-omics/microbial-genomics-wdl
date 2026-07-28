version 1.0

# task_annotation.wdl
#
# Functional annotation of the pangenome with eggNOG-mapper, which supplies the
# COG/GO/KEGG terms every enrichment and annotated-summary task depends on.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: sample_whole_pangenome,
# download_eggnog, annotate_pangenome, annotate_reference.
#
# COST WARNING: the eggNOG database is tens of GB. `download_eggnog` is
# deliberately a separate task producing a tarball output, so it can be run
# once, parked in a GCS bucket, and passed back in as `eggnog_db_tarball` on
# every subsequent run. Do not let it re-download per iteration.

task download_eggnog {

  meta {
    version: "0.1.0"
    description: "Download the eggNOG annotation database once, for reuse across runs. Translated from the microGWAS `download_eggnog` rule."
  }

  input {
    # config["eggnogdb"]: NCBI tax ID of the HMM database.
    # "2" = all Bacteria (the production default).
    # "1236" = Gammaproteobacteria, which is what upstream's own test harness
    # substitutes to keep the download small.
    String eggnog_taxid = "2"

    Int memory    = 16
    Int cpu       = 4
    Int disk_size = 200
    String docker_image = "aarvani1/microgwas-eggnog:2.1.13"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    emapper.py --version 2>&1 | tee VERSION

    # Patches the stale database URL baked into eggNOG-mapper releases.
    # Without this the download below 404s.
    eggnog-mapper-fixurl

    mkdir -p eggnog-mapper
    download_eggnog_data.py --data_dir eggnog-mapper -y -H -d ~{eggnog_taxid}

    test -f eggnog-mapper/eggnog.db

    # Tarred so the whole directory can round-trip through a single File output
    tar -czf eggnog-mapper.tar.gz eggnog-mapper
    du -sh eggnog-mapper.tar.gz | tee DB_SIZE
  >>>

  output {
    String date            = read_string("DATE")
    String emapper_version = read_string("VERSION")
    String db_size         = read_string("DB_SIZE")
    File eggnog_db_tarball = "eggnog-mapper.tar.gz"
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

task sample_whole_pangenome {

  meta {
    version: "0.1.0"
    description: "Pick one representative protein sequence per gene cluster, preferring the focus reference strains, to give eggNOG-mapper a compact input. Translated from the microGWAS `sample_whole_pangenome` rule."
  }

  input {
    File gene_presence_absence_csv
    File gene_data

    # config["annotation_references"], the raw flag string, e.g.
    # "--focus-strain 536 --focus-strain CFT073 --focus-strain IAI39 ..."
    String annotation_references

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    python3 "$SCRIPTS/sample_pangenome.py" \
      ~{gene_presence_absence_csv} \
      ~{gene_data} \
      ~{annotation_references} \
      > pangenome_sample.faa

    grep -c '^>' pangenome_sample.faa | tee N_PROTEINS
  >>>

  output {
    String date    = read_string("DATE")
    Int n_proteins = read_int("N_PROTEINS")
    File pangenome_faa = "pangenome_sample.faa"
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

task annotate_pangenome {

  meta {
    version: "0.1.0"
    description: "Annotate the sampled pangenome proteins with eggNOG-mapper, producing the COG/GO/KEGG assignments used downstream. Translated from the microGWAS `annotate_pangenome` rule."
  }

  input {
    File pangenome_faa
    File eggnog_db_tarball

    Int memory    = 64
    Int cpu       = 16
    Int disk_size = 300
    String docker_image = "aarvani1/microgwas-eggnog:2.1.13"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    emapper.py --version 2>&1 | tee VERSION

    mkdir -p eggnog_data
    tar -xzf ~{eggnog_db_tarball} -C eggnog_data --strip-components=1

    # Upstream fails the rule outright if the database directory is missing,
    # rather than letting emapper.py produce a confusing partial result
    if [[ ! -f eggnog_data/eggnog.db ]]; then
      echo "ERROR: eggnog.db not found in the supplied database tarball" >&2
      ls -R eggnog_data >&2
      exit 1
    fi

    emapper.py -i ~{pangenome_faa} -o pangenome \
      --cpu ~{cpu} \
      --target_orthologs one2one \
      --go_evidence all \
      --tax_scope Bacteria \
      --pfam_realign none \
      --override \
      --data_dir eggnog_data

    wc -l < pangenome.emapper.annotations | tee N_ANNOTATIONS
  >>>

  output {
    String date            = read_string("DATE")
    String emapper_version = read_string("VERSION")
    String eggnog_docker_image = docker_image
    Int n_annotation_lines = read_int("N_ANNOTATIONS")
    File annotations       = "pangenome.emapper.annotations"
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

task annotate_reference {

  meta {
    version: "0.1.0"
    description: "Build the reference-strain annotation background that functional enrichment tests are scored against. Translated from the microGWAS `annotate_reference` rule."
  }

  input {
    File gene_presence_absence_csv
    File annotations

    # config["enrichment_reference"], a single strain name e.g. "IAI39"
    String enrichment_reference

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # Restrict the pangenome annotations to gene clusters present in the
    # enrichment reference - this is the universe the Fisher tests draw from
    python3 "$SCRIPTS/subset_annotations.py" ~{annotations} \
      --pangenome ~{gene_presence_absence_csv} \
      --focus-strain ~{enrichment_reference} \
      --only-focus \
      > reference.emapper.annotations

    # --no-summary: emit only the annotation columns, with no association
    # statistics joined on (there is no summary table for the background set)
    python3 "$SCRIPTS/enhance_summary.py" /dev/null reference.emapper.annotations \
      --no-summary \
      > annotated_reference.tsv

    wc -l < annotated_reference.tsv | tee N_REFERENCE_GENES
  >>>

  output {
    String date = read_string("DATE")
    Int n_reference_genes = read_int("N_REFERENCE_GENES")
    File annotated_reference = "annotated_reference.tsv"
    File reference_annotations = "reference.emapper.annotations"
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
