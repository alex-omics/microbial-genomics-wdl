version 1.0

# task_variants.wdl
#
# Variant/feature generation for multimodal_gwas: the tasks that turn
# assemblies into the genetic elements that get tested for association.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: unitigs, pangenome, panfeed_kmers,
# ggcaller.
#
# Every task here rebuilds its own path manifest from its Array[File] inputs -
# see the header of task_setup.wdl for why.

task unitigs {

  meta {
    version: "0.1.0"
    description: "Build the de Bruijn graph unitig presence/absence matrix across all genomes with unitig-counter. Translated from the microGWAS `unitigs` rule."
  }

  input {
    # Sample assemblies followed by reference assemblies, matching the order
    # bootstrap.sh writes into out/unitigs_input.tsv (samples first, then the
    # references appended). Reference genomes are deliberately included: they
    # participate in unitig counting upstream.
    Array[File] sample_fastas
    Array[File] reference_fastas = []

    Int memory    = 32
    Int cpu       = 16
    Int disk_size = 250
    String docker_image = "quay.io/biocontainers/unitig-counter:1.1.0--h56fc30b_0"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    # Rebuild out/unitigs_input.tsv: header, then "ID<TAB>path" per genome.
    # ID is the basename minus extension, matching aid_bootstrap.py's
    # SAMPLE.fasta convention.
    printf 'ID\tPath\n' > unitigs_input.tsv
    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    for f in "${fastas[@]}"; do
      id=$(basename "$f")
      id="${id%.*}"
      printf '%s\t%s\n' "$id" "$f" >> unitigs_input.tsv
    done

    echo "Genomes for unitig counting: $(( $(wc -l < unitigs_input.tsv) - 1 ))"

    rm -rf unitigs_out
    unitig-counter \
      -strains unitigs_input.tsv \
      -output unitigs_out \
      -gzip \
      -nb-cores ~{cpu}

    # Upstream verifies the gzip stream rather than trusting the exit code,
    # because unitig-counter has been seen to emit a truncated archive
    gzip -t unitigs_out/unitigs.txt.gz

    zcat unitigs_out/unitigs.txt.gz | wc -l | tee N_UNITIGS
  >>>

  output {
    String date       = read_string("DATE")
    String unitigs_docker_image = docker_image
    Int n_unitigs     = read_int("N_UNITIGS")

    File unitigs      = "unitigs_out/unitigs.txt.gz"
    File unitigs_rtab = "unitigs_out/unitigs.unique_rows.Rtab.gz"
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

task pangenome {

  meta {
    version: "0.1.0"
    description: "Build the pangenome and core genome alignment with Panaroo in strict clean mode. Translated from the microGWAS `pangenome` rule."
  }

  input {
    # Sample GFFs followed by reference GFFs, matching out/panaroo_input.txt
    Array[File] sample_gffs
    Array[File] reference_gffs = []

    # Upstream hardcodes these in the rule's shell command rather than exposing
    # them as config keys; surfaced here as inputs with the upstream defaults.
    String clean_mode    = "strict"
    String alignment     = "core"

    Int memory    = 128
    Int cpu       = 24
    Int disk_size = 500
    String docker_image = "quay.io/biocontainers/panaroo:1.5.0--pyhdfd78af_0"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    panaroo --version 2>&1 | tee PANAROO_VERSION_RAW
    sed 's/panaroo //' PANAROO_VERSION_RAW | tee VERSION

    # Rebuild out/panaroo_input.txt: one GFF path per line
    gffs=(~{sep=' ' sample_gffs} ~{sep=' ' reference_gffs})
    printf '%s\n' "${gffs[@]}" > panaroo_input.txt
    echo "GFFs for pangenome: $(wc -l < panaroo_input.txt)"

    panaroo \
      -t ~{cpu} \
      -i panaroo_input.txt \
      -o panaroo_tmp \
      --clean-mode ~{clean_mode} \
      -a ~{alignment}

    # Upstream moves only the five outputs the rest of the DAG consumes out of
    # the temp directory and deletes the rest; kept as-is so downstream file
    # expectations match.
    mkdir -p panaroo
    mv panaroo_tmp/gene_presence_absence.Rtab   panaroo/
    mv panaroo_tmp/gene_presence_absence.csv    panaroo/
    mv panaroo_tmp/gene_data.csv                panaroo/
    mv panaroo_tmp/struct_presence_absence.Rtab panaroo/
    mv panaroo_tmp/core_gene_alignment.aln      panaroo/

    # Keep the summary for QC even though upstream discards it
    if [[ -f panaroo_tmp/summary_statistics.txt ]]; then
      mv panaroo_tmp/summary_statistics.txt panaroo/
    fi

    rm -rf panaroo_tmp

    grep -c '^>' panaroo/core_gene_alignment.aln | tee N_CORE_GENOMES
  >>>

  output {
    String date            = read_string("DATE")
    String panaroo_version = read_string("VERSION")
    String panaroo_docker_image = docker_image
    Int n_genomes_in_core  = read_int("N_CORE_GENOMES")

    File gene_presence_absence_rtab = "panaroo/gene_presence_absence.Rtab"
    File gene_presence_absence_csv  = "panaroo/gene_presence_absence.csv"
    File gene_data                  = "panaroo/gene_data.csv"
    File struct_presence_absence    = "panaroo/struct_presence_absence.Rtab"
    File core_genome_aln            = "panaroo/core_gene_alignment.aln"
    File? summary_statistics        = "panaroo/summary_statistics.txt"
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

task ggcaller {

  meta {
    version: "0.1.0"
    description: "Call genes across all assemblies with ggCaller, producing the per-sample GFF annotations Panaroo and panfeed consume. Translated from the microGWAS `ggcaller` rule."
  }

  input {
    Array[File] sample_fastas

    Int memory    = 64
    Int cpu       = 16
    Int disk_size = 300
    String docker_image = "quay.io/biocontainers/ggcaller:1.5.1--py311h9c8ac4a_0"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    # Rebuild out/ggcaller_input.txt: one assembly path per line
    fastas=(~{sep=' ' sample_fastas})
    printf '%s\n' "${fastas[@]}" > ggcaller_input.txt

    rm -rf ggcaller_out
    ggcaller \
      --refs ggcaller_input.txt \
      --out ggcaller_out \
      --threads ~{cpu}

    ls -1 ggcaller_out/GFF/*.gff | wc -l | tee N_GFFS
  >>>

  output {
    String date  = read_string("DATE")
    String ggcaller_docker_image = docker_image
    Int n_gffs   = read_int("N_GFFS")

    Array[File] sample_gffs = glob("ggcaller_out/GFF/*.gff")
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

task panfeed_kmers {

  meta {
    version: "0.1.0"
    description: "First-pass panfeed run: generate gene-cluster-specific k-mer presence/absence patterns across all gene clusters. Translated from the microGWAS `panfeed_kmers` rule."
  }

  input {
    # panfeed needs matched GFF and FASTA sets for the *samples plus
    # references*, since the pangenome it is keyed against includes both
    Array[File] sample_gffs
    Array[File] reference_gffs = []
    Array[File] sample_fastas
    Array[File] reference_fastas = []

    File gene_presence_absence_csv

    Int upstream_bp   = 250
    Int downstream_bp = 100

    Int memory    = 32
    Int cpu       = 4
    Int disk_size = 250
    String docker_image = "aarvani1/microgwas-panfeed:1.6.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    panfeed --version 2>&1 | tee VERSION

    gffs=(~{sep=' ' sample_gffs} ~{sep=' ' reference_gffs})
    printf '%s\n' "${gffs[@]}" > gff_list.txt

    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    printf '%s\n' "${fastas[@]}" > fasta_list.txt

    rm -rf panfeed_out
    panfeed \
      --gff gff_list.txt \
      --fasta fasta_list.txt \
      -o panfeed_out \
      -p ~{gene_presence_absence_csv} \
      --upstream ~{upstream_bp} --downstream ~{downstream_bp} \
      --no-filter \
      -v \
      --cores ~{cpu}

    wc -l < panfeed_out/hashes_to_patterns.tsv | tee N_PATTERNS
  >>>

  output {
    String date            = read_string("DATE")
    String panfeed_version = read_string("VERSION")
    String panfeed_docker_image = docker_image
    Int n_patterns         = read_int("N_PATTERNS")

    File panfeed_patterns   = "panfeed_out/hashes_to_patterns.tsv"
    File panfeed_conversion = "panfeed_out/kmers_to_hashes.tsv"
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
