version 1.0

# task_amr.wdl
#
# AMR and virulence gene detection across the study genomes.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rule covered: find_amr_vag.
#
# Scope note: this is a side branch. Its output feeds no association task -
# upstream calls it only as an explicit target for context alongside the GWAS
# results, and it is optional in the workflow here for the same reason.
#
# VERSION FLAG: workflow/envs/abritamr.yaml pins `abritamr=1.0.14`, but no such
# version exists on bioconda - published builds go 1.0.9 -> 1.1.0, and
# quay.io/biocontainers has no 1.0.14 tag. Per the milestone brief's
# version-tag fallback rule this is surfaced rather than silently substituted.
# The default below is the nearest published release in the 1.x line; because
# nothing downstream consumes this task's output, the mismatch cannot affect
# Tier 2 fidelity. Override `docker_image` if you need a specific build.

task find_amr_vag {

  meta {
    version: "0.1.0"
    description: "Detect AMR and virulence-associated genes across all study genomes with abritamr. Translated from the microGWAS `find_amr_vag` rule."
  }

  input {
    Array[File] sample_fastas

    # config["species_amr"], e.g. "Escherichia", "Klebsiella",
    # "Pseudomonas_aeruginosa"
    String species

    Int memory    = 64
    Int cpu       = 16
    Int disk_size = 200

    # See the VERSION FLAG note in this file's header before changing this
    String docker_image = "quay.io/biocontainers/abritamr:1.0.9--hdfd78af_0"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    abritamr --version 2>&1 | tee VERSION

    # abritamr wants a two-column "name<TAB>absolute path" manifest
    : > abritamr_input.tsv
    for f in ~{sep=' ' sample_fastas}; do
      id=$(basename "$f"); id="${id%.*}"
      printf '%s\t%s\n' "$id" "$(realpath "$f")" >> abritamr_input.tsv
    done
    echo "Genomes to screen: $(wc -l < abritamr_input.tsv)"

    mkdir -p abritamr_out
    cd abritamr_out
    abritamr run \
      --contigs ../abritamr_input.tsv \
      --species ~{species} \
      --jobs ~{cpu}
  >>>

  output {
    String date            = read_string("DATE")
    String abritamr_version = read_string("VERSION")
    String abritamr_docker_image = docker_image
    File summary_matches   = "abritamr_out/summary_matches.txt"
    File summary_partials  = "abritamr_out/summary_partials.txt"
    File summary_virulence = "abritamr_out/summary_virulence.txt"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB" # TES
    preemptible: 0
    maxRetries:  1
  }
}
