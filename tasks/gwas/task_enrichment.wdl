version 1.0

# task_enrichment.wdl
#
# Functional enrichment of hit gene sets against the reference background:
# Fisher's exact tests over COG categories, GO terms and KEGG pathways, with
# FDR correction, plus the bar plots of the results.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: run_enrich, run_enrich_alt,
# run_enrich_wg, run_enrichment_plots, run_enrichment_plots_alt,
# run_enrichment_plots_wg.
#
# Consolidation note: the three `run_enrich*` rules run the same script with
# the same flags, differing only in which annotated summary they read and where
# they write; likewise the three `run_enrichment_plots*` rules. Each trio
# collapses to one parameterised task here.

task run_enrich {

  meta {
    version: "0.1.0"
    description: "Test whether a hit gene set is enriched for COG categories, GO terms or KEGG pathways relative to the reference background. Covers the microGWAS `run_enrich`, `run_enrich_alt` and `run_enrich_wg` rules."
  }

  input {
    String phenotype

    # Distinguishes this enrichment from the others in the same run; used only
    # in output filenames. "" for the main unitig summary, otherwise one of
    # gpa | rare | panfeed | ridge | lasso.
    String label

    File annotated_summary
    File annotated_reference
    File go_obo

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-enrich:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # Positional output order is COG, GO, KEGG - upstream relies on Snakemake
    # expanding its `output:` block in declaration order, so it is spelled out
    # explicitly here.
    python3 "$SCRIPTS/functional_enrichment.py" ~{annotated_summary} \
      ~{annotated_reference} \
      ~{go_obo} \
      COG_~{label}.tsv \
      GO_~{label}.tsv \
      KEGG_~{label}.tsv

    for f in COG GO KEGG; do
      echo "$f: $(( $(wc -l < ${f}_~{label}.tsv) - 1 )) terms tested"
    done
  >>>

  output {
    String date = read_string("DATE")
    File cog  = "COG_~{label}.tsv"
    File go   = "GO_~{label}.tsv"
    File kegg = "KEGG_~{label}.tsv"
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

task run_enrichment_plots {

  meta {
    version: "0.1.0"
    description: "Plot COG, GO and KEGG enrichment results, in both raster and vector formats. Covers the microGWAS `run_enrichment_plots`, `_alt` and `_wg` rules."
  }

  input {
    String phenotype
    String label

    File cog
    File go
    File kegg

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # Upstream scatters over {format}; both formats are produced here in one
    # task since the script is cheap and re-reading the inputs is trivial
    for fmt in png svg; do
      python3 "$SCRIPTS/enrich_plots.py" \
        ~{cog}  "COG_~{label}.$fmt" \
        ~{go}   "GO_~{label}.$fmt" \
        ~{kegg} "KEGG_~{label}.$fmt"
    done
  >>>

  output {
    String date = read_string("DATE")
    File cog_png  = "COG_~{label}.png"
    File cog_svg  = "COG_~{label}.svg"
    File go_png   = "GO_~{label}.png"
    File go_svg   = "GO_~{label}.svg"
    File kegg_png = "KEGG_~{label}.png"
    File kegg_svg = "KEGG_~{label}.svg"
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
