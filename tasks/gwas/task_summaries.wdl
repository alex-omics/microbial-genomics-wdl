version 1.0

# task_summaries.wdl
#
# Per-gene summarisation of association hits, and joining eggNOG annotations
# onto those summaries. These are the tables you actually read to decide which
# loci came out - the Tier 2 pass criteria (HPI, aerobactin, sitABCD, zinT,
# mtfA, shiA) are judged from here and from the Manhattan plots.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: run_map_summary,
# run_map_summary_panfeed, run_map_summary_wg, run_gpa_summary,
# run_rare_summary, run_annotate_summary, run_annotate_rare_summary,
# run_annotate_gpa_summary, run_annotate_panfeed_summary,
# run_annotate_summary_wg, run_annotate_panfeed_small.
#
# Consolidation note: upstream has five separate `run_annotate_*` rules whose
# shell bodies are character-for-character identical apart from input and
# output paths. They collapse into the single parameterised `annotate_summary`
# task below rather than five near-duplicate WDL tasks. The `map_summary`
# variants genuinely differ in their flags, so those stay separate.

task map_summary {

  meta {
    version: "0.1.0"
    description: "Collapse mapped significant unitigs into a per-gene summary table, filtered by unitig length and hit counts. Translated from the microGWAS `run_map_summary` rule."
  }

  input {
    String phenotype

    File mapped
    File phenotypes
    File filtered_variants
    File gene_presence_absence_rtab
    File gene_presence_absence_csv
    Array[File] reference_gffs

    # config["summary_references"], raw flag string:
    # "--reference 536 --reference CFT073 --reference IAI39 ..."
    String summary_references

    # Spurious-hit filters. Upstream's config defaults are deliberately loose;
    # the paper's own E. coli analysis used length 30 / min_hits 9 / max_genes 10.
    Int length      = 10
    Int min_hits    = 1
    Int max_genes   = 25

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    mkdir -p ref_gffs
    for f in ~{sep=' ' reference_gffs}; do ln -s "$f" "ref_gffs/$(basename "$f")"; done

    python3 "$SCRIPTS/mapped_summary.py" ~{mapped} \
      ~{phenotypes} ~{phenotype} ~{filtered_variants} \
      --pangenome ~{gene_presence_absence_rtab} \
      --pangenome-genes ~{gene_presence_absence_csv} \
      --length ~{length} --minimum-hits ~{min_hits} --maximum-genes ~{max_genes} \
      ~{summary_references} \
      --gff-dir ref_gffs \
      --unique --sort avg-lrt-pvalue \
      > summary.tsv

    echo $(( $(wc -l < summary.tsv) - 1 )) | tee N_GENES
  >>>

  output {
    String date = read_string("DATE")
    Int n_genes = read_int("N_GENES")
    File summary = "summary.tsv"
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

task map_summary_panfeed {

  meta {
    version: "0.1.0"
    description: "Per-gene summary of significant panfeed k-mer hits. Translated from the microGWAS `run_map_summary_panfeed` rule."
  }

  input {
    String phenotype

    File panfeed_kmers            # panfeed_kmers.tsv.gz
    File phenotypes
    File gene_presence_absence_rtab
    File gene_presence_absence_csv
    Array[File] reference_gffs

    String summary_references

    Int min_hits  = 1
    Int max_genes = 25

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    mkdir -p ref_gffs
    for f in ~{sep=' ' reference_gffs}; do ln -s "$f" "ref_gffs/$(basename "$f")"; done

    # Two differences from the unitig variant, both upstream's:
    #   - /dev/null for the filtered table: the k-mer file is already filtered
    #   - no --length: k-mer length is fixed by panfeed, so the filter is moot
    #     (mapped_summary.py ignores --length under --panfeed anyway)
    python3 "$SCRIPTS/mapped_summary.py" ~{panfeed_kmers} \
      ~{phenotypes} ~{phenotype} /dev/null \
      --pangenome ~{gene_presence_absence_rtab} \
      --pangenome-genes ~{gene_presence_absence_csv} \
      --minimum-hits ~{min_hits} --maximum-genes ~{max_genes} \
      ~{summary_references} \
      --gff-dir ref_gffs \
      --unique --sort avg-lrt-pvalue \
      --panfeed \
      > summary_panfeed.tsv

    echo $(( $(wc -l < summary_panfeed.tsv) - 1 )) | tee N_GENES
  >>>

  output {
    String date = read_string("DATE")
    Int n_genes = read_int("N_GENES")
    File summary = "summary_panfeed.tsv"
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

task map_summary_wg {

  meta {
    version: "0.1.0"
    description: "Per-gene summary of the unitigs a whole-genome model selected, ranked by effect size rather than p-value. Translated from the microGWAS `run_map_summary_wg` rule."
  }

  input {
    String phenotype
    String model                  # "ridge" or "lasso"

    File mapped
    File phenotypes
    File model_variants           # <model>.tsv
    File gene_presence_absence_rtab
    File gene_presence_absence_csv
    Array[File] reference_gffs

    String summary_references

    Int length    = 10
    Int min_hits  = 1
    Int max_genes = 25

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    mkdir -p ref_gffs
    for f in ~{sep=' ' reference_gffs}; do ln -s "$f" "ref_gffs/$(basename "$f")"; done

    # --sort avg-beta, not avg-lrt-pvalue: elastic net output has no p-values,
    # so hits are ranked by mean absolute effect size
    python3 "$SCRIPTS/mapped_summary.py" ~{mapped} \
      ~{phenotypes} ~{phenotype} ~{model_variants} \
      --pangenome ~{gene_presence_absence_rtab} \
      --pangenome-genes ~{gene_presence_absence_csv} \
      --length ~{length} --minimum-hits ~{min_hits} --maximum-genes ~{max_genes} \
      ~{summary_references} \
      --gff-dir ref_gffs \
      --unique --sort avg-beta \
      > summary_~{model}.tsv
  >>>

  output {
    String date = read_string("DATE")
    File summary = "summary_~{model}.tsv"
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

task gpa_summary {

  meta {
    version: "0.1.0"
    description: "Summarise significant gene presence/absence hits against the pangenome and reference annotations. Translated from the microGWAS `run_gpa_summary` rule."
  }

  input {
    String phenotype

    File gpa_filtered
    File gene_presence_absence_rtab
    File gene_presence_absence_csv
    Array[File] reference_gffs

    String summary_references

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    mkdir -p ref_gffs
    for f in ~{sep=' ' reference_gffs}; do ln -s "$f" "ref_gffs/$(basename "$f")"; done

    python3 "$SCRIPTS/gpa_summary.py" ~{gpa_filtered} \
      --pangenome ~{gene_presence_absence_rtab} \
      --pangenome-genes ~{gene_presence_absence_csv} \
      ~{summary_references} \
      --gff-dir ref_gffs \
      --sort lrt-pvalue \
      > gpa_summary.tsv

    echo $(( $(wc -l < gpa_summary.tsv) - 1 )) | tee N_GENES
  >>>

  output {
    String date = read_string("DATE")
    Int n_genes = read_int("N_GENES")
    File summary = "gpa_summary.tsv"
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

task rare_summary {

  meta {
    version: "0.1.0"
    description: "Summarise significant rare-variant burden hits against the pangenome and reference annotations. Translated from the microGWAS `run_rare_summary` rule."
  }

  input {
    String phenotype

    File rare_filtered
    File gene_presence_absence_rtab
    File gene_presence_absence_csv
    Array[File] reference_gffs

    String summary_references
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

    mkdir -p ref_gffs
    for f in ~{sep=' ' reference_gffs}; do ln -s "$f" "ref_gffs/$(basename "$f")"; done

    # Burden-test hits are keyed by reference locus tag, so this summariser
    # takes the enrichment reference as a positional argument
    python3 "$SCRIPTS/rare_summary.py" ~{rare_filtered} ~{enrichment_reference} \
      --pangenome ~{gene_presence_absence_rtab} \
      --pangenome-genes ~{gene_presence_absence_csv} \
      ~{summary_references} \
      --gff-dir ref_gffs \
      --sort lrt-pvalue \
      > rare_summary.tsv
  >>>

  output {
    String date = read_string("DATE")
    File summary = "rare_summary.tsv"
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

task annotate_summary {

  meta {
    version: "0.1.0"
    description: "Join eggNOG COG/GO/KEGG annotations onto a hit summary table. Covers all five of the microGWAS `run_annotate_*` rules, which share an identical command body."
  }

  input {
    String phenotype

    # Label distinguishing this summary from the others in the same run, used
    # only in output filenames: "", "rare", "gpa", "panfeed", "ridge", "lasso"
    String label

    File summary
    File annotations

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # Narrow the full pangenome annotation set to just the gene clusters that
    # appear in this summary, then splice those columns onto it
    python3 "$SCRIPTS/subset_annotations.py" ~{annotations} \
      --groups ~{summary} \
      > subset.emapper.annotations

    python3 "$SCRIPTS/enhance_summary.py" ~{summary} subset.emapper.annotations \
      > annotated_summary_~{label}.tsv

    echo $(( $(wc -l < annotated_summary_~{label}.tsv) - 1 )) | tee N_ROWS
  >>>

  output {
    String date = read_string("DATE")
    Int n_rows  = read_int("N_ROWS")
    File annotated_summary = "annotated_summary_~{label}.tsv"
    File subset_annotations = "subset.emapper.annotations"
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

task annotate_panfeed_small {

  meta {
    version: "0.1.0"
    description: "Expand significant panfeed k-mer hashes back into their k-mer sequences and gene-cluster coordinates. Translated from the microGWAS `run_annotate_panfeed_small` rule."
  }

  input {
    String phenotype

    File panfeed_associations       # panfeed.tsv
    File panfeed_conversion         # kmers_to_hashes.tsv
    File panfeed_patterns           # hashes_to_patterns.tsv
    File unitigs_patterns

    Int memory    = 32
    Int cpu       = 2
    Int disk_size = 150
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    THRESHOLD=$(python3 "$SCRIPTS/count_patterns.py" --threshold ~{unitigs_patterns})
    echo "$THRESHOLD" | tee THRESHOLD

    python3 "$SCRIPTS/combine_panfeed.py" \
      ~{panfeed_associations} ~{panfeed_conversion} ~{panfeed_patterns} \
      --threshold "$THRESHOLD" \
      | gzip > panfeed_kmers.tsv.gz

    gzip -t panfeed_kmers.tsv.gz
  >>>

  output {
    String date = read_string("DATE")
    File panfeed_kmers = "panfeed_kmers.tsv.gz"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB" # TES
    preemptible: 1
    maxRetries:  2
  }
}

task panfeed_downstream {

  meta {
    version: "0.1.0"
    description: "Second-pass panfeed run restricted to significant gene clusters, producing annotated k-mers and per-cluster k-mer plots. Translated from the microGWAS `panfeed_downstream` rule."
  }

  input {
    String phenotype

    File phenotypes
    File panfeed_associations       # panfeed.tsv
    File panfeed_conversion         # kmers_to_hashes.tsv
    File unitigs_patterns
    File gene_presence_absence_csv

    Array[File] sample_gffs
    Array[File] reference_gffs = []
    Array[File] sample_fastas
    Array[File] reference_fastas = []

    Int upstream_bp   = 250
    Int downstream_bp = 100

    Int memory    = 32
    Int cpu       = 4
    Int disk_size = 200
    String docker_image = "aarvani1/microgwas-panfeed:1.6.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    gffs=(~{sep=' ' sample_gffs} ~{sep=' ' reference_gffs})
    printf '%s\n' "${gffs[@]}" > gff_list.txt
    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    printf '%s\n' "${fastas[@]}" > fasta_list.txt

    THRESHOLD=$(python3 "$SCRIPTS/count_patterns.py" --threshold ~{unitigs_patterns})
    echo "$THRESHOLD" | tee THRESHOLD

    # Which gene clusters cleared the threshold - the second pass only needs
    # to regenerate k-mers for these, not the whole pangenome
    panfeed-get-clusters \
      -a ~{panfeed_associations} \
      -p ~{panfeed_conversion} \
      -t "$THRESHOLD" \
      > panfeed_clusters.txt

    # Restrict to the study samples (references are excluded from the plots)
    for f in ~{sep=' ' sample_fastas}; do
      b=$(basename "$f"); echo "${b%.*}"
    done | sort | uniq > panfeed_targets.txt

    rm -rf panfeed_second_pass
    panfeed \
      --gff gff_list.txt \
      --fasta fasta_list.txt \
      -o panfeed_second_pass \
      -p ~{gene_presence_absence_csv} \
      --upstream ~{upstream_bp} --downstream ~{downstream_bp} \
      --no-filter \
      -v \
      --genes panfeed_clusters.txt \
      --targets panfeed_targets.txt \
      --cores ~{cpu}

    panfeed-get-kmers \
      -a ~{panfeed_associations} \
      -p panfeed_second_pass/kmers_to_hashes.tsv \
      -k panfeed_second_pass/kmers.tsv \
      | gzip > panfeed_annotated_kmers.tsv.gz

    gzip -t panfeed_annotated_kmers.tsv.gz

    mkdir -p panfeed_plots
    panfeed-plot \
      -k "$(realpath panfeed_annotated_kmers.tsv.gz)" \
      -p "$(realpath ~{phenotypes})" \
      --phenotype-column ~{phenotype} \
      -t "$THRESHOLD" \
      --output-directory panfeed_plots

    rm -rf panfeed_second_pass

    ls -1 panfeed_plots | wc -l | tee N_PLOTS
  >>>

  output {
    String date = read_string("DATE")
    Int n_plots = read_int("N_PLOTS")
    File panfeed_annotated_kmers = "panfeed_annotated_kmers.tsv.gz"
    File panfeed_clusters        = "panfeed_clusters.txt"
    Array[File] panfeed_plots    = glob("panfeed_plots/*")
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
