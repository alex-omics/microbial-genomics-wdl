version 1.0

# task_pyseer.wdl
#
# The association-testing core: everything backed by microGWAS's
# workflow/envs/pyseer.yaml.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: prepare_pyseer, prepare_wg,
# run_pyseer, run_pyseer_vcf, run_panfeed, run_wg, run_wg_metrics, run_qq_plot.
#
# All of these are scattered per phenotype (`{target}` upstream). The
# `covariates` input is the raw pyseer flag string upstream substitutes
# verbatim from config["covariates"][target], e.g. "--use-covariates 6q 7"
# (1-based column indices into the phenotype file, "q" suffix = quantitative).
# The Tier 2 E. coli validation set has no covariate columns, so it stays empty.
#
# Shell-safety note: upstream's `shell:` blocks run without `pipefail`, and a
# few of them depend on that leniency - `grep -v 'cv.glmnet'` legitimately
# filters every line when glmnet is quiet, and `head -1` on a large file
# SIGPIPEs its producer. Those pipelines are decomposed into separate steps
# here so `set -euxo pipefail` can stay on without changing behaviour.

task prepare_pyseer {

  meta {
    version: "0.1.0"
    description: "Subset and align phenotype, similarity, distance and lineage inputs to the samples with data for one phenotype. Translated from the microGWAS `prepare_pyseer` rule."
  }

  input {
    String phenotype

    # Sample and reference assemblies, used only to recover the genome ID list
    Array[File] sample_fastas
    Array[File] reference_fastas = []

    File phenotypes_tsv
    File all_similarities
    File all_distances
    File all_lineages

    Int memory    = 8
    Int cpu       = 1
    Int disk_size = 20
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # prepare_pyseer.py reads this file only for its index (the genome list)
    printf 'ID\tPath\n' > unitigs_input.tsv
    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    for f in "${fastas[@]}"; do
      id=$(basename "$f"); id="${id%.*}"
      printf '%s\t%s\n' "$id" "$f" >> unitigs_input.tsv
    done

    mkdir -p prepared

    # Drops samples with a NaN phenotype, then intersects the sample sets of
    # the genome list, phenotypes, similarity, distances and lineages, writing
    # a matched bundle. Reference genomes fall out here naturally - they carry
    # no phenotype.
    python3 "$SCRIPTS/prepare_pyseer.py" \
      unitigs_input.tsv \
      ~{phenotypes_tsv} \
      ~{all_similarities} \
      ~{all_distances} \
      ~{all_lineages} \
      prepared \
      ~{phenotype} | tee PREPARE_LOG

    # Fail loudly rather than letting pyseer run on an empty intersection,
    # which it reports only as an obscure downstream error
    n=$(( $(wc -l < prepared/phenotypes.tsv) - 1 ))
    echo "$n" | tee N_SHARED
    if [[ "$n" -lt 1 ]]; then
      echo "ERROR: no samples shared across phenotype/similarity/distance/lineage inputs for '~{phenotype}'" >&2
      exit 1
    fi
  >>>

  output {
    String date       = read_string("DATE")
    String prepare_log = read_string("PREPARE_LOG")
    Int n_shared_samples = read_int("N_SHARED")

    File phenotypes = "prepared/phenotypes.tsv"
    File similarity = "prepared/similarity.tsv"
    File distances  = "prepared/distances.tsv"
    File lineages   = "prepared/lineages.tsv"
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

task run_pyseer {

  meta {
    version: "0.1.0"
    description: "Run the pyseer LMM against unitigs, gene presence/absence and structural variants for one phenotype, and apply the pattern-counted significance threshold. Translated from the microGWAS `run_pyseer` rule."
  }

  input {
    String phenotype

    File unitigs
    File gene_presence_absence_rtab
    File struct_presence_absence

    File phenotypes
    File similarity
    File distances
    File lineages

    String covariates = ""

    Int memory    = 32
    Int cpu       = 2
    Int disk_size = 200
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    pyseer --version 2>&1 | tee VERSION

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # --- Lineage-effect pass -------------------------------------------------
    # Upstream runs this first against a 10-line slice of the unitig file. The
    # slice is deliberate: the run exists to produce unitigs_lineage.txt (the
    # per-lineage effect table), not association statistics, so pyseer is given
    # only enough input to initialise.
    zcat ~{unitigs} | head -n 10 > small.txt || true

    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --kmers small.txt \
           --uncompressed \
           --cpu ~{cpu} \
           --lineage --lineage-clusters ~{lineages} \
           --lineage-file unitigs_lineage.txt \
           --distances ~{distances} \
           ~{covariates} \
           > lineage_pass.log 2>&1

    # --- Unitigs -------------------------------------------------------------
    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --kmers ~{unitigs} \
           --similarity ~{similarity} \
           --lmm \
           --output-patterns unitigs_patterns.txt \
           --cpu ~{cpu} \
           ~{covariates} \
           > unitigs.tsv

    # The significance threshold is 0.05 / (number of unique presence/absence
    # patterns), which is far less conservative than counting every variant.
    THRESHOLD=$(python3 "$SCRIPTS/count_patterns.py" --threshold unitigs_patterns.txt)
    echo "$THRESHOLD" | tee UNITIGS_THRESHOLD

    head -1 unitigs.tsv > unitigs_filtered.tsv
    LC_ALL=C awk -v pval="$THRESHOLD" '$4<pval {print $0}' unitigs.tsv >> unitigs_filtered.tsv

    # --- Gene presence/absence ----------------------------------------------
    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --pres ~{gene_presence_absence_rtab} \
           --similarity ~{similarity} \
           --lmm --uncompressed \
           --output-patterns gpa_patterns.txt \
           --cpu ~{cpu} \
           ~{covariates} \
           > gpa.tsv

    # NOTE: upstream filters gpa and struct against the *unitigs* threshold,
    # not their own pattern counts, even though it writes gpa_patterns.txt and
    # struct_patterns.txt. Preserved verbatim - changing it would shift which
    # hits pass and break comparability with the published results.
    head -1 gpa.tsv > gpa_filtered.tsv
    LC_ALL=C awk -v pval="$THRESHOLD" '$4<pval {print $0}' gpa.tsv >> gpa_filtered.tsv

    # --- Structural variants -------------------------------------------------
    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --pres ~{struct_presence_absence} \
           --similarity ~{similarity} \
           --lmm --uncompressed \
           --output-patterns struct_patterns.txt \
           --cpu ~{cpu} \
           ~{covariates} \
           > struct.tsv

    head -1 struct.tsv > struct_filtered.tsv
    LC_ALL=C awk -v pval="$THRESHOLD" '$4<pval {print $0}' struct.tsv >> struct_filtered.tsv

    # Hit counts, minus the header line, for at-a-glance QC in Terra
    echo $(( $(wc -l < unitigs_filtered.tsv) - 1 )) | tee N_UNITIG_HITS
    echo $(( $(wc -l < gpa_filtered.tsv) - 1 ))     | tee N_GPA_HITS
    echo $(( $(wc -l < struct_filtered.tsv) - 1 ))  | tee N_STRUCT_HITS
  >>>

  output {
    String date            = read_string("DATE")
    String pyseer_version  = read_string("VERSION")
    String pyseer_docker_image = docker_image
    String significance_threshold = read_string("UNITIGS_THRESHOLD")

    Int n_unitig_hits = read_int("N_UNITIG_HITS")
    Int n_gpa_hits    = read_int("N_GPA_HITS")
    Int n_struct_hits = read_int("N_STRUCT_HITS")

    File unitigs_patterns = "unitigs_patterns.txt"
    File unitigs_results  = "unitigs.tsv"
    File unitigs_filtered = "unitigs_filtered.tsv"
    File unitigs_lineage  = "unitigs_lineage.txt"
    File gpa_results      = "gpa.tsv"
    File gpa_filtered     = "gpa_filtered.tsv"
    File struct_results   = "struct.tsv"
    File struct_filtered  = "struct_filtered.tsv"
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

task run_pyseer_vcf {

  meta {
    version: "0.1.0"
    description: "Run the pyseer LMM against common SNPs, and optionally the rare-variant gene burden test. Translated from the microGWAS `run_pyseer_vcf` rule."
  }

  input {
    String phenotype

    File common_snps
    File common_snps_index

    # Rare variant set plus the CDS region file that defines burden groupings.
    # Both optional - see task_rare_variants.wdl for why the rare path is
    # off by default in Milestone 1.
    File? rare_snps
    File? rare_snps_index
    File? regions

    File phenotypes
    File similarity
    # As in run_panfeed, upstream's `lineages` input here is DAG ordering only
    # and is never read by the command; dropped.

    # Threshold source: upstream reuses the unitig pattern count here too
    File unitigs_patterns

    String covariates = ""

    Int memory    = 32
    Int cpu       = 1
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    THRESHOLD=$(python3 "$SCRIPTS/count_patterns.py" --threshold ~{unitigs_patterns})
    echo "$THRESHOLD" | tee THRESHOLD

    # Co-locate VCFs with their indexes for pysam
    mkdir -p vcf_in
    ln -s ~{common_snps} vcf_in/common.vcf.gz
    ln -s ~{common_snps_index} vcf_in/common.vcf.gz.csi

    # --- Rare variants (gene burden test) ------------------------------------
    ~{if defined(rare_snps) then "RUN_RARE=1" else "RUN_RARE=0"}
    if [[ "$RUN_RARE" == "1" ]]; then
      ln -s ~{default="/dev/null" rare_snps} vcf_in/rare.vcf.gz
      ln -s ~{default="/dev/null" rare_snps_index} vcf_in/rare.vcf.gz.csi

      pyseer --phenotypes ~{phenotypes} \
             --phenotype-column ~{phenotype} \
             --vcf vcf_in/rare.vcf.gz \
             --burden ~{default="" regions} \
             --similarity ~{similarity} \
             --lmm \
             --output-patterns rare_patterns.txt \
             --cpu ~{cpu} \
             ~{covariates} \
             > rare.tsv

      head -1 rare.tsv > rare_filtered.tsv
      LC_ALL=C awk -v pval="$THRESHOLD" '$4<pval {print $0}' rare.tsv >> rare_filtered.tsv
    else
      echo "No rare variant VCF supplied - skipping the gene burden test." | tee RARE_SKIPPED
    fi

    # --- Common variants -----------------------------------------------------
    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --vcf vcf_in/common.vcf.gz \
           --similarity ~{similarity} \
           --lmm \
           --output-patterns vcf_patterns.txt \
           --cpu ~{cpu} \
           ~{covariates} \
           > vcf.tsv

    head -1 vcf.tsv > vcf_filtered.tsv
    LC_ALL=C awk -v pval="$THRESHOLD" '$4<pval {print $0}' vcf.tsv >> vcf_filtered.tsv

    # Attach the VCF's own annotation fields to the significant common hits
    python3 "$SCRIPTS/vcf_subset.py" vcf_filtered.tsv vcf_in/common.vcf.gz > annotated_vcf.tsv

    echo $(( $(wc -l < vcf_filtered.tsv) - 1 )) | tee N_COMMON_HITS
  >>>

  output {
    String date = read_string("DATE")
    Int n_common_hits = read_int("N_COMMON_HITS")

    File common_results   = "vcf.tsv"
    File common_filtered  = "vcf_filtered.tsv"
    File common_annotated = "annotated_vcf.tsv"
    File common_patterns  = "vcf_patterns.txt"

    File? rare_results  = "rare.tsv"
    File? rare_filtered = "rare_filtered.tsv"
    File? rare_patterns = "rare_patterns.txt"
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

task run_panfeed {

  meta {
    version: "0.1.0"
    description: "Run the pyseer LMM against panfeed's gene-cluster-specific k-mer patterns. Translated from the microGWAS `run_panfeed` rule."
  }

  input {
    String phenotype

    File panfeed_patterns

    File phenotypes
    File similarity
    # NOTE: upstream also declares distances and lineages as inputs to this
    # rule, but its shell body never reads them - they are there only to order
    # the Snakemake DAG. Cromwell orders on the explicit call graph instead, so
    # they are dropped rather than localized for nothing.
    File unitigs_patterns

    String covariates = ""

    Int memory    = 32
    Int cpu       = 2
    Int disk_size = 150
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    THRESHOLD=$(python3 "$SCRIPTS/count_patterns.py" --threshold ~{unitigs_patterns})
    echo "$THRESHOLD" | tee THRESHOLD

    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --pres ~{panfeed_patterns} \
           --similarity ~{similarity} \
           --lmm --uncompressed \
           --output-patterns panfeed_patterns_out.txt \
           --cpu ~{cpu} \
           ~{covariates} \
           > panfeed.tsv

    head -1 panfeed.tsv > panfeed_filtered.tsv
    LC_ALL=C awk -v pval="$THRESHOLD" '$4<pval {print $0}' panfeed.tsv >> panfeed_filtered.tsv

    echo $(( $(wc -l < panfeed_filtered.tsv) - 1 )) | tee N_PANFEED_HITS
  >>>

  output {
    String date = read_string("DATE")
    Int n_panfeed_hits = read_int("N_PANFEED_HITS")

    File panfeed_results  = "panfeed.tsv"
    File panfeed_filtered = "panfeed_filtered.tsv"
    File panfeed_patterns_out = "panfeed_patterns_out.txt"
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

task prepare_wg {

  meta {
    version: "0.1.0"
    description: "Prepare per-phenotype inputs for the whole-genome elastic net and cache the parsed unitig matrix as a pickle. Translated from the microGWAS `prepare_wg` rule."
  }

  input {
    String phenotype

    Array[File] sample_fastas
    Array[File] reference_fastas = []

    File unitigs
    File phenotypes_tsv
    File all_similarities
    File all_distances
    File all_lineages

    String covariates = ""

    Int memory    = 64
    Int cpu       = 4
    Int disk_size = 300
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    printf 'ID\tPath\n' > unitigs_input.tsv
    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    for f in "${fastas[@]}"; do
      id=$(basename "$f"); id="${id%.*}"
      printf '%s\t%s\n' "$id" "$f" >> unitigs_input.tsv
    done

    mkdir -p prepared

    python3 "$SCRIPTS/prepare_pyseer.py" \
      unitigs_input.tsv \
      ~{phenotypes_tsv} \
      ~{all_similarities} \
      ~{all_distances} \
      ~{all_lineages} \
      prepared \
      ~{phenotype}

    # This pyseer run exists solely to parse and cache the unitig matrix.
    # --save-vars writes <prefix>.pkl, which run_wg then --load-vars back in,
    # so ridge and lasso each avoid re-parsing the full unitig file.
    # --cor-filter 0.0 disables correlation pre-filtering at this stage.
    pyseer --phenotypes prepared/phenotypes.tsv \
           --phenotype-column ~{phenotype} \
           --kmers ~{unitigs} \
           --wg enet --save-vars prepared/variants \
           --cor-filter 0.0 \
           ~{covariates} \
           --cpu ~{cpu}

    ls -l prepared/variants.pkl
  >>>

  output {
    String date = read_string("DATE")

    File phenotypes = "prepared/phenotypes.tsv"
    File similarity = "prepared/similarity.tsv"
    File distances  = "prepared/distances.tsv"
    File lineages   = "prepared/lineages.tsv"
    File variants_pkl = "prepared/variants.pkl"
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

task run_wg {

  meta {
    version: "0.1.0"
    description: "Fit the whole-genome elastic net models (ridge alpha=0.01, lasso alpha=1) over the full unitig set. Translated from the microGWAS `run_wg` rule."
  }

  input {
    String phenotype

    File unitigs
    File variants_pkl
    File phenotypes
    # Upstream declares `similarity` here too, but run_wg's pyseer invocations
    # pass --distances and --lineage-clusters, never --similarity; dropped.
    File distances
    File lineages

    # Upstream's alpha values: 0.01 is ridge-like, 1 is pure lasso. The paper's
    # headline whole-genome R^2 (~0.48 for the E. coli validation) is the lasso.
    Float ridge_alpha = 0.01
    Float lasso_alpha = 1

    String covariates = ""

    Int memory    = 64
    Int cpu       = 4
    Int disk_size = 300
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    # --load-vars takes the prefix, not the .pkl path
    mkdir -p vars
    cp ~{variants_pkl} vars/variants.pkl

    # --- Ridge ---------------------------------------------------------------
    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --kmers ~{unitigs} \
           --distances ~{distances} \
           --wg enet \
           --load-vars vars/variants \
           --save-model ridge \
           --alpha ~{ridge_alpha} \
           --sequence-reweighting \
           --lineage-clusters ~{lineages} \
           --cor-filter 0.0 \
           --save-predictions ridge_predictions.tsv \
           ~{covariates} \
           > ridge_raw.tsv 2> ridge.log

    # glmnet chatters onto stdout; upstream strips it. Decoupled from the pipe
    # so a fully-filtered file does not trip pipefail.
    grep -v 'cv.glmnet' ridge_raw.tsv > ridge.tsv || true

    # --- Lasso ---------------------------------------------------------------
    pyseer --phenotypes ~{phenotypes} \
           --phenotype-column ~{phenotype} \
           --kmers ~{unitigs} \
           --distances ~{distances} \
           --wg enet \
           --load-vars vars/variants \
           --save-model lasso \
           --alpha ~{lasso_alpha} \
           --sequence-reweighting \
           --lineage-clusters ~{lineages} \
           --cor-filter 0.0 \
           --save-predictions lasso_predictions.tsv \
           ~{covariates} \
           > lasso_raw.tsv 2> lasso.log

    grep -v 'cv.glmnet' lasso_raw.tsv > lasso.tsv || true

    rm -f ridge_raw.tsv lasso_raw.tsv
  >>>

  output {
    String date = read_string("DATE")

    File ridge             = "ridge.tsv"
    File ridge_predictions = "ridge_predictions.tsv"
    File ridge_model       = "ridge.pkl"
    File ridge_log         = "ridge.log"
    File lasso             = "lasso.tsv"
    File lasso_predictions = "lasso_predictions.tsv"
    File lasso_model       = "lasso.pkl"
    File lasso_log         = "lasso.log"
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

task run_wg_metrics {

  meta {
    version: "0.1.0"
    description: "Compute predictive performance metrics (including R^2) for one whole-genome model. Translated from the microGWAS `run_wg_metrics` rule."
  }

  input {
    String phenotype
    String model          # "ridge" or "lasso"
    File predictions

    Int memory    = 8
    Int cpu       = 1
    Int disk_size = 20
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    python3 "$SCRIPTS/prediction_metrics.py" ~{predictions} > metrics_~{model}.tsv

    cat metrics_~{model}.tsv
  >>>

  output {
    String date = read_string("DATE")
    File metrics = "metrics_~{model}.tsv"
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

task run_qq_plot {

  meta {
    version: "0.1.0"
    description: "Draw a QQ plot of observed vs expected p-values for one variant class, to check for residual population-structure inflation. Translated from the microGWAS `run_qq_plot` rule."
  }

  input {
    String phenotype
    String variant_class   # unitigs | gpa | rare | vcf | panfeed
    File association_results

    Int memory    = 8
    Int cpu       = 1
    Int disk_size = 20
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    python3 "$SCRIPTS/qq_plot.py" ~{association_results} \
      --output qq_~{variant_class}.png
  >>>

  output {
    String date = read_string("DATE")
    File qq_plot = "qq_~{variant_class}.png"
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
