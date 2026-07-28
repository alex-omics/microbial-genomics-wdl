version 1.0

# task_population_structure.wdl
#
# Population-structure and relatedness tasks: the kinship/distance/lineage
# inputs that let pyseer's LMM separate real associations from lineage effects.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: mash_sketch, distance, lineage_st,
# tree, variant_sites, aln2vcf, similarity.
#
# Sample-name safety note: every tool here derives sample IDs from file
# basenames (square_mash does os.path.split(...)[-1].split('.')[0];
# sanitize_mlst.py strips the path and final extension), so Cromwell's
# arbitrary localization paths do not leak into the matrices that
# prepare_pyseer later intersects on sample ID.

task mash_sketch {

  meta {
    version: "0.1.0"
    description: "Sketch all genomes with mash for pairwise distance estimation. Translated from the microGWAS `mash_sketch` rule."
  }

  input {
    Array[File] sample_fastas
    Array[File] reference_fastas = []

    Int sketch_size = 10000

    Int memory    = 16
    Int cpu       = 5
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-mash:2.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    mash --version | tee VERSION

    # Rebuild out/mash_input.txt: one genome path per line
    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    printf '%s\n' "${fastas[@]}" > mash_input.txt
    echo "Genomes to sketch: $(wc -l < mash_input.txt)"

    mash sketch -p ~{cpu} -s ~{sketch_size} -o sketches -l mash_input.txt
  >>>

  output {
    String date         = read_string("DATE")
    String mash_version = read_string("VERSION")
    String mash_docker_image = docker_image
    File sketches       = "sketches.msh"
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

task distance {

  meta {
    version: "0.1.0"
    description: "Compute the square mash distance matrix used as pyseer's population-structure covariate. Translated from the microGWAS `distance` rule."
  }

  input {
    File sketches

    Int memory    = 16
    Int cpu       = 5
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-mash:2.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    # square_mash is pyseer's helper (pyseer.mash:main); it reshapes mash's
    # long-format all-vs-all output into a square matrix and reduces each
    # genome path to its basename stem for the row/column labels.
    mash dist -p ~{cpu} ~{sketches} ~{sketches} | square_mash > distances.tsv

    head -1 distances.tsv | awk '{print NF}' | tee N_GENOMES
  >>>

  output {
    String date    = read_string("DATE")
    Int n_genomes  = read_int("N_GENOMES")
    File distances = "distances.tsv"
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

task lineage_st {

  meta {
    version: "0.1.0"
    description: "Assign MLST sequence types to every genome, used as the lineage/clustering covariate. Translated from the microGWAS `lineage_st` rule."
  }

  input {
    Array[File] sample_fastas
    Array[File] reference_fastas = []

    # config["mlst_scheme"], e.g. "ecoli"
    String mlst_scheme

    Int memory    = 16
    Int cpu       = 16
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-mlst:2.16.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    mlst --version | tee VERSION

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    fastas=(~{sep=' ' sample_fastas} ~{sep=' ' reference_fastas})
    printf '%s\n' "${fastas[@]}" > mash_input.txt

    mkdir -p mlst_tmp

    # Verbatim from the rule: build a job file of one `mlst` call per genome,
    # fan out with GNU parallel, then concatenate. `cut -f1,3` keeps the
    # filename and the ST call, dropping the scheme column.
    for i in $(cat mash_input.txt); do
      echo "mlst --scheme ~{mlst_scheme} --threads 1 \"$i\" | cut -f1,3 > mlst_tmp/$(basename "$i").mlst"
    done > mlst_tmp/mlst_jobs.txt

    parallel -j ~{cpu} --progress < mlst_tmp/mlst_jobs.txt

    cat mlst_tmp/*.mlst > mlst_tmp/all_mlst.txt

    # sanitize_mlst.py strips paths and the file extension so the lineage table
    # is keyed by bare sample ID
    python3 "$SCRIPTS/sanitize_mlst.py" mlst_tmp/all_mlst.txt > lineages_mlst.txt

    rm -rf mlst_tmp

    wc -l < lineages_mlst.txt | tee N_TYPED
    cut -f2 lineages_mlst.txt | sort -u | wc -l | tee N_LINEAGES
  >>>

  output {
    String date         = read_string("DATE")
    String mlst_version = read_string("VERSION")
    String mlst_docker_image = docker_image
    Int n_typed         = read_int("N_TYPED")
    Int n_lineages      = read_int("N_LINEAGES")
    File lineages       = "lineages_mlst.txt"
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

task tree {

  meta {
    version: "0.1.0"
    description: "Build an approximately-maximum-likelihood core genome phylogeny with FastTree. Translated from the microGWAS `tree` rule."
  }

  input {
    File core_genome_aln

    Int memory    = 32
    Int cpu       = 4
    Int disk_size = 100
    String docker_image = "staphb/fasttree:2.1.11"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    FastTree -nt -gtr < ~{core_genome_aln} > core_gene_alignment.aln.treefile
  >>>

  output {
    String date = read_string("DATE")
    String fasttree_docker_image = docker_image
    File treefile = "core_gene_alignment.aln.treefile"
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

task variant_sites {

  meta {
    version: "0.1.0"
    description: "Reduce the core genome alignment to variable sites only. Translated from the microGWAS `variant_sites` rule."
  }

  input {
    File core_genome_aln

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-snpsites:2.5.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    snp-sites -V 2>&1 | tee VERSION

    snp-sites -o core_gene_alignment_variable.aln ~{core_genome_aln}
  >>>

  output {
    String date             = read_string("DATE")
    String snp_sites_version = read_string("VERSION")
    File variant_alignment  = "core_gene_alignment_variable.aln"
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

task aln2vcf {

  meta {
    version: "0.1.0"
    description: "Convert the core genome alignment to a normalised, indexed VCF of core SNPs. Translated from the microGWAS `aln2vcf` rule."
  }

  input {
    File core_genome_aln

    Int memory    = 16
    Int cpu       = 16
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-snpsites:2.5.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    snp-sites -o core.vcf -v ~{core_genome_aln}

    # `norm -m -` splits multiallelic records into biallelic ones, which is
    # what similarity.py's per-variant presence/absence logic assumes
    bcftools norm -m - -O z --threads ~{cpu} core.vcf > core_gene_alignment.vcf.gz
    bcftools index core_gene_alignment.vcf.gz

    rm core.vcf
  >>>

  output {
    String date          = read_string("DATE")
    File core_genome_vcf = "core_gene_alignment.vcf.gz"
    File core_genome_vcf_index = "core_gene_alignment.vcf.gz.csi"
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

task similarity {

  meta {
    version: "0.1.0"
    description: "Compute the core-SNP-derived kinship matrix that pyseer's LMM uses to correct for population structure. Translated from the microGWAS `similarity` rule."
  }

  input {
    File core_genome_aln
    File core_genome_vcf
    File core_genome_vcf_index

    Int memory    = 32
    Int cpu       = 2
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-pyseer:1.4.2"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    pyseer --version 2>&1 | tee PYSEER_VERSION

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # The sample list comes from the alignment headers, so the kinship matrix
    # is keyed by the same strain names Panaroo used
    grep '>' ~{core_genome_aln} | cut -c 2- > similarity_samples.txt

    # Co-locate the VCF and its index; pysam requires them side by side and
    # Cromwell may localize them into separate directories
    mkdir -p vcf_in
    ln -s ~{core_genome_vcf} vcf_in/core.vcf.gz
    ln -s ~{core_genome_vcf_index} vcf_in/core.vcf.gz.csi

    python3 "$SCRIPTS/similarity.py" similarity_samples.txt \
      --vcf vcf_in/core.vcf.gz \
      > similarity.tsv

    rm similarity_samples.txt

    head -1 similarity.tsv | awk '{print NF}' | tee N_SAMPLES
  >>>

  output {
    String date           = read_string("DATE")
    String pyseer_version = read_string("PYSEER_VERSION")
    String pyseer_docker_image = docker_image
    Int n_samples         = read_int("N_SAMPLES")
    File similarities     = "similarity.tsv"
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
