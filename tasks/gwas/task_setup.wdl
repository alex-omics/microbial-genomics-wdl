version 1.0

# task_setup.wdl
#
# Input preparation tasks for the multimodal_gwas workflow.
#
# Translated from microGWAS (https://github.com/microbial-pangenomes-lab/microGWAS,
# MIT-licensed, Copyright (c) 2022, Marco Galardini), pinned at commit 1307250.
#
# Covers `bootstrap.sh` (which upstream runs by hand before Snakemake starts)
# plus the `prepare_regions` and `download_obo` rules.
#
# A note on manifests, which is the single biggest structural difference
# between the Snakemake and WDL versions of this pipeline:
#
#   Upstream passes tools "file-of-files" manifests (out/mash_input.txt,
#   out/unitigs_input.tsv, out/panaroo_input.txt, ...) containing *paths*,
#   written once by bootstrap.sh and read by many rules. That works because
#   every Snakemake rule runs in the same working directory.
#
#   Under Cromwell each task runs in its own container with its own localized
#   paths, so a manifest written by one task is meaningless to the next. Every
#   task that needs a manifest therefore rebuilds it in its own command block
#   from the Array[File] it was given. Sample IDs come from the file basename,
#   which is safe here because microGWAS already requires SAMPLE.fasta /
#   SAMPLE.gff naming (enforced upstream by aid_bootstrap.py).

task bootstrap_references {

  meta {
    version: "0.1.0"
    description: "Download and normalise reference genomes from NCBI RefSeq, producing the per-strain FASTA/GFF/GBK/FAA sets microGWAS expects. Translated from bootstrap.sh."
  }

  input {
    # NCBI RefSeq assembly accessions, e.g. ["GCF_000013305.1", "GCF_000007445.1"]
    Array[String] assembly_accessions

    # Genus/species, used only to locate the ncbi-genome-download output tree
    String genus
    String species

    # Strain name of the primary reference used for rare-variant calling and
    # enrichment. Must match one of the downloaded assemblies' strain names
    # (e.g. "IAI39") or one of the local reference basenames.
    String reference_strain

    # Optional pre-downloaded references (the `--local-dirs` path in
    # bootstrap.sh). Provide matched triples; the FASTA basename is taken as
    # the strain name, so 536.fasta / 536.gff / 536.gbk.
    Array[File] local_reference_fastas = []
    Array[File] local_reference_gffs   = []
    Array[File] local_reference_gbks   = []

    Int memory    = 8
    Int cpu       = 2
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    ncbi-genome-download --version | tee NGD_VERSION

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    mkdir -p data/references/fastas \
             data/references/gffs \
             data/references/gbks \
             data/references/faas

    HR="data/references/human_readable/refseq/bacteria/~{genus}/~{species}"
    mkdir -p "$HR"

    # --- Local references (bootstrap.sh --local-dirs branch) ------------------
    # Upstream expects genome.fasta/genome.gff/genome.gbk inside a directory
    # whose basename is the strain name. Here the strain name comes from the
    # FASTA basename instead, since WDL hands us files rather than directories.
    local_fastas=(~{sep=' ' local_reference_fastas})
    local_gffs=(~{sep=' ' local_reference_gffs})
    local_gbks=(~{sep=' ' local_reference_gbks})

    if [[ ${#local_fastas[@]} -gt 0 ]]; then
      if [[ ${#local_fastas[@]} -ne ${#local_gffs[@]} || ${#local_fastas[@]} -ne ${#local_gbks[@]} ]]; then
        echo "ERROR: local reference FASTA/GFF/GBK arrays must be the same length" >&2
        exit 1
      fi
      for i in "${!local_fastas[@]}"; do
        ref=$(basename "${local_fastas[$i]}")
        ref="${ref%.*}"
        mkdir -p "$HR/$ref"
        gzip -c "${local_fastas[$i]}" > "$HR/$ref/${ref}_genomic.fna.gz"
        gzip -c "${local_gffs[$i]}"   > "$HR/$ref/${ref}_genomic.gff.gz"
        gzip -c "${local_gbks[$i]}"   > "$HR/$ref/${ref}_genomic.gbff.gz"
        # bootstrap.sh generates the proteome from the GenBank file when no
        # .faa was supplied alongside the local assembly
        python3 "$SCRIPTS/gbk2faa.py" "${local_gbks[$i]}" > "${ref}.faa"
        gzip -c "${ref}.faa" > "$HR/$ref/${ref}_protein.faa.gz"
      done
    fi

    # --- RefSeq references ----------------------------------------------------
    accessions="~{sep=',' assembly_accessions}"
    if [[ -n "$accessions" ]]; then
      ncbi-genome-download -H \
        -F gff,fasta,genbank,protein-fasta \
        -A "$accessions" \
        -p 1 \
        -o data/references \
        bacteria
    fi

    if [[ -z "$(ls -A "$HR" 2>/dev/null)" ]]; then
      echo "ERROR: no reference genomes were obtained (neither local nor RefSeq)" >&2
      exit 1
    fi

    # --- Normalise every reference into the flat per-type directories ---------
    # Verbatim from bootstrap.sh, including the RefSeq -> Prokka GFF conversion
    # that makes reference annotations comparable with the sample annotations.
    for ref in $(ls "$HR"); do
      echo "Normalising reference: $ref"
      zcat "$HR/$ref/"*_genomic.fna.gz > "data/references/fastas/$ref.fasta"
      tempfile=$(mktemp)
      zcat "$HR/$ref/"*_genomic.gff.gz > "$tempfile"
      python3 "$SCRIPTS/convert_refseq_to_prokka_gff.py" \
        -g "$tempfile" \
        -f "data/references/fastas/$ref.fasta" \
        -o "data/references/gffs/$ref.gff"
      zcat "$HR/$ref/"*_genomic.gbff.gz  > "data/references/gbks/$ref.gbk"
      zcat "$HR/$ref/"*_protein.faa.gz   > "data/references/faas/$ref.faa"
    done

    # --- Primary reference for rare variants / enrichment ---------------------
    if [[ ! -f "data/references/gbks/~{reference_strain}.gbk" ]]; then
      echo "ERROR: reference_strain '~{reference_strain}' not found among downloaded references:" >&2
      ls data/references/gbks/ >&2
      exit 1
    fi

    python3 "$SCRIPTS/gbk2faa.py" "data/references/gbks/~{reference_strain}.gbk" > reference.faa
    cp "data/references/gffs/~{reference_strain}.gff" reference.gff
    cp "data/references/gbks/~{reference_strain}.gbk" reference.gbk

    ls -1 data/references/fastas/*.fasta | wc -l | tee N_REFERENCES
  >>>

  output {
    String date                     = read_string("DATE")
    String ncbi_genome_download_version = read_string("NGD_VERSION")
    String bootstrap_docker_image   = docker_image
    Int n_references                = read_int("N_REFERENCES")

    Array[File] reference_fastas = glob("data/references/fastas/*.fasta")
    Array[File] reference_gffs   = glob("data/references/gffs/*.gff")
    Array[File] reference_gbks   = glob("data/references/gbks/*.gbk")
    Array[File] reference_faas   = glob("data/references/faas/*.faa")

    # The single primary reference, used by get_snps and prepare_regions
    File snps_reference_faa = "reference.faa"
    File snps_reference_gff = "reference.gff"
    File snps_reference_gbk = "reference.gbk"
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

task prepare_regions {

  meta {
    version: "0.1.0"
    description: "Extract CDS regions from the reference GFF into the burden-test region file consumed by pyseer --burden. Translated from the microGWAS `prepare_regions` rule."
  }

  input {
    File reference_gff

    Int memory    = 4
    Int cpu       = 1
    Int disk_size = 10
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    # Verbatim from the `prepare_regions` rule. The `sed 's/.1:/:/g'` strips
    # the RefSeq accession version suffix so contig names match those snippy
    # emits into the VCF.
    grep CDS ~{reference_gff} | \
    awk -F '\t' '{print $1":"$4"-"$5"\t"$9}' | \
    sed 's/.1:/:/g' | awk -F ';' '{print $1}' | \
    sed 's/ID=//g' | awk '{print $2"\t"$1}' > regions.tsv

    wc -l < regions.tsv | tee N_REGIONS
  >>>

  output {
    String date        = read_string("DATE")
    Int n_regions      = read_int("N_REGIONS")
    File regions       = "regions.tsv"
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

task download_obo {

  meta {
    version: "0.1.0"
    description: "Fetch the basic Gene Ontology DAG used by the functional enrichment tasks. Translated from the microGWAS `download_obo` rule."
  }

  input {
    String obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo"

    Int memory    = 2
    Int cpu       = 1
    Int disk_size = 10
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    wget -O go-basic.obo '~{obo_url}'

    # A truncated or HTML error page would otherwise surface much later as a
    # confusing goatools parse failure inside the enrichment tasks
    grep -q '^format-version:' go-basic.obo
    grep -c '^\[Term\]' go-basic.obo | tee N_TERMS
  >>>

  output {
    String date  = read_string("DATE")
    Int n_terms  = read_int("N_TERMS")
    File go_obo  = "go-basic.obo"
  }

  runtime {
    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " HDD"
    disk:        disk_size + " GB" # TES
    preemptible: 1
    maxRetries:  3
  }
}
