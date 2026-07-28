version 1.0

# task_rare_variants.wdl
#
# Reference-based variant calling and the rare/common variant split.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: get_snps, prepare_vcf_variants.
#
# Scope note: upstream gates the *rare* variant set behind Sequence UNET
# deleteriousness predictions (vcf2deleterious.py hard-errors without either
# --unet or --sift). Sequence UNET is out of scope for Milestone 1 - the
# microGWAS paper itself reports the gene burden test as the one modality that
# did not reproduce the earlier validation results. So `deleteriousness_dir` is
# optional here: without it this task still produces the *common* variant set,
# which the post-publication microGWAS added as its 6th association approach
# and which Tier 2 does exercise.

task get_snps {

  meta {
    version: "0.1.0"
    description: "Call variants for one assembly against the reference with snippy. Scattered per sample. Translated from the microGWAS `get_snps` rule."
  }

  input {
    # Sample ID drives snippy's --outdir, which is what snippy writes into the
    # VCF's sample column - so it must be the bare strain name, not a path.
    String samplename
    File assembly_fasta
    File reference_gbk

    Int snippy_ram = 8

    Int memory    = 16
    Int cpu       = 2
    Int disk_size = 50
    String docker_image = "aarvani1/microgwas-snippy:4.6.0"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    snippy --version 2>&1 | tee VERSION

    # --outdir basename becomes the VCF sample name; keep it as the strain ID
    # so `bcftools merge` downstream produces correctly-labelled columns.
    snippy \
      --force \
      --outdir "~{samplename}" \
      --ref ~{reference_gbk} \
      --ctgs ~{assembly_fasta} \
      --cpus ~{cpu} \
      --ram ~{snippy_ram}

    # Confirm the sample really is labelled with the strain ID before this
    # propagates into a merged multi-sample VCF where it is hard to unpick
    bcftools query -l "~{samplename}/snps.vcf.gz" | tee VCF_SAMPLE
  >>>

  output {
    String date           = read_string("DATE")
    String snippy_version = read_string("VERSION")
    String vcf_sample     = read_string("VCF_SAMPLE")

    File snps_vcf       = "~{samplename}/snps.vcf.gz"
    File snps_vcf_index = "~{samplename}/snps.vcf.gz.csi"
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

task prepare_vcf_variants {

  meta {
    version: "0.1.0"
    description: "Merge per-sample VCFs and split them into rare (deleterious, <5% AF) and common (>1% AF) variant sets. Translated from the microGWAS `prepare_vcf_variants` rule."
  }

  input {
    Array[File] snps_vcfs
    Array[File] snps_vcf_indexes

    # Sequence UNET per-protein deleteriousness predictions, as a tarball of
    # the directory vcf2deleterious.py expects. Omit to skip the rare-variant
    # set entirely (see file header).
    File? deleteriousness_dir_tarball

    # bcftools view -Q: keep sites *below* this AF for the rare set
    Float rare_max_af = 0.05
    # bcftools view -q: keep sites *above* this AF for the common set
    Float common_min_af = 0.01

    Int memory    = 32
    Int cpu       = 16
    Int disk_size = 200
    String docker_image = "aarvani1/microgwas-bcftools:1.13"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE
    bcftools --version | head -1 | tee VERSION

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # Co-locate each VCF with its index. Cromwell localizes array members into
    # separate shard directories, so a bare list of VCF paths would leave
    # bcftools unable to find the .csi sitting next to each one.
    mkdir -p vcfs
    vcfs=(~{sep=' ' snps_vcfs})
    idxs=(~{sep=' ' snps_vcf_indexes})
    if [[ ${#vcfs[@]} -ne ${#idxs[@]} ]]; then
      echo "ERROR: VCF and index arrays differ in length" >&2
      exit 1
    fi
    : > bcftools_input.txt
    for i in "${!vcfs[@]}"; do
      # Name by VCF sample so collisions surface immediately
      s=$(bcftools query -l "${vcfs[$i]}" | head -1)
      ln -s "${vcfs[$i]}" "vcfs/${s}.vcf.gz"
      ln -s "${idxs[$i]}" "vcfs/${s}.vcf.gz.csi"
      echo "vcfs/${s}.vcf.gz" >> bcftools_input.txt
    done
    echo "Samples to merge: $(wc -l < bcftools_input.txt)"

    # -0 treats missing genotypes as reference, which is the assumption the
    # downstream AF filters are built on
    bcftools merge -l bcftools_input.txt -0 -O z --threads ~{cpu} > merged.vcf.gz
    bcftools norm -m - -O z --threads ~{cpu} merged.vcf.gz > norm.vcf.gz

    # Common set: sites above the minimum AF
    bcftools view -q ~{common_min_af} norm.vcf.gz -O z --threads ~{cpu} > common.vcf.gz
    bcftools index common.vcf.gz

    ~{if defined(deleteriousness_dir_tarball) then "RUN_RARE=1" else "RUN_RARE=0"}

    if [[ "$RUN_RARE" == "1" ]]; then
      mkdir -p unet
      tar -xf ~{default="/dev/null" deleteriousness_dir_tarball} -C unet --strip-components=1

      # Rare set: sites below the maximum AF, then filtered to those predicted
      # deleterious (UNET score >= 0.5)
      bcftools view -Q ~{rare_max_af} norm.vcf.gz -O z --threads ~{cpu} > filtered.vcf.gz
      python3 "$SCRIPTS/vcf2deleterious.py" filtered.vcf.gz --unet unet | bgzip > rare.vcf.gz
      bcftools index rare.vcf.gz
      rm -f filtered.vcf.gz
    else
      echo "No deleteriousness predictions supplied - skipping the rare variant set." | tee RARE_SKIPPED
    fi

    rm -f merged.vcf.gz norm.vcf.gz

    bcftools index -n common.vcf.gz | tee N_COMMON
  >>>

  output {
    String date     = read_string("DATE")
    String bcftools_version = read_string("VERSION")
    Int n_common_variants = read_int("N_COMMON")

    File common_snps       = "common.vcf.gz"
    File common_snps_index = "common.vcf.gz.csi"

    File? rare_snps        = "rare.vcf.gz"
    File? rare_snps_index  = "rare.vcf.gz.csi"
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
