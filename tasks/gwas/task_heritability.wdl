version 1.0

# task_heritability.wdl
#
# Narrow-sense heritability estimation, kinship-based and lineage-based, with
# ALBI confidence intervals.
#
# Translated from microGWAS (MIT-licensed, Copyright (c) 2022, Marco Galardini),
# pinned at commit 1307250. Rules covered: lineages2covariance,
# run_heritability, combine_heritability.
#
# Carried-forward upstream caveat: the heritability point estimate assumes
# normally distributed errors. The microGWAS maintainers themselves flag this
# as potentially inappropriate for binary phenotypes - which is exactly what
# the Tier 2 E. coli virulence phenotype is. Treat these numbers as indicative,
# not as a validation criterion.

task lineages2covariance {

  meta {
    version: "0.1.0"
    description: "Convert the per-sample lineage assignments into a variance/covariance matrix for lineage-based heritability. Translated from the microGWAS `lineages2covariance` rule."
  }

  input {
    String phenotype
    File lineages

    Int memory    = 8
    Int cpu       = 1
    Int disk_size = 20
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    # One-hot encodes lineage membership, then takes the covariance of the
    # transposed indicator matrix
    python3 "$SCRIPTS/lineage2covar.py" ~{lineages} > lineages_covariance.tsv
  >>>

  output {
    String date = read_string("DATE")
    File lineages_covariance = "lineages_covariance.tsv"
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

task run_heritability {

  meta {
    version: "0.1.0"
    description: "Estimate narrow-sense heritability from the kinship matrix and from lineage covariance, with ALBI confidence intervals. Translated from the microGWAS `run_heritability` rule."
  }

  input {
    String phenotype

    File phenotypes
    File similarity
    File lineages_covariance

    Int memory    = 32
    Int cpu       = 2
    Int disk_size = 100
    String docker_image = "aarvani1/microgwas-limix:3.0.4"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"
    ALBI="${MICROGWAS_ALBI:-/opt/microgwas/albi}"

    # --- Kinship-based -------------------------------------------------------
    python3 "$SCRIPTS/estimate_heritability.py" ~{phenotypes} ~{similarity} \
      -p ~{phenotype} \
      > heritability.tsv

    tempfile=$(mktemp -u ./fiesta_kinship)
    python3 "$SCRIPTS/prepare_fiesta.py" ~{phenotypes} ~{similarity} \
      -p ~{phenotype} \
      --prefix "$tempfile"

    # ALBI takes a *pair* of estimates delimiting the interval to search; the
    # rule writes the same normal-model point estimate twice, which upstream
    # relies on. Preserved as-is.
    grep normal heritability.tsv | awk '{print $3}'  > "$tempfile.estimates.txt"
    grep normal heritability.tsv | awk '{print $3}' >> "$tempfile.estimates.txt"

    # `|| true`: ALBI legitimately fails to converge on some phenotypes, and
    # upstream treats an empty CI as acceptable rather than fatal
    python3 "$ALBI/albi.py" \
      -k "${tempfile}_values.txt" \
      -f "$tempfile.estimates.txt" \
      | grep -v "Estimating" | grep -v '#' > heritability.ci.tsv || true

    # --- Lineage-based -------------------------------------------------------
    python3 "$SCRIPTS/estimate_heritability.py" ~{phenotypes} ~{lineages_covariance} \
      -p ~{phenotype} \
      > heritability_lineages.tsv

    tempfile2=$(mktemp -u ./fiesta_lineages)
    python3 "$SCRIPTS/prepare_fiesta.py" ~{phenotypes} ~{lineages_covariance} \
      -p ~{phenotype} \
      --prefix "$tempfile2"

    grep normal heritability_lineages.tsv | awk '{print $3}'  > "$tempfile2.estimates.txt"
    grep normal heritability_lineages.tsv | awk '{print $3}' >> "$tempfile2.estimates.txt"

    python3 "$ALBI/albi.py" \
      -k "${tempfile2}_values.txt" \
      -f "$tempfile2.estimates.txt" \
      | grep -v "Estimat" | grep -v '#' >> heritability.ci.tsv || true

    # heritability.ci.tsv must exist even when both ALBI runs bailed
    touch heritability.ci.tsv
  >>>

  output {
    String date = read_string("DATE")
    File heritability           = "heritability.tsv"
    File heritability_ci        = "heritability.ci.tsv"
    File heritability_lineages  = "heritability_lineages.tsv"
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

task combine_heritability {

  meta {
    version: "0.1.0"
    description: "Merge the kinship, lineage and confidence-interval heritability estimates into one table. Translated from the microGWAS `combine_heritability` rule."
  }

  input {
    String phenotype

    File heritability
    File heritability_lineages
    File heritability_ci

    Int memory    = 4
    Int cpu       = 1
    Int disk_size = 10
    String docker_image = "aarvani1/microgwas-base:0.9.1"
  }

  command <<<
    set -euxo pipefail

    date | tee DATE

    SCRIPTS="${MICROGWAS_SCRIPTS:-/opt/microgwas/workflow/scripts}"

    python3 "$SCRIPTS/combine_heritability.py" \
      ~{heritability} \
      ~{heritability_lineages} \
      ~{heritability_ci} \
      > heritability_all.tsv

    cat heritability_all.tsv
  >>>

  output {
    String date = read_string("DATE")
    File heritability_all = "heritability_all.tsv"
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
