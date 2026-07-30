version 1.0

task modkit_pileup {

    input {
        File     aligned_bam
        File     aligned_bam_index
        String   sample_name
        File     reference_fasta
        Boolean  no_filtering    = false
        Float?   filter_threshold
        Int      min_coverage    = 10
        Float    min_percent     = 50.0
        Int      cpu             = 8
        Int      mem_gb          = 16
        Int      disk_gb         = 100
        String   docker          = "quay.io/biocontainers/ont-modkit:0.6.4--h7f49ad2_0@sha256:003ed784097737910c19abbe71fa2dffc32f888d228e5b5b12c4ec2adb799e0f"
    }

    parameter_meta {
        aligned_bam:       "Coordinate-sorted BAM with MM/ML tags intact, aligned to reference_fasta"
        aligned_bam_index: "The .bai for aligned_bam. WDL localizes files individually, so the index must be passed explicitly or modkit will not find it."
        sample_name:       "Some identifier for naming outputs"
        reference_fasta:   "The same reference the BAM was aligned to"
        no_filtering:      "Keep every modification call regardless of confidence. Off by default; modkit otherwise estimates a per-mod confidence threshold by sampling reads (default = false)"
        filter_threshold:  "Override modkit's estimated confidence threshold with a fixed value (0-1). Leave unset to let modkit estimate."
        min_coverage:      "Minimum Nvalid_cov for a site to count as confidently observed in the summary stats (default = 10)"
        min_percent:       "Minimum percent-modified for a site to be counted as methylated in the summary stats (default = 50.0)"
        cpu:               "Number of CPUs delegated to task (default = 8)"
        mem_gb:            "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:           "Amount of disk space in GB delegated to task (default = 100)"
        docker:            "Container image"
    }

    command <<<
        set -euo pipefail

        # modkit reads the index from alongside the BAM. Localization can land
        # the two in separate directories, so stage them together.
        mkdir -p bam
        cp ~{aligned_bam}       bam/~{sample_name}.bam
        cp ~{aligned_bam_index} bam/~{sample_name}.bam.bai

        EXTRA_ARGS=()
        ~{if no_filtering then "EXTRA_ARGS+=(--no-filtering)" else ""}
        ~{if defined(filter_threshold) then "EXTRA_ARGS+=(--filter-threshold " + filter_threshold + ")" else ""}

        modkit pileup \
            bam/~{sample_name}.bam \
            ~{sample_name}.bedmethyl.bed \
            --ref ~{reference_fasta} \
            --threads ~{cpu} \
            --log-filepath ~{sample_name}_modkit.log \
            "${EXTRA_ARGS[@]}"

        # bedMethyl is one row per (position, strand, modification code). The
        # codes that matter for bacteria are 6mA and 4mC; 5mC in CpG context is
        # the vertebrate story and is largely incidental here. Splitting by code
        # up front keeps every downstream intersect honest about which
        # modification it is actually talking about.
        #
        #   a      6mA
        #   m      5mC
        #   21839  4mC   (ChEBI identifier, not a typo)
        #
        # Columns: 4 = mod code, 10 = Nvalid_cov, 11 = percent modified.
        awk -F'\t' '$4=="a"'     ~{sample_name}.bedmethyl.bed > ~{sample_name}.6mA.bed  || true
        awk -F'\t' '$4=="21839"' ~{sample_name}.bedmethyl.bed > ~{sample_name}.4mC.bed  || true
        awk -F'\t' '$4=="m"'     ~{sample_name}.bedmethyl.bed > ~{sample_name}.5mC.bed  || true

        # Per-code counts of confidently methylated sites. A site is counted
        # only if it clears both coverage and percent-modified floors: a 100%
        # modified position at 3x is noise, not a methylation call.
        count_sites() {
            local f="$1"
            if [ ! -s "${f}" ]; then echo "0"; return; fi
            awk -F'\t' -v mincov=~{min_coverage} -v minpct=~{min_percent} \
                '$10>=mincov && $11>=minpct {n++} END {print n+0}' "${f}"
        }

        mean_pct() {
            local f="$1"
            if [ ! -s "${f}" ]; then echo "NA"; return; fi
            awk -F'\t' -v mincov=~{min_coverage} \
                '$10>=mincov {s+=$11; n++} END {if (n>0) printf "%.2f\n", s/n; else print "NA"}' "${f}"
        }

        count_sites ~{sample_name}.6mA.bed > N_6MA
        count_sites ~{sample_name}.4mC.bed > N_4MC
        count_sites ~{sample_name}.5mC.bed > N_5MC

        mean_pct ~{sample_name}.6mA.bed > PCT_6MA
        mean_pct ~{sample_name}.4mC.bed > PCT_4MC
        mean_pct ~{sample_name}.5mC.bed > PCT_5MC

        # Total positions modkit was able to evaluate at all, which is the
        # denominator for judging whether coverage was adequate genome-wide.
        awk -F'\t' -v mincov=~{min_coverage} '$10>=mincov {n++} END {print n+0}' \
            ~{sample_name}.bedmethyl.bed > N_COVERED

        gzip -c ~{sample_name}.bedmethyl.bed > ~{sample_name}.bedmethyl.bed.gz

        modkit --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1 > VERSION
    >>>

    output {
        File    bedmethyl           = "~{sample_name}.bedmethyl.bed"
        File    bedmethyl_gz        = "~{sample_name}.bedmethyl.bed.gz"
        File    bed_6ma             = "~{sample_name}.6mA.bed"
        File    bed_4mc             = "~{sample_name}.4mC.bed"
        File    bed_5mc             = "~{sample_name}.5mC.bed"
        File    modkit_log          = "~{sample_name}_modkit.log"
        Int     n_sites_6ma         = read_int("N_6MA")
        Int     n_sites_4mc         = read_int("N_4MC")
        Int     n_sites_5mc         = read_int("N_5MC")
        Int     n_positions_covered = read_int("N_COVERED")
        String  mean_percent_6ma    = read_string("PCT_6MA")
        String  mean_percent_4mc    = read_string("PCT_4MC")
        String  mean_percent_5mc    = read_string("PCT_5MC")
        String  modkit_version      = read_string("VERSION")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}


task modkit_find_motifs {

    input {
        File     bedmethyl
        String   sample_name
        File     reference_fasta
        Int      cpu     = 8
        Int      mem_gb  = 32
        Int      disk_gb = 100
        String   docker  = "quay.io/biocontainers/ont-modkit:0.6.4--h7f49ad2_0@sha256:003ed784097737910c19abbe71fa2dffc32f888d228e5b5b12c4ec2adb799e0f"
    }

    parameter_meta {
        bedmethyl:       "Uncompressed bedMethyl from modkit_pileup"
        sample_name:     "Some identifier for naming outputs"
        reference_fasta: "The same reference the pileup was called against"
        cpu:             "Number of CPUs delegated to task (default = 8)"
        mem_gb:          "Amount of memory in GB delegated to task (default = 32)"
        disk_gb:         "Amount of disk space in GB delegated to task (default = 100)"
        docker:          "Container image"
    }

    meta {
        description: "De novo discovery of methylated sequence motifs. In bacteria these correspond to restriction-modification system recognition sites, so the motif inventory is effectively a readout of which MTases the isolate carries — which is the control you need before attributing any cross-isolate methylation difference to regulation."
    }

    command <<<
        set -euo pipefail

        modkit find-motifs \
            --in-bedmethyl ~{bedmethyl} \
            --ref ~{reference_fasta} \
            --out-table ~{sample_name}_motifs.tsv \
            --threads ~{cpu} \
            > ~{sample_name}_find_motifs.log 2>&1 || {
                echo "modkit find-motifs failed; see log" >&2
                cat ~{sample_name}_find_motifs.log >&2
                exit 1
            }

        # One line per discovered motif, tab separated. Collapse to a single
        # comma-joined string so it can ride along in the summary table.
        if [ -s ~{sample_name}_motifs.tsv ]; then
            awk -F'\t' 'NR>1 {printf "%s%s", sep, $1; sep=","} END {print ""}' \
                ~{sample_name}_motifs.tsv > MOTIFS
            awk 'NR>1 {n++} END {print n+0}' ~{sample_name}_motifs.tsv > N_MOTIFS
        else
            echo "NA" > MOTIFS
            echo "0"  > N_MOTIFS
        fi
    >>>

    output {
        File    motifs_tsv  = "~{sample_name}_motifs.tsv"
        File    motifs_log  = "~{sample_name}_find_motifs.log"
        String  motifs       = read_string("MOTIFS")
        Int     n_motifs     = read_int("N_MOTIFS")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}
