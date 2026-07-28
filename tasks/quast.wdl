version 1.0

task quast {

    input {
        File    assembly
        String  sample_name
        File?   reference_genome
        File?   reference_annotation
        String  genome_type    = "prokaryote"
        Int     min_contig     = 500
        Int?    est_ref_size
        Int     cpu            = 8
        Int     mem_gb         = 8
        Int     disk_gb        = 50
        String  docker         = "staphb/quast:5.2.0"
    }

    parameter_meta {
        assembly:               "FASTA file containing assembled genome to check assembly statistics of"
        sample_name:            "Some identifier for naming outputs"
        reference_genome:       "Reference genome to produce comparative statistics of supplied assembly"
        reference_annotation:   "Annotation of reference genome features for detailed comparative statistics"
        genome_type:            "Taxonomic domain of supplied assembly: prokaryote, eukaryote, or fungus (default = prokaryote)"
        min_contig:             "Minimum size of contigs to be assessed by tool, in bp (default = 500)"
        est_ref_size:           "Approximate size of supplied reference genome, in bp"
        cpu:                    "Number of CPUs delegated to task (default = 8)"
        mem_gb:                 "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:                "Amount of disk space in GB delegated to task (default = 50)"
        docker:                 "Container image"
    }

    command <<<
        set -euo pipefail

        # Build the domain flag as an array. A quoted scalar would expand to an
        # empty positional argument for prokaryotes, which QUAST reads as an
        # extra (empty) input filename and rejects.
        GENOME_TYPE_FLAG=()
        case "~{genome_type}" in
            eukaryote)  GENOME_TYPE_FLAG=(--eukaryote) ;;
            fungus)     GENOME_TYPE_FLAG=(--fungus) ;;
            prokaryote) ;;  # QUAST's default; no flag needed
            *)
                echo "ERROR: unrecognized genome_type '~{genome_type}'" >&2
                echo "Expected one of: prokaryote, eukaryote, fungus" >&2
                exit 1
                ;;
        esac

        # The staphb image exposes quast.py on PATH; there is no bare `quast`.
        quast.py \
            ~{assembly} \
            ~{"-r " + reference_genome} \
            ~{"--features " + reference_annotation} \
            --min-contig ~{min_contig} \
            ~{"--est-ref-size " + est_ref_size} \
            "${GENOME_TYPE_FLAG[@]}" \
            -o ~{sample_name}_quast \
            -t ~{cpu}

        # Lift the headline assembly stats out of report.tsv so they land in the
        # merged summary as real values instead of a file the user has to open.
        # report.tsv is two columns: metric name <tab> value.
        pull() {
            awk -F'\t' -v key="$1" '$1==key {print $2; found=1} END {if (!found) print "NA"}' \
                ~{sample_name}_quast/report.tsv | head -n1
        }

        pull "# contigs"                  > N_CONTIGS
        pull "Total length"               > TOTAL_LENGTH
        pull "Largest contig"             > LARGEST_CONTIG
        pull "N50"                        > N50
        pull "L50"                        > L50
        pull "GC (%)"                     > GC_PERCENT
        pull "# N's per 100 kbp"          > NS_PER_100KBP

        quast.py --version | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1 > VERSION
    >>>

    output {
        File    report_txt      = "~{sample_name}_quast/report.txt"
        File    report_tsv      = "~{sample_name}_quast/report.tsv"
        File    report_html     = "~{sample_name}_quast/report.html"
        Int     n_contigs       = read_int("N_CONTIGS")
        Int     total_length    = read_int("TOTAL_LENGTH")
        Int     largest_contig  = read_int("LARGEST_CONTIG")
        Int     n50             = read_int("N50")
        Int     l50             = read_int("L50")
        Float   gc_percent      = read_float("GC_PERCENT")
        String  ns_per_100kbp   = read_string("NS_PER_100KBP")
        String  quast_version   = read_string("VERSION")
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
