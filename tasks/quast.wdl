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
    }

    parameter_meta {
        assembly:               "FASTA file containing assembled genome to check assembly statistics of"
        sample_name:            "Some identifier for naming outputs"
        reference_genome:       "Reference genome to produce comparative statistcis of supplied assembly"
        reference_annotation:   "Annotation of reference genome features for detailed comparative statistics"
        genome_type:            "Taxonomic domain of supplied assembly (default = prokaryote)"
        min_contig:             "Minimum size of contigs to be assessed by tool, in bp (default = 500)"
        est_ref_size:           "Approximate size of supplied reference genome, in bp"
        cpu:                    "Number of CPUs deletaged to task (default = 8)"
        mem_gb:                 "Amount of memory in GB deletaged to task (default = 8)"
        disk_gb:                "Amount of disk space in GB deletaged to task (default = 50)"
    }

    command <<<
        set -euo pipefail

        # Set genome type flag
        GENOME_TYPE_FLAG=""
        if [ "~{genome_type}" = "eukaryote" ]; then GENOME_TYPE_FLAG="--eukaryote"; fi
        if [ "~{genome_type}" = "fungus" ]; then GENOME_TYPE_FLAG="--fungus"; fi

        # Quast command
        quast \
            ~{assembly} \
            ~{"-r " + reference_genome} \
            ~{"--features " + reference_annotation} \
            --min-contig ~{min_contig} \
            ~{"--est-ref-size " + est_ref_size} \
            "${GENOME_TYPE_FLAG}" \
            -o ~{sample_name}_quast \
            -t ~{cpu}

        # Emit tool version
        quast --version | tee VERSION

    >>>

    output {
        File   report_txt     = "~{sample_name}_quast/report.txt"
        File   report_tsv     = "~{sample_name}_quast/report.tsv"
        String quast_version  = read_string("VERSION")
    }

    runtime {
        docker:         "staphb/quast:5.2.0"
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}