version 1.0

task mark_duplicates {

    input {
        File     input_bam
        File     input_bai
        String   sample_name
        Boolean  remove_duplicates = false
        Int      cpu               = 2
        Int      mem_gb            = 16
        Int      disk_gb           = 50
        String   docker            = "quay.io/biocontainers/picard@sha256:4da4010d1a328f6dec66cf4ccc087480f13d0a96ae21b3bba63169ac9f6164fb"
    }

    parameter_meta {
        input_bam:         "Coordinate-sorted BAM"
        input_bai:         "Index for input_bam"
        sample_name:       "Sample identifier for naming outputs"
        remove_duplicates: "Remove rather than flag duplicates (default = false). Not recommended for RNA-seq."
        cpu:               "Number of CPUs delegated to task (default = 2)"
        mem_gb:            "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:           "Amount of disk space in GB delegated to task (default = 50)"
        docker:            "Container image"
    }

    meta {
        description: "Flag (or optionally remove) PCR/optical duplicates with Picard MarkDuplicates."
    }

    command <<<
        set -euo pipefail

        # The image's picard wrapper is happiest with an explicit heap.
        picard -Xmx~{mem_gb - 2}g MarkDuplicates \
            I=~{input_bam} \
            O="~{sample_name}.markdup.bam" \
            M="~{sample_name}.dup_metrics.txt" \
            REMOVE_DUPLICATES=~{remove_duplicates} \
            ASSUME_SORT_ORDER=coordinate \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT

        # Picard writes X.bai; every other tool in this repo expects X.bam.bai.
        if [ -f "~{sample_name}.markdup.bai" ]; then
            mv "~{sample_name}.markdup.bai" "~{sample_name}.markdup.bam.bai"
        fi

        # PERCENT_DUPLICATION lives in the column-headed row of the metrics file.
        awk -F'\t' '/^LIBRARY/{for(i=1;i<=NF;i++) if($i=="PERCENT_DUPLICATION") c=i; next} c && NF>=c && $1!="" && !/^#/{print $c; exit}' \
            "~{sample_name}.dup_metrics.txt" > PCT_DUP
        [ -s PCT_DUP ] || echo "NA" > PCT_DUP
    >>>

    output {
        File   markdup_bam  = "~{sample_name}.markdup.bam"
        File   markdup_bai  = "~{sample_name}.markdup.bam.bai"
        File   metrics_file = "~{sample_name}.dup_metrics.txt"
        String pct_duplication = read_string("PCT_DUP")
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
