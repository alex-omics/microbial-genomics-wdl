version 1.0

task busco {
    input {
        File    assembly
        String  sample_name
        String  busco_lineage   = "auto"
        Int     cpu             = 8
        Int     mem_gb          = 16
        Int     disk_gb         = 100
    }

    parameter_meta {
        assembly:       "FASTA file containing assembled genome to assess completeness of"
        sample_name:    "Some identifier for naming outputs"
        cpu:            "Number of CPUs deletaged to task (default = 8)"
        mem_gb:         "Amount of memory in GB deletaged to task (default = 16)"
        disk_gb:        "Amount of disk space in GB deletaged to task (default = 100)"
    }

    command <<<
        set -euo pipefail

        # Determine run mode: Select specific lineage, or allow to auto-detect
        if [ "~{busco_lineage}" = "auto" ]; then
            echo "No lineage provided; running auto-lineage-prok mode..."
            LINEAGE_ARGS="--auto-lineage-prok"

        else
            echo "Running BUSCO with provided lineage DB: ~{busco_lineage}..."

            # Verify the lineage exists in the image before running
            if [ ! -d "/data/busco_downloads/lineages/~{busco_lineage}" ]; then
                echo "ERROR: lineage '~{busco_lineage}' not found in image." >&2
                echo "Available lineages:" >&2
                ls /data/busco_downloads/lineages/ >&2
                exit 1
            fi

            LINEAGE_ARGS="--lineage_dataset ~{busco_lineage} --offline"
        fi

        # Command execution
        busco \
            --in ~{assembly} \
            --out ~{sample_name} \
            --mode genome \
            --cpu ~{cpu} \
            --download_path /data/busco_downloads \
            "${LINEAGE_ARGS}"

        # Copy summary file to working directory for output capture
        cp ~{sample_name}_busco/short_summary*.txt \
           ~{sample_name}_busco_summary.txt

        # Extract the one-line result string for output
        grep "C:" ~{sample_name}_busco_summary.txt | tr -d '\t' > BUSCO_RESULT

        # Record lineage actually used (for when auto-detect mode is launched)
        grep "Lineage dataset used: ~{sample_name}_busco_summary.txt" | awk '{print $NF}' > LINEAGE_USED

        # Emit tool version too
        busco --version | tee VERSION
    >>>

    output {
        File    summary_txt     =   "~{sample_name}_busco_summary.txt"
        String busco_result     =   read_string("BUSCO_RESULT")
        String lineage_used     =   read_string("LINEAGE_USED")
        String busco_version    =   read_string("VERSION")
    }

    runtime {
        docker:         "aarvani1/busco-prokaryota:5.7.1"
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disk:           "local_disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}