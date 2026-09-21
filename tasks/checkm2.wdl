version 1.0

task checkm2 {

    input {
        File    assembly
        String  sample_name
        File?   checkm2_db
        Boolean force_general_model = false
        Boolean low_memory          = false
        Int     cpu                 = 8
        Int     mem_gb              = 16
        Int     disk_gb             = 50
        String  docker              = "staphb/checkm2:1.1.0"
    }

    parameter_meta {
        assembly:            "FASTA file containing assembled genome to assess completeness and contamination of"
        sample_name:         "Some identifier for naming outputs"
        checkm2_db:          "CheckM2 DIAMOND database (uniref100.KO.1.dmnd, ~3 GB). Recommended: stage it once in your bucket and pass it here. Otherwise the task downloads it at runtime, which is slower and depends on an external host."
        force_general_model: "Force the general (gradient boost) model instead of letting CheckM2 choose (default = false)"
        low_memory:          "Reduce DIAMOND block size to cut RAM use at the cost of runtime (default = false)"
        cpu:                 "Number of CPUs delegated to task (default = 8)"
        mem_gb:              "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:             "Amount of disk space in GB delegated to task (default = 50)"
        docker:              "Container image"
    }

    command <<<
        set -euo pipefail

        # staphb/checkm2 ships no DIAMOND database, so one must be supplied or
        # fetched. Passing it as a File is the Terra-native path: GCS -> VM
        # localization is fast and avoids a ~3 GB image.
        DB_PATH="~{select_first([checkm2_db, ''])}"
        if [ -n "${DB_PATH}" ]; then
            echo "Using supplied CheckM2 database: ${DB_PATH}"
        else
            echo "WARNING: no checkm2_db supplied; downloading at runtime." >&2
            checkm2 database --download --path ./checkm2_db
            DB_PATH="$(find ./checkm2_db -name '*.dmnd' | head -n1)"
            if [ -z "${DB_PATH}" ]; then
                echo "ERROR: database download produced no .dmnd file" >&2
                exit 1
            fi
        fi

        # CheckM2 scans a directory and keys results off the filename, so stage
        # the assembly under a name we control. Decompress if needed; CheckM2
        # does not read gzipped FASTA.
        mkdir -p input_bins
        case "~{assembly}" in
            *.gz) gunzip -c ~{assembly} > "input_bins/~{sample_name}.fna" ;;
            *)    cp ~{assembly} "input_bins/~{sample_name}.fna" ;;
        esac

        EXTRA_ARGS=()
        ~{if force_general_model then "EXTRA_ARGS+=(--general)" else ""}
        ~{if low_memory then "EXTRA_ARGS+=(--lowmem)" else ""}

        checkm2 predict \
            --input input_bins \
            --output-directory ~{sample_name}_checkm2 \
            --extension .fna \
            --threads ~{cpu} \
            --database_path "${DB_PATH}" \
            --force \
            "${EXTRA_ARGS[@]}"

        cp ~{sample_name}_checkm2/quality_report.tsv ~{sample_name}_checkm2_report.tsv

        # quality_report.tsv is a header row plus one row per bin. We run one
        # assembly per task, so read row 2 by column name.
        pull() {
            awk -F'\t' -v key="$1" '
                NR==1 { for (i=1; i<=NF; i++) if ($i==key) col=i; next }
                NR==2 { print (col ? $col : "NA"); found=1 }
                END   { if (!found) print "NA" }
            ' ~{sample_name}_checkm2_report.tsv
        }

        # CheckM2 names the completeness column after the model mode: 'Completeness'
        # in auto mode, but 'Completeness_General' under --general and
        # 'Completeness_Specific' under --specific. Take whichever is present, or
        # read_float would fail the task whenever a model is forced.
        pull_first() {
            local v
            for key in "$@"; do
                v="$(pull "${key}")"
                if [ "${v}" != "NA" ]; then
                    echo "${v}"
                    return 0
                fi
            done
            echo "NA"
        }

        pull_first "Completeness" "Completeness_General" "Completeness_Specific" > COMPLETENESS
        pull "Contamination"                                                     > CONTAMINATION
        # Absent when a model is forced; fall back to the mode we asked for.
        pull_first "Completeness_Model_Used" > MODEL_USED
        if [ "$(cat MODEL_USED)" = "NA" ]; then
            ~{if force_general_model then "echo general" else "echo auto"} > MODEL_USED
        fi
        pull "Coding_Density"           > CODING_DENSITY
        pull "Contig_N50"               > CONTIG_N50
        pull "Genome_Size"              > GENOME_SIZE
        pull "GC_Content"               > GC_CONTENT
        pull "Total_Coding_Sequences"   > TOTAL_CDS

        checkm2 --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1 > VERSION
    >>>

    output {
        File    report_tsv          = "~{sample_name}_checkm2_report.tsv"
        Float   completeness        = read_float("COMPLETENESS")
        Float   contamination       = read_float("CONTAMINATION")
        String  model_used          = read_string("MODEL_USED")
        String  coding_density      = read_string("CODING_DENSITY")
        String  contig_n50          = read_string("CONTIG_N50")
        String  genome_size         = read_string("GENOME_SIZE")
        String  gc_content          = read_string("GC_CONTENT")
        String  total_cds           = read_string("TOTAL_CDS")
        String  checkm2_version     = read_string("VERSION")
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
