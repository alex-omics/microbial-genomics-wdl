version 1.0

# Small, workflow-agnostic helpers shared across workflows. Anything here must
# stay generic: if a helper starts needing to know what it is summarising, it
# belongs in a tool-specific task file instead.

task name_summary {

    input {
        Array[Array[String]] rows
        String               basename
        String               docker = "ubuntu:22.04"
    }

    parameter_meta {
        rows:     "Table contents including the header row"
        basename: "Basename for the emitted TSV, without extension"
        docker:   "Container image"
    }

    meta {
        description: "Write a table to a recognisably named file. write_tsv alone emits a temp-named file that lands in the bucket as something like 'tmpzsaunsvb', so routing it through a task is what gives the deliverable a name you can find."
    }

    command <<<
        set -euo pipefail
        cp ~{write_tsv(rows)} "~{basename}.tsv"
    >>>

    output {
        File summary = "~{basename}.tsv"
    }

    runtime {
        docker:         docker
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     2
    }
}


task concat_tables {

    input {
        Array[File]  tables
        String       basename
        Int          disk_gb = 50
        Int          mem_gb  = 8
        String       docker  = "ubuntu:22.04"
    }

    parameter_meta {
        tables:   "TSVs sharing an identical header. Each file's header is dropped except the first."
        basename: "Basename for the merged TSV, without extension"
        disk_gb:  "Amount of disk space in GB delegated to task (default = 50)"
        mem_gb:   "Amount of memory in GB delegated to task (default = 8)"
        docker:   "Container image"
    }

    meta {
        description: "Concatenate same-schema TSVs into one long-format table, keeping a single header."
    }

    command <<<
        set -euo pipefail

        head -n1 ~{tables[0]} > "~{basename}.tsv"

        # Guard against silently stitching together tables that do not share a
        # schema; the result would look fine and be wrong.
        EXPECTED="$(head -n1 ~{tables[0]})"
        for f in ~{sep=' ' tables}; do
            if [ "$(head -n1 "${f}")" != "${EXPECTED}" ]; then
                echo "ERROR: header mismatch in ${f}" >&2
                echo "  expected: ${EXPECTED}" >&2
                echo "  found:    $(head -n1 "${f}")" >&2
                exit 1
            fi
            awk 'NR>1' "${f}" >> "~{basename}.tsv"
        done

        gzip -c "~{basename}.tsv" > "~{basename}.tsv.gz"
        awk 'NR>1' "~{basename}.tsv" | wc -l > N_ROWS
    >>>

    output {
        File merged     = "~{basename}.tsv"
        File merged_gz  = "~{basename}.tsv.gz"
        Int  n_rows     = read_int("N_ROWS")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            1
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}
