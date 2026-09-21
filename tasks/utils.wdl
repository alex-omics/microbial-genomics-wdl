version 1.0

# Small, workflow-agnostic helpers shared across workflows. Anything here must
# stay generic: if a helper starts needing to know what it is summarising, it
# belongs in a tool-specific task file instead.

task name_summary {

    input {
        Array[Array[String]] rows
        String               basename
        String               docker = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        rows:     "Table contents including the header row"
        basename: "Basename for the emitted TSV, without extension"
        docker:   "Container image"
    }

    meta {
        description: "Write a table to a named file. write_tsv alone produces a temp-named file (e.g. 'tmpzsaunsvb') in the bucket; routing it through a task gives the output a findable name."
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


task validate_panel {

    input {
        Array[String]   derived_names
        Int             n_primary
        Array[Int]      companion_counts
        Array[String]?  sample_names
        String          docker = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        derived_names:    "Names derived from the primary input filenames, used when sample_names is not supplied"
        n_primary:        "Length of the primary input array"
        companion_counts: "Lengths of every array that must be positionally matched to the primary one"
        sample_names:     "Optional user-supplied labels"
        docker:           "Container image"
    }

    meta {
        description: "Fail fast on a mis-specified panel and emit sanitised sample IDs. Positional matching across several arrays is easy to get wrong and otherwise surfaces as an index error inside a scatter, after the expensive stages have run. IDs are also checked for characters that would break Bakta locus tags or Panaroo column headers."
    }

    command <<<
        set -euo pipefail

        N=~{n_primary}
        for c in ~{sep=' ' companion_counts}; do
            if [ "${c}" -ne "${N}" ]; then
                echo "ERROR: input arrays are not the same length." >&2
                echo "       primary array has ${N} entries; a companion array has ${c}." >&2
                echo "       modbams, assemblies and sample_names are matched positionally." >&2
                exit 1
            fi
        done

        # write_lines materialises the arrays as files, which keeps shell
        # quoting out of the picture entirely -- names may contain characters
        # that would need escaping if interpolated into a command string.
        SUPPLIED="~{if defined(sample_names) then write_lines(select_first([sample_names, []])) else ''}"
        if [ -n "${SUPPLIED}" ]; then
            cp "${SUPPLIED}" names.txt
        else
            cp ~{write_lines(derived_names)} names.txt
        fi

        if [ "$(wc -l < names.txt)" -ne "${N}" ]; then
            echo "ERROR: resolved $(wc -l < names.txt) sample names for ${N} inputs." >&2
            exit 1
        fi

        # These IDs become Bakta locus-tag prefixes, output filenames, and the
        # isolate column headers Panaroo derives from GFF filenames. Anything
        # outside [A-Za-z0-9._-] breaks at least one of those, usually silently
        # until the ortholog join cannot match a column.
        if grep -qE '[^A-Za-z0-9._-]' names.txt; then
            echo "ERROR: sample names contain characters that will break locus tags" >&2
            echo "       or Panaroo column headers. Use only letters, digits, dot," >&2
            echo "       underscore and hyphen. Offending names:" >&2
            grep -nE '[^A-Za-z0-9._-]' names.txt >&2
            exit 1
        fi

        if [ "$(sort names.txt | uniq -d | wc -l)" -ne 0 ]; then
            echo "ERROR: duplicate sample names would collide in the pangenome:" >&2
            sort names.txt | uniq -d >&2
            exit 1
        fi

        echo "Validated ${N} samples"
        cat names.txt
    >>>

    output {
        Array[String] sample_ids = read_lines("names.txt")
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
        String       docker  = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
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
