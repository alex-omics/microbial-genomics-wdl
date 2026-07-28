version 1.0

task busco {
    input {
        File    assembly
        String  sample_name
        String  busco_lineage   = "auto"
        File?   lineage_tarball
        Int     cpu             = 8
        Int     mem_gb          = 16
        Int     disk_gb         = 100
        Int     boot_disk_gb    = 30
        String  docker          = "aarvani1/busco-prokaryota:5.8.0"
    }

    parameter_meta {
        assembly:        "FASTA file containing assembled genome to assess completeness of"
        sample_name:     "Some identifier for naming outputs"
        busco_lineage:   "Lineage dataset name (e.g. bacteria_odb10, spirochaetales_odb10), or 'auto' to let BUSCO place the assembly in the prokaryote tree (default = auto)"
        lineage_tarball: "Optional .tar.gz of a lineage dataset not baked into the image; overrides busco_lineage"
        cpu:             "Number of CPUs delegated to task (default = 8)"
        mem_gb:          "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:         "Amount of disk space in GB delegated to task (default = 100)"
        boot_disk_gb:    "Boot disk size in GB. Must exceed the unpacked image (~10 GB); Cromwell's 10 GB default is too small (default = 30)"
        docker:          "Container image; must carry an offline /busco_downloads tree"
    }

    command <<<
        set -euo pipefail

        # The image consolidates lineages, placement files, and file_versions.tsv
        # under a single path. BUSCO's offline resolver needs all three together.
        BUSCO_DL=/busco_downloads

        if [ ! -f "${BUSCO_DL}/file_versions.tsv" ]; then
            echo "ERROR: ${BUSCO_DL}/file_versions.tsv missing; image is not offline-ready." >&2
            exit 1
        fi

        # A user-supplied dataset wins over the baked-in set, so lineages we
        # haven't packaged (odb12, eukaryotes) can still be run without a rebuild.
        LINEAGE="~{busco_lineage}"
        LINEAGE_TARBALL="~{select_first([lineage_tarball, ''])}"
        if [ -n "${LINEAGE_TARBALL}" ]; then
            echo "Installing user-supplied lineage from ${LINEAGE_TARBALL}"
            tar -xzf "${LINEAGE_TARBALL}" -C "${BUSCO_DL}/lineages/"
            # Dataset name is the top-level directory inside the archive.
            LINEAGE="$(tar -tzf "${LINEAGE_TARBALL}" | head -n1 | cut -d/ -f1)"
            echo "Resolved supplied lineage to: ${LINEAGE}"
        fi

        if [ "${LINEAGE}" = "auto" ]; then
            echo "Running BUSCO in offline auto-lineage mode over the prokaryote tree..."
            LINEAGE_ARGS="--auto-lineage-prok"
        else
            echo "Running BUSCO with lineage dataset: ${LINEAGE}"
            if [ ! -d "${BUSCO_DL}/lineages/${LINEAGE}" ]; then
                echo "ERROR: lineage '${LINEAGE}' not found in image." >&2
                echo "Available lineages:" >&2
                ls "${BUSCO_DL}/lineages/" >&2
                exit 1
            fi
            LINEAGE_ARGS="--lineage_dataset ${LINEAGE}"
        fi

        OUTDIR="~{sample_name}_busco"

        # LINEAGE_ARGS is deliberately unquoted: it carries multiple words.
        busco \
            --in ~{assembly} \
            --out "${OUTDIR}" \
            --out_path . \
            --mode genome \
            --cpu ~{cpu} \
            --download_path "${BUSCO_DL}" \
            --offline \
            ${LINEAGE_ARGS}

        # Auto-lineage emits one summary per dataset it tried (generic root plus
        # the placed specific lineage). The specific result is the one to keep;
        # fall back to whatever exists for the explicit-lineage case.
        SUMMARY="$(ls -1 "${OUTDIR}"/short_summary.specific.*.txt 2>/dev/null | head -n1 || true)"
        if [ -z "${SUMMARY}" ]; then
            SUMMARY="$(ls -1 "${OUTDIR}"/short_summary.*.txt 2>/dev/null | head -n1 || true)"
        fi
        if [ -z "${SUMMARY}" ]; then
            echo "ERROR: BUSCO produced no short summary in ${OUTDIR}" >&2
            ls -R "${OUTDIR}" >&2
            exit 1
        fi
        cp "${SUMMARY}" ~{sample_name}_busco_summary.txt

        # Which dataset BUSCO actually used — the point of auto-lineage mode.
        # Header line reads: "# The lineage dataset is: NAME (Creation date: ...)"
        grep -m1 'The lineage dataset is:' ~{sample_name}_busco_summary.txt \
            | sed -E 's/.*The lineage dataset is:[[:space:]]*([^[:space:]]+).*/\1/' \
            > LINEAGE_USED
        USED="$(cat LINEAGE_USED)"

        # Keep the per-marker table too; it's what you want when a genome looks
        # incomplete and you need to know which BUSCOs are actually missing.
        # Auto-lineage leaves a run_ dir for every dataset it tried, so select the
        # one matching the dataset actually reported — not whatever sorts first,
        # which is the generic root run.
        FULL="${OUTDIR}/run_${USED}/full_table.tsv"
        if [ ! -f "${FULL}" ]; then
            FULL="$(find "${OUTDIR}" -name full_table.tsv | head -n1 || true)"
        fi
        if [ -n "${FULL}" ] && [ -f "${FULL}" ]; then
            cp "${FULL}" ~{sample_name}_busco_full_table.tsv
        else
            echo "# no full_table.tsv produced" > ~{sample_name}_busco_full_table.tsv
        fi

        # The one-line BUSCO notation, e.g. C:99.2%[S:98.4%,D:0.8%],F:0.0%,M:0.8%,n:124
        grep -m1 -oE 'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+' \
            ~{sample_name}_busco_summary.txt > BUSCO_RESULT
        RESULT="$(cat BUSCO_RESULT)"

        # Break the notation into individually addressable outputs so the summary
        # table carries real numbers rather than a string to re-parse downstream.
        # Each key occurs exactly once, so grep the key and strip it; sed word
        # boundaries (\b) are not portable ERE and silently pass input through.
        parse_metric() {
            echo "${RESULT}" | grep -oE "$1[0-9.]+" | head -n1 | sed "s/^$1//"
        }
        parse_metric 'C:' > COMPLETE_PCT
        parse_metric 'S:' > SINGLE_PCT
        parse_metric 'D:' > DUP_PCT
        parse_metric 'F:' > FRAG_PCT
        parse_metric 'M:' > MISSING_PCT
        parse_metric 'n:' > N_MARKERS

        busco --version | awk '{print $NF}' > VERSION
    >>>

    output {
        File    summary_txt         = "~{sample_name}_busco_summary.txt"
        File    full_table_tsv      = "~{sample_name}_busco_full_table.tsv"
        String  busco_result        = read_string("BUSCO_RESULT")
        String  lineage_used        = read_string("LINEAGE_USED")
        Float   complete_pct        = read_float("COMPLETE_PCT")
        Float   single_copy_pct     = read_float("SINGLE_PCT")
        Float   duplicated_pct      = read_float("DUP_PCT")
        Float   fragmented_pct      = read_float("FRAG_PCT")
        Float   missing_pct         = read_float("MISSING_PCT")
        Int     n_markers           = read_int("N_MARKERS")
        String  busco_version       = read_string("VERSION")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        # The offline lineage set makes this image ~10 GB unpacked, and the image
        # lands on the boot disk, not local-disk. Cromwell's 10 GB boot default is
        # not enough and the pull fails before the task ever starts.
        bootDiskSizeGb: boot_disk_gb
        preemptible:    1
        maxRetries:     2
    }
}
