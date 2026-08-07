version 1.0

task fasterq_dump {
    input {
        String accession
        Int    disk_gb   = 50
        Int    cpu       = 4
        Int    memory_gb = 8
    }

    parameter_meta {
        accession: "Single-run SRA/ENA/DDBJ accession (SRR/ERR/DRR) — wire to this.SRA_ID. Not SRX/SRS/SRP/PRJ*: those identify an experiment, sample, study, or project, each of which can span multiple runs, and this task fetches exactly one"
        disk_gb:   "Disk space in GB (default=50; fine for MiSeq-scale Illumina, bump substantially for a full ONT/PacBio flowcell)"
        cpu:       "CPUs for fasterq-dump (default=4)"
        memory_gb: "Memory in GB (default=8)"
    }

    command <<<
        set -euo pipefail

        # ── Validate accession shape ──────────────────────────────
        # A container accession (SRX/SRS/SRP/PRJ*) can resolve to many runs;
        # feeding one to prefetch either fails outright or silently grabs only
        # its first run. Fail loud here instead of a confusing error three
        # minutes into prefetch, or a quietly wrong single-run result.
        case "~{accession}" in
            SRR*|ERR*|DRR*) ;;
            *)
                echo "ERROR: '~{accession}' is not a single-run accession (expected SRR/ERR/DRR)." >&2
                echo "This task fetches exactly one run. Experiment/sample/study/project" >&2
                echo "accessions (SRX/SRS/SRP/PRJ*) can span multiple runs — resolve to a" >&2
                echo "run accession first, e.g. via the ENA filereport API." >&2
                exit 1
                ;;
        esac

        # ── Log version ───────────────────────────────────────────
        fasterq-dump --version 2>&1 | head -n2 | tee SRA_TOOLS_VERSION

        # ── Prefetch: cache .sra locally before dumping ───────────
        # More reliable than streaming directly on cloud VMs —
        # prefetch handles retries internally if the connection
        # drops mid-download.
        prefetch ~{accession} --output-directory .

        # ── Dump: split into R1 + R2 ─────────────────────────────
        # --split-3 matches pairs by read name not position, fixing the
        # BWA "paired reads have different names" error seen with --split-files
        # on older SRA submissions. Unpaired reads go to a third file.
        fasterq-dump ~{accession} \
            --split-3 \
            --threads ~{cpu} \
            --outdir . \
            --temp .

        # ── Detect what --split-3 actually produced ──────────────
        # Paired-end runs emit _1/_2; single-end runs (common for ONT/PacBio
        # long-read bacterial data) emit only the unsuffixed file — there is
        # no _2 to gzip in that case. A run nominally paired can also leave a
        # few reads that couldn't be matched in the unsuffixed file alongside
        # _1/_2; keep that rather than leaving it uncollected on disk.
        R1="~{accession}_1.fastq"
        R2="~{accession}_2.fastq"
        UNPAIRED="~{accession}.fastq"

        if [ -f "${R1}" ] && [ -f "${R2}" ]; then
            echo "paired" > LAYOUT
            gzip "${R1}"
            gzip "${R2}"
            if [ -f "${UNPAIRED}" ]; then
                gzip "${UNPAIRED}"
            fi
        elif [ -f "${UNPAIRED}" ]; then
            echo "single" > LAYOUT
            gzip "${UNPAIRED}"
            # Land single-end output in the same read1 slot paired-end uses,
            # so the output block below doesn't need a layout-conditional glob.
            mv "${UNPAIRED}.gz" "${R1}.gz"
        else
            echo "ERROR: fasterq-dump produced neither paired (_1/_2) nor single output for ~{accession}." >&2
            ls -la >&2
            exit 1
        fi

        # ── Best-effort platform lookup ───────────────────────────
        # ENA's filereport mirrors SRA/DDBJ metadata and needs no auth. This
        # is informational only — never fail the task over it, since the
        # fastq extraction above is already done and is the actual point.
        PLATFORM="unknown"
        INSTRUMENT="unknown"
        ENA_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=~{accession}&result=read_run&fields=instrument_platform,instrument_model&format=tsv"
        if command -v curl >/dev/null 2>&1; then
            ENA_TSV="$(curl -fsSL --max-time 30 "${ENA_URL}" 2>/dev/null || true)"
        elif command -v wget >/dev/null 2>&1; then
            ENA_TSV="$(wget -qO- --timeout=30 "${ENA_URL}" 2>/dev/null || true)"
        else
            ENA_TSV=""
        fi
        # Second line is the data row; first is the header.
        DATA_LINE="$(printf '%s\n' "${ENA_TSV}" | sed -n '2p')"
        if [ -n "${DATA_LINE}" ]; then
            PARSED_PLATFORM="$(printf '%s' "${DATA_LINE}" | cut -f2)"
            PARSED_INSTRUMENT="$(printf '%s' "${DATA_LINE}" | cut -f3)"
            [ -n "${PARSED_PLATFORM}" ]   && PLATFORM="${PARSED_PLATFORM}"
            [ -n "${PARSED_INSTRUMENT}" ] && INSTRUMENT="${PARSED_INSTRUMENT}"
        fi
        echo "${PLATFORM}" > PLATFORM
        echo "${INSTRUMENT}" > INSTRUMENT
        echo "${PLATFORM}_$(cat LAYOUT)" > READ_FORMAT
    >>>

    output {
        File   read1             = "~{accession}_1.fastq.gz"
        File?  read2             = "~{accession}_2.fastq.gz"
        File?  orphan_reads      = "~{accession}.fastq.gz"
        String layout            = read_string("LAYOUT")
        String platform          = read_string("PLATFORM")
        String instrument_model  = read_string("INSTRUMENT")
        String read_format       = read_string("READ_FORMAT")
        String sra_tools_version = read_string("SRA_TOOLS_VERSION")
    }

    runtime {
        docker:      "staphb/sratoolkit:3.3.0"
        memory:      memory_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " SSD"
        preemptible: 1
        maxRetries:  2
    }
}
