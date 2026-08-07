version 1.0

import "../../tasks/sra_tools.wdl" as sra_tools

# workflow fetch_reads_from_sra
#
# Fetches a single SRA/ENA/DDBJ run accession and emits gzip-compressed
# FASTQ, split into R1/R2 for paired-end runs or a single file for
# single-end runs (e.g. ONT/PacBio long-read data). Designed for Terra data
# table use: point "accession" at this.SRA_ID, select all rows, and Terra
# scatters across rows automatically — outputs (including read1, read2, and
# sra_tools_version) write straight back to the table, as validated in
# production against a live 300+-row data table.
#
# Deliberately kept to one accession per call rather than an Array input:
# the per-row Terra scatter is the whole point of this shape, and an
# Array[String] version would trade that away for a batch-run pattern this
# workflow doesn't need. Similarly, this only ever fetches a single run
# (SRR/ERR/DRR) — it does not resolve SRX/SRS/SRP/PRJ* container accessions
# to their constituent runs. See fasterq_dump's parameter_meta and the
# README for both.

workflow fetch_reads_from_sra {

    meta {
        description: "Fetches a single SRA/ENA/DDBJ run accession and emits gzip-compressed FASTQ - split R1/R2 for paired-end, a single file for single-end - along with the platform and layout actually extracted."
        author: "Alex Arvanitis"
    }

    input {
        String accession
        Int    disk_gb   = 50
        Int    cpu       = 4
        Int    memory_gb = 8
    }

    parameter_meta {
        accession: "Single-run SRA/ENA/DDBJ accession (SRR, ERR, or DRR) - wire to this.SRA_ID"
        disk_gb:   "Disk space in GB (default=50; fine for MiSeq-scale Illumina, bump substantially for a full ONT/PacBio flowcell)"
        cpu:       "CPUs for fasterq-dump (default=4)"
        memory_gb: "Memory in GB (default=8)"
    }

    call sra_tools.fasterq_dump {
        input:
            accession = accession,
            disk_gb   = disk_gb,
            cpu       = cpu,
            memory_gb = memory_gb
    }

    output {
        File   read1             = fasterq_dump.read1              # -> this.read1
        File?  read2             = fasterq_dump.read2              # -> this.read2
        File?  orphan_reads      = fasterq_dump.orphan_reads        # -> this.orphan_reads
        String layout            = fasterq_dump.layout              # -> this.layout
        String platform          = fasterq_dump.platform            # -> this.platform
        String instrument_model  = fasterq_dump.instrument_model    # -> this.instrument_model
        String sra_tools_version = fasterq_dump.sra_tools_version   # -> this.sra_tools_version
    }
}
