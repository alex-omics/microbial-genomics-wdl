version 1.0

import "../../tasks/listt.wdl" as listt

# workflow call_ospA_type
#
# Assigns the OspA in silico type (IST) of one Borrelia burgdorferi s.l.
# isolate with LISTT (Lee et al. 2024, Microb Genom 10:001252), from its
# assembly and, optionally, from its Illumina paired-end reads. When both are
# supplied they are typed independently and reconciled into one call and one
# flag.
#
# One sample per execution, with no Array inputs. Terra scatters across data
# table rows, so a thousand isolates are a thousand independent executions that
# succeed and fail separately. Point assembly_fasta at this.assembly_fasta (and
# read1/read2 at this.read1/this.read2) and select every row.
#
# LISTT's code is unmodified. The identity thresholds that the paper describes
# but the released code never applies are enforced here, and ISTs above
# max_trusted_ist (default 17: those in the paper) are reported as provisional
# rather than as calls. See the README for the flags and for what is and is not
# validated.

workflow call_ospA_type {

    meta {
        description: "Types the OspA in silico type (IST) of one Borrelia burgdorferi s.l. isolate with LISTT from its assembly and, optionally, its Illumina paired-end reads, and reports one call with a confidence flag."
        author: "Alex Arvanitis"
    }

    input {
        String  sample_name
        File    assembly_fasta
        File?   read1
        File?   read2
        Int     max_trusted_ist = 17

        Int     cpu         = 2
        Int     mem_gb      = 4
        Int?    disk_gb
        Int     preemptible = 1
        String  docker      = "aarvani1/listt:1.2.0@sha256:c016417e25c3e4b96d4fc98bedbf4d77e15608674e5e5320cdaff1df6a1a38b2"
    }

    parameter_meta {
        sample_name:     "Sample label, used in the outputs - wire to this.<table>_id"
        assembly_fasta:  "Genome assembly FASTA (short- or long-read; plain or gzip). Required: it is the one input every isolate has"
        read1:           "Optional Illumina paired-end R1 FASTQ. LISTT's read pileup is built for Illumina paired-end data, which also makes it an independent check on the assembly. Give read1 and read2 together or neither"
        read2:           "Optional Illumina paired-end R2 FASTQ"
        max_trusted_ist: "Highest IST number reported as a call (default = 17, the ISTs described in the paper). A closest match to a higher IST, or to a variant that has no IST, is reported as 'Unknown' (flag UK) with the closest IST in ospA_details"
        cpu:             "CPUs per LISTT run (default = 2)"
        mem_gb:          "Memory in GB per LISTT run (default = 4)"
        disk_gb:         "Disk in GB per LISTT run. Default is 10 plus twice the input size"
        preemptible:     "Preemptible attempts per LISTT run (default = 1)"
        docker:          "Container image"
    }

    call listt.validate_ospA_inputs {
        input:
            sample_name     = sample_name,
            has_read1       = defined(read1),
            has_read2       = defined(read2),
            max_trusted_ist = max_trusted_ist,
            docker          = docker
    }

    call listt.bbsl_ospA_listt as from_assembly {
        input:
            sample_name     = sample_name,
            mode            = "assembly",
            assembly_fasta  = assembly_fasta,
            max_trusted_ist = max_trusted_ist,
            validated       = validate_ospA_inputs.validated,
            cpu             = cpu,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible,
            docker          = docker
    }

    if (defined(read1) && defined(read2)) {
        call listt.bbsl_ospA_listt as from_reads {
            input:
                sample_name     = sample_name,
                mode            = "reads",
                read1           = read1,
                read2           = read2,
                max_trusted_ist = max_trusted_ist,
                validated       = validate_ospA_inputs.validated,
                cpu             = cpu,
                mem_gb          = mem_gb,
                disk_gb         = disk_gb,
                preemptible     = preemptible,
                docker          = docker
        }
    }

    call listt.reconcile_ospA_calls {
        input:
            sample_name      = sample_name,
            assembly_tsv     = from_assembly.mode_tsv,
            reads_tsv        = from_reads.mode_tsv,
            assembly_protein = from_assembly.protein_fasta,
            reads_protein    = from_reads.protein_fasta,
            assembly_raw     = from_assembly.raw_tar,
            reads_raw        = from_reads.raw_tar,
            listt_version    = from_assembly.listt_version,
            max_trusted_ist  = max_trusted_ist,
            docker           = docker
    }

    output {
        String ospA_call          = reconcile_ospA_calls.final_call           # -> this.ospA_call
        String ospA_flag          = reconcile_ospA_calls.final_flag           # -> this.ospA_flag
        String ospA_species       = reconcile_ospA_calls.final_species        # -> this.ospA_species
        String ospA_identity      = reconcile_ospA_calls.final_identity       # -> this.ospA_identity
        String ospA_details       = reconcile_ospA_calls.final_details        # -> this.ospA_details
        File   ospA_results_tsv   = reconcile_ospA_calls.results_tsv    # -> this.ospA_results_tsv
        File   ospA_protein_fasta = reconcile_ospA_calls.protein_fasta  # -> this.ospA_protein_fasta
        File   ospA_raw_tar       = reconcile_ospA_calls.raw_tar        # -> this.ospA_raw_tar
        String ospA_listt_version = from_assembly.listt_version         # -> this.ospA_listt_version
    }
}
