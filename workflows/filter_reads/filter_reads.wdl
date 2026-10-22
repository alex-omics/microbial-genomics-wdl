version 1.0

import "../../tasks/filter_reads.wdl" as filter_tasks

workflow filter_reads {

    meta {
        description: "Filter one sample's ONT reads on length and mean quality. Takes an unaligned BAM (including a modBAM, whose MM/ML modification tags are carried through untouched) or FASTQ, and returns reads in the same format with a before/after table of read counts, bases, N50 and quality. FASTQ reads can also be cropped at either end. Filters on the reads' own properties only, with no taxonomic or reference-based filtering. One sample per run."
        author: "Alex Arvanitis"
    }

    input {
        String  sample_name

        # Reads: a BAM or a FASTQ
        File?   reads_bam
        File?   reads_fastq

        # Filters. Unset means no limit; at least one filter or crop must be set.
        Int?    min_length
        Int?    max_length
        Float?  min_mean_q

        # FASTQ only
        Int?    head_crop
        Int?    tail_crop

        String  qscore_source   = "auto"

        # Fail rather than hand an empty or near-empty file to the next step
        Int     min_reads_out   = 1
        Int     min_bases_out   = 0

        Int     cpu             = 8
        Int     mem_gb          = 16
        Int?    disk_gb
        Int     preemptible     = 1
        String  docker          = "aarvani1/read-filter:1.0.0@sha256:a59effac6d94b62d9fe1de7c250d972f888c8df8d011baf31a4b8e4b6f40b11b"
    }

    parameter_meta {
        sample_name:   "Label for the outputs"
        reads_bam:     "Unaligned BAM, such as a dorado modBAM. Use this or reads_fastq"
        reads_fastq:   "FASTQ, plain or gzipped. Use this or reads_bam"
        min_length:    "Keep reads of at least this many bp (inclusive)"
        max_length:    "Keep reads of at most this many bp (inclusive)"
        min_mean_q:    "Keep reads with mean quality of at least this. The right floor depends on the chemistry and basecalling model; look at the run's own quality distribution before choosing, and see the README"
        head_crop:     "FASTQ only. Bases removed from the start of each read, before the length filter"
        tail_crop:     "FASTQ only. Bases removed from the end of each read, before the length filter"
        qscore_source: "BAM only. auto, qs_tag or computed; see the README (default = auto)"
        min_reads_out: "Fail if fewer reads than this survive (default = 1, which fails on an empty result)"
        min_bases_out: "Fail if fewer bases than this survive (default = 0)"
        disk_gb:       "Disk for the filtering step; default is 3x the input plus 20"
        docker:        "Image with samtools, chopper and seqkit"
    }

    Boolean is_bam = defined(reads_bam)

    call filter_tasks.validate_filter_inputs {
        input:
            sample_name   = sample_name,
            has_bam       = defined(reads_bam),
            has_fastq     = defined(reads_fastq),
            min_length    = min_length,
            max_length    = max_length,
            min_mean_q    = min_mean_q,
            head_crop     = head_crop,
            tail_crop     = tail_crop,
            min_reads_out = min_reads_out,
            min_bases_out = min_bases_out
    }

    # Everything below takes the name from the validation step, so none of it starts until the inputs pass
    String name = validate_filter_inputs.validated_name

    # The statistics of the input run alongside the filter, not after it
    call filter_tasks.read_stats as stats_before {
        input:
            reads       = select_first([reads_bam, reads_fastq]),
            is_bam      = is_bam,
            sample_name = name,
            stage       = "input",
            preemptible = preemptible,
            docker      = docker
    }

    if (is_bam) {
        call filter_tasks.filter_reads_bam {
            input:
                reads_bam     = select_first([reads_bam]),
                sample_name   = name,
                min_length    = min_length,
                max_length    = max_length,
                min_mean_q    = min_mean_q,
                qscore_source = qscore_source,
                cpu           = cpu,
                mem_gb        = mem_gb,
                disk_gb       = disk_gb,
                preemptible   = preemptible,
                docker        = docker
        }
    }

    if (!is_bam) {
        call filter_tasks.filter_reads_fastq {
            input:
                reads_fastq = select_first([reads_fastq]),
                sample_name = name,
                min_length  = min_length,
                max_length  = max_length,
                min_mean_q  = min_mean_q,
                head_crop   = head_crop,
                tail_crop   = tail_crop,
                cpu         = cpu,
                mem_gb      = mem_gb,
                disk_gb     = disk_gb,
                preemptible = preemptible,
                docker      = docker
        }
    }

    File filtered = select_first([filter_reads_bam.filtered_bam, filter_reads_fastq.filtered_fastq])

    call filter_tasks.read_stats as stats_after {
        input:
            reads       = filtered,
            is_bam      = is_bam,
            sample_name = name,
            stage       = "filtered",
            preemptible = preemptible,
            docker      = docker
    }

    call filter_tasks.summarize_filtering {
        input:
            stats_before  = stats_before.stats_tsv,
            stats_after   = stats_after.stats_tsv,
            sample_name   = name,
            reads_in      = stats_before.n_reads,
            reads_out     = stats_after.n_reads,
            bases_in      = stats_before.n_bases,
            bases_out     = stats_after.n_bases,
            min_reads_out = min_reads_out,
            min_bases_out = min_bases_out
    }

    output {
        File    filtered_reads            = filtered
        File    filter_stats_tsv          = summarize_filtering.filter_stats_tsv
        Int     filter_reads_in           = stats_before.n_reads
        Int     filter_reads_out          = stats_after.n_reads
        Int     filter_bases_in           = stats_before.n_bases
        Int     filter_bases_out          = stats_after.n_bases
        Float   filter_pct_reads_retained = summarize_filtering.pct_reads_retained
        Float   filter_pct_bases_retained = summarize_filtering.pct_bases_retained
        Int     filter_n50_in             = stats_before.n50
        Int     filter_n50_out            = stats_after.n50
        String  filter_qscore_source      = select_first([filter_reads_bam.qscore_source_used, if defined(min_mean_q) then "computed" else "none"])
    }
}
