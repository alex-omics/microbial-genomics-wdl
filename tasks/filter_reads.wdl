version 1.0

# Filter ONT reads on length and mean quality, optionally cropping their ends.
#
# Two filters, because the input can be two things:
#
#   filter_reads_bam    A BAM (typically an unaligned modBAM) is filtered record
#                       for record with samtools, so every tag, MM/ML included,
#                       is carried through untouched. A BAM is never cropped:
#                       MM/ML are positional, and cutting a read's ends without
#                       rewriting them would leave modification calls pointing at
#                       the wrong bases.
#   filter_reads_fastq  FASTQ goes through chopper, which can also crop.
#
# read_stats and summarize_filtering give the same before/after numbers for both.
#
# samtools, chopper and seqkit ship together in docker/read-filter.

task validate_filter_inputs {

    input {
        String  sample_name
        Boolean has_bam
        Boolean has_fastq
        Int?    min_length
        Int?    max_length
        Float?  min_mean_q
        Int?    head_crop
        Int?    tail_crop
        Int     min_reads_out
        Int     min_bases_out
        String  docker = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        sample_name:   "Becomes part of the output filenames"
        has_bam:       "Whether reads_bam was supplied"
        has_fastq:     "Whether reads_fastq was supplied"
        min_length:    "Minimum read length"
        max_length:    "Maximum read length"
        min_mean_q:    "Minimum mean read quality"
        head_crop:     "Bases to remove from the start of each read"
        tail_crop:     "Bases to remove from the end of each read"
        min_reads_out: "Fail if fewer reads than this survive"
        min_bases_out: "Fail if fewer bases than this survive"
        docker:        "Container image"
    }

    meta {
        description: "Fail within seconds if the inputs are contradictory or would do nothing, before any reads are read. Emits the validated sample name; the workflow's other tasks take it from here, which is what makes them wait for this check."
    }

    Int min_len  = select_first([min_length, 0])
    Int max_len  = select_first([max_length, 0])
    Int head     = select_first([head_crop, 0])
    Int tail     = select_first([tail_crop, 0])
    Float min_q  = select_first([min_mean_q, 0.0])
    Boolean any_filter = defined(min_length) || defined(max_length) || defined(min_mean_q)

    command <<<
        set -euo pipefail

        # The name goes into filenames and shell commands
        if ! printf '%s' "~{sample_name}" | grep -qE '^[A-Za-z0-9._-]+$'; then
            echo "ERROR: sample_name '~{sample_name}' must be non-empty and use only letters," >&2
            echo "       digits, dot, underscore and hyphen." >&2; exit 1
        fi

        if [ "~{has_bam}" = "true" ] && [ "~{has_fastq}" = "true" ]; then
            echo "ERROR: supply reads_bam or reads_fastq, not both." >&2; exit 1
        fi
        if [ "~{has_bam}" = "false" ] && [ "~{has_fastq}" = "false" ]; then
            echo "ERROR: supply reads_bam or reads_fastq." >&2; exit 1
        fi

        # An unset filter would copy the reads and report 100% retained, which
        # reads as a successful filtering run.
        if [ "~{any_filter}" = "false" ] && [ ~{head} -eq 0 ] && [ ~{tail} -eq 0 ]; then
            echo "ERROR: nothing to do. Set at least one of min_length, max_length," >&2
            echo "       min_mean_q, head_crop or tail_crop." >&2; exit 1
        fi

        if [ "~{has_bam}" = "true" ] && { [ ~{head} -gt 0 ] || [ ~{tail} -gt 0 ]; }; then
            echo "ERROR: head_crop and tail_crop are not supported for a BAM. A modBAM's" >&2
            echo "       MM/ML tags are positional, so cropping a read would leave its" >&2
            echo "       modification calls on the wrong bases. Filter the BAM, or crop" >&2
            echo "       FASTQ reads (which carry no modification calls)." >&2; exit 1
        fi

        if [ "~{defined(min_length)}" = "true" ] && [ "~{defined(max_length)}" = "true" ] && [ ~{min_len} -gt ~{max_len} ]; then
            echo "ERROR: min_length (~{min_len}) is greater than max_length (~{max_len})." >&2; exit 1
        fi
        if [ ~{min_len} -lt 0 ] || [ ~{max_len} -lt 0 ] || [ ~{head} -lt 0 ] || [ ~{tail} -lt 0 ]; then
            echo "ERROR: lengths and crops cannot be negative." >&2; exit 1
        fi
        if awk -v q="~{min_q}" 'BEGIN {exit !(q < 0)}'; then
            echo "ERROR: min_mean_q cannot be negative." >&2; exit 1
        fi
        if [ ~{min_reads_out} -lt 0 ] || [ ~{min_bases_out} -lt 0 ]; then
            echo "ERROR: min_reads_out and min_bases_out cannot be negative." >&2; exit 1
        fi
        echo "Inputs validated"
    >>>

    output {
        String validated_name = sample_name
    }

    runtime {
        docker:         docker
        memory:         "1 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     0
    }
}


task filter_reads_bam {

    input {
        File    reads_bam
        String  sample_name
        Int?    min_length
        Int?    max_length
        Float?  min_mean_q
        String  qscore_source   = "auto"
        Int     cpu             = 8
        Int     mem_gb          = 16
        Int?    disk_gb
        Int     preemptible     = 1
        String  docker          = "aarvani1/read-filter:1.0.0@sha256:a59effac6d94b62d9fe1de7c250d972f888c8df8d011baf31a4b8e4b6f40b11b"
    }

    parameter_meta {
        reads_bam:     "Unaligned BAM, typically a modBAM. Every record is kept or dropped whole, with all of its tags."
        sample_name:   "Names the output"
        min_length:    "Drop reads shorter than this (bp). Inclusive: a read of exactly this length is kept"
        max_length:    "Drop reads longer than this (bp). Inclusive"
        min_mean_q:    "Drop reads whose mean quality is below this. Ties are decided by rounding, so do not rely on a read landing exactly on the value"
        qscore_source: "auto, qs_tag or computed. qs_tag reads the mean quality the basecaller wrote to each record's qs tag (dorado does); computed derives it from the quality string with chopper. auto uses qs_tag if every record has one, otherwise computed (default = auto)"
        cpu:           "Number of CPUs delegated to task (default = 8)"
        mem_gb:        "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:       "Disk in GB; default is 3x the input plus 20"
        preemptible:   "Preemptible attempts (default = 1)"
        docker:        "Container image with samtools and chopper"
    }

    meta {
        description: "Filter a BAM on read length and mean quality, record for record, so MM/ML modification tags are carried through. Refuses an aligned BAM, because the sequence of a clipped or supplementary record is not the read."
    }

    Boolean has_min_len = defined(min_length)
    Boolean has_max_len = defined(max_length)
    Boolean has_min_q   = defined(min_mean_q)
    Int     min_len     = select_first([min_length, 0])
    Int     max_len     = select_first([max_length, 0])
    Float   min_q       = select_first([min_mean_q, 0.0])
    Int     disk        = select_first([disk_gb, ceil(size(reads_bam, "GB") * 3) + 20])

    command <<<
        set -euo pipefail

        case "~{qscore_source}" in
            auto|qs_tag|computed) ;;
            *) echo "ERROR: qscore_source must be auto, qs_tag or computed, not '~{qscore_source}'." >&2; exit 1 ;;
        esac

        # Length is read from SEQ, which is the read only if the BAM is
        # unaligned. A mapped record can be hard-clipped or supplementary, and
        # its SEQ is then a fragment.
        N_SQ=$(samtools view -H ~{reads_bam} | grep -c '^@SQ' || true)
        if [ "${N_SQ}" -gt 0 ]; then
            echo "ERROR: ~{reads_bam} has ${N_SQ} @SQ header lines, so it is aligned." >&2
            echo "       Filter the reads before aligning them, or filter the FASTQ." >&2
            exit 1
        fi

        # Did the input carry modification tags? Checked again on the output.
        IN_MM=$(samtools view ~{reads_bam} 2>/dev/null | head -n 10000 | grep -c 'MM:Z:' || true)

        # Length-only filtering needs no quality at all.
        SOURCE="none"
        if [ "~{has_min_q}" = "true" ]; then
            SOURCE="~{qscore_source}"
            if [ "${SOURCE}" != "computed" ]; then
                # A record without qs would fail "[qs]>=Q" and be dropped without
                # a word, so count them over the whole file rather than a sample.
                N_NO_QS=$(samtools view -@ ~{cpu} -c -e '!exists([qs])' ~{reads_bam})
                if [ "${N_NO_QS}" -eq 0 ]; then
                    SOURCE="qs_tag"
                elif [ "${SOURCE}" = "qs_tag" ]; then
                    echo "ERROR: ${N_NO_QS} records in ~{reads_bam} have no qs tag." >&2
                    echo "       Use qscore_source=auto or computed to derive quality from the quality string." >&2
                    exit 1
                else
                    SOURCE="computed"
                fi
            fi
        fi
        echo "Quality source: ${SOURCE}"
        echo "${SOURCE}" > QSCORE_SOURCE

        if [ "${SOURCE}" = "computed" ]; then
            # Read quality is computed from the quality string by chopper, the
            # same arithmetic the FASTQ filter uses. chopper cannot read a BAM
            # and cannot write one without losing the tags, so it decides which
            # reads pass and samtools extracts exactly those records.
            CHOPPER_ARGS=(-t ~{cpu} -q ~{min_q})
            if [ "~{has_min_len}" = "true" ]; then CHOPPER_ARGS+=(--minlength ~{min_len}); fi
            if [ "~{has_max_len}" = "true" ]; then CHOPPER_ARGS+=(--maxlength ~{max_len}); fi

            samtools fastq -@ ~{cpu} -n ~{reads_bam} 2>/dev/null \
                | chopper "${CHOPPER_ARGS[@]}" 2>chopper.log \
                | awk 'NR%4==1 {sub(/^@/, "", $1); print $1}' > keep.txt
            cat chopper.log >&2

            # An empty -N list selects nothing (not everything); the output is
            # then a header-only BAM and the summary step reports it.
            samtools view -@ ~{cpu} -b -N keep.txt -o ~{sample_name}.filtered.bam ~{reads_bam}
        else
            EXPR=""
            add() { EXPR="${EXPR:+${EXPR} && }$1"; }
            if [ "~{has_min_len}" = "true" ]; then add "length(seq)>=~{min_len}"; fi
            if [ "~{has_max_len}" = "true" ]; then add "length(seq)<=~{max_len}"; fi
            if [ "~{has_min_q}"   = "true" ]; then add "[qs]>=~{min_q}"; fi
            echo "samtools filter expression: ${EXPR}"
            samtools view -@ ~{cpu} -b -e "${EXPR}" -o ~{sample_name}.filtered.bam ~{reads_bam}
        fi

        if [ "${IN_MM}" -gt 0 ]; then
            N_OUT=$(samtools view -c ~{sample_name}.filtered.bam)
            if [ "${N_OUT}" -gt 0 ]; then
                OUT_MM=$(samtools view ~{sample_name}.filtered.bam 2>/dev/null | head -n 10000 | grep -c 'MM:Z:' || true)
                if [ "${OUT_MM}" -eq 0 ]; then
                    echo "ERROR: the input carried MM tags but the filtered BAM does not." >&2
                    exit 1
                fi
            fi
        fi
    >>>

    output {
        File   filtered_bam   = "~{sample_name}.filtered.bam"
        String qscore_source_used = read_string("QSCORE_SOURCE")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk} SSD"
        preemptible:    preemptible
        maxRetries:     1
    }
}


task filter_reads_fastq {

    input {
        File    reads_fastq
        String  sample_name
        Int?    min_length
        Int?    max_length
        Float?  min_mean_q
        Int?    head_crop
        Int?    tail_crop
        Int     cpu             = 8
        Int     mem_gb          = 16
        Int?    disk_gb
        Int     preemptible     = 1
        String  docker          = "aarvani1/read-filter:1.0.0@sha256:a59effac6d94b62d9fe1de7c250d972f888c8df8d011baf31a4b8e4b6f40b11b"
    }

    parameter_meta {
        reads_fastq: "FASTQ, plain or gzipped"
        sample_name: "Names the output"
        min_length:  "Drop reads shorter than this (bp), measured after cropping. Inclusive"
        max_length:  "Drop reads longer than this (bp), measured after cropping. Inclusive"
        min_mean_q:  "Drop reads whose mean quality is below this. Mean quality is taken over error probabilities, not over the scores. It is judged on the whole read before cropping, so cropping a bad end off a read does not rescue it. Ties are decided by rounding"
        head_crop:   "Bases removed from the start of every read"
        tail_crop:   "Bases removed from the end of every read"
        cpu:         "Number of CPUs delegated to task (default = 8)"
        mem_gb:      "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:     "Disk in GB; default is 3x the input plus 20"
        preemptible: "Preemptible attempts (default = 1)"
        docker:      "Container image with chopper and pigz"
    }

    meta {
        description: "Filter FASTQ on read length and mean quality with chopper, optionally cropping read ends first. Output is gzipped."
    }

    Boolean has_min_len = defined(min_length)
    Boolean has_max_len = defined(max_length)
    Int     min_len     = select_first([min_length, 0])
    Int     max_len     = select_first([max_length, 0])
    Float   min_q       = select_first([min_mean_q, 0.0])
    Int     head        = select_first([head_crop, 0])
    Int     tail        = select_first([tail_crop, 0])
    Int     disk        = select_first([disk_gb, ceil(size(reads_fastq, "GB") * 3) + 20])

    command <<<
        set -euo pipefail

        CHOPPER_ARGS=(-t ~{cpu} -q ~{min_q})
        if [ "~{has_min_len}" = "true" ]; then CHOPPER_ARGS+=(--minlength ~{min_len}); fi
        if [ "~{has_max_len}" = "true" ]; then CHOPPER_ARGS+=(--maxlength ~{max_len}); fi

        # chopper ignores --headcrop and --tailcrop (with a warning, not an
        # error) unless told how to trim, so a crop that silently did nothing
        # would otherwise pass for a crop that worked.
        if [ ~{head} -gt 0 ] || [ ~{tail} -gt 0 ]; then
            CHOPPER_ARGS+=(--trim-approach fixed-crop --headcrop ~{head} --tailcrop ~{tail})
        fi
        echo "chopper ${CHOPPER_ARGS[*]}"

        # gzip -dcf passes plain FASTQ through unchanged
        gzip -dcf ~{reads_fastq} \
            | chopper "${CHOPPER_ARGS[@]}" 2>chopper.log \
            | pigz -p ~{cpu} > ~{sample_name}.filtered.fastq.gz
        cat chopper.log >&2

        if grep -qi "ignored" chopper.log; then
            echo "ERROR: chopper ignored an option (see above)." >&2; exit 1
        fi
    >>>

    output {
        File filtered_fastq = "~{sample_name}.filtered.fastq.gz"
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk} SSD"
        preemptible:    preemptible
        maxRetries:     1
    }
}


task read_stats {

    input {
        File    reads
        Boolean is_bam
        String  sample_name
        String  stage
        Int     cpu         = 4
        Int     mem_gb      = 8
        Int?    disk_gb
        Int     preemptible = 1
        String  docker      = "aarvani1/read-filter:1.0.0@sha256:a59effac6d94b62d9fe1de7c250d972f888c8df8d011baf31a4b8e4b6f40b11b"
    }

    parameter_meta {
        reads:       "BAM or FASTQ (plain or gzipped)"
        is_bam:      "Whether reads is a BAM"
        sample_name: "Value for the sample column of the stats table"
        stage:       "Value for the stage column of the stats table (input or filtered)"
        cpu:         "Number of CPUs delegated to task (default = 4)"
        mem_gb:      "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:     "Disk in GB; default is the input size plus 10"
        preemptible: "Preemptible attempts (default = 1)"
        docker:      "Container image with samtools and seqkit"
    }

    meta {
        description: "Read count, bases, length distribution and quality of a read set, with seqkit. A BAM is streamed through samtools fastq, so nothing is written to disk."
    }

    Int disk = select_first([disk_gb, ceil(size(reads, "GB")) + 10])

    command <<<
        set -euo pipefail

        if [ "~{is_bam}" = "true" ]; then
            samtools fastq -@ ~{cpu} -n ~{reads} 2>/dev/null | seqkit stats -a -T -j ~{cpu} - > raw_stats.tsv
        else
            gzip -dcf ~{reads} | seqkit stats -a -T -j ~{cpu} - > raw_stats.tsv
        fi

        # seqkit's first column is the file name, which here is "-" (stdin).
        # Swap it for sample and stage columns so the tables can be stacked.
        awk -F'\t' -v OFS='\t' -v s="~{sample_name}" -v g="~{stage}" \
            'NR==1 {$1="sample\tstage"} NR>1 {$1=s "\t" g} {print}' raw_stats.tsv > stats.tsv

        # num_seqs and sum_len, by column name rather than position
        awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) col[$i]=i; next}
                    {print $col["num_seqs"] > "N_READS"; print $col["sum_len"] > "N_BASES";
                     print $col["N50"] > "N50"}' stats.tsv
    >>>

    output {
        File stats_tsv = "stats.tsv"
        Int  n_reads   = read_int("N_READS")
        Int  n_bases   = read_int("N_BASES")
        Int  n50       = read_int("N50")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk} SSD"
        preemptible:    preemptible
        maxRetries:     1
    }
}


task summarize_filtering {

    input {
        File    stats_before
        File    stats_after
        String  sample_name
        Int     reads_in
        Int     reads_out
        Int     bases_in
        Int     bases_out
        Int     min_reads_out
        Int     min_bases_out
        String  docker = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        stats_before:  "read_stats table of the input"
        stats_after:   "read_stats table of the filtered reads"
        sample_name:   "Names the output"
        reads_in:      "Reads before filtering"
        reads_out:     "Reads after filtering"
        bases_in:      "Bases before filtering"
        bases_out:     "Bases after filtering"
        min_reads_out: "Fail if fewer reads than this survive"
        min_bases_out: "Fail if fewer bases than this survive"
        docker:        "Container image"
    }

    meta {
        description: "Join the before and after tables, work out what was retained, and fail if the filter left too little. An over-tight threshold otherwise ends as an empty or near-empty file that the next workflow trips over."
    }

    command <<<
        set -euo pipefail

        # One table, input row first. The after table's header is dropped.
        { cat ~{stats_before}; tail -n +2 ~{stats_after}; } > ~{sample_name}.filter_stats.tsv

        awk -v a=~{reads_in} -v b=~{reads_out} 'BEGIN {printf "%.2f\n", (a>0 ? 100*b/a : 0)}' > PCT_READS
        awk -v a=~{bases_in} -v b=~{bases_out} 'BEGIN {printf "%.2f\n", (a>0 ? 100*b/a : 0)}' > PCT_BASES
        echo "~{sample_name}: ~{reads_in} -> ~{reads_out} reads ($(cat PCT_READS)% retained), ~{bases_in} -> ~{bases_out} bases ($(cat PCT_BASES)% retained)"

        if [ ~{reads_in} -eq 0 ]; then
            echo "ERROR: the input has no reads." >&2; exit 1
        fi
        if [ ~{reads_out} -lt ~{min_reads_out} ]; then
            echo "ERROR: ~{reads_out} reads survived filtering, fewer than min_reads_out (~{min_reads_out})." >&2
            echo "       The thresholds may be too tight for this run; compare ~{sample_name}.filter_stats.tsv." >&2
            exit 1
        fi
        if [ ~{bases_out} -lt ~{min_bases_out} ]; then
            echo "ERROR: ~{bases_out} bases survived filtering, fewer than min_bases_out (~{min_bases_out})." >&2
            exit 1
        fi
    >>>

    output {
        File  filter_stats_tsv       = "~{sample_name}.filter_stats.tsv"
        Float pct_reads_retained = read_float("PCT_READS")
        Float pct_bases_retained = read_float("PCT_BASES")
    }

    runtime {
        docker:         docker
        memory:         "1 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     0
    }
}
