version 1.0

task align_modbam {

    input {
        File    modbam
        String  sample_name
        File    reference_fasta
        String  minimap2_preset     = "map-ont"
        Int     min_read_length     = 500
        Float   min_mapped_percent  = 0
        Int     cpu                 = 8
        Int     mem_gb              = 32
        Int     disk_gb             = 150
        String  docker              = "nanozoo/minimap2:2.31--c2b4c91@sha256:0b4e1d2204aa10863ae7508f52992665161d14c557e6cfede03e5a0739cf5dc4"
    }

    parameter_meta {
        modbam:             "Modified-basecalled BAM from dorado / ont_basecall_client, carrying MM and ML tags. Aligned or unaligned both work; an aligned input is reduced back to primary records in original orientation before realignment."
        sample_name:        "Some identifier for naming outputs"
        reference_fasta:    "Sequence to align against. The intended use is each isolate's OWN assembly: self-mapping keeps every modification call in the isolate's native coordinates and, critically, means motif discovery reads sequence context from the isolate's real sequence rather than a reference's."
        minimap2_preset:    "minimap2 -x preset (default = map-ont)"
        min_read_length:    "Drop reads shorter than this many bp before alignment, not after. A raw modBAM commonly carries a tail of very short fragments (adapter remnants, truncated translocations) that are mechanically incapable of a confident minimap2 placement; left in, they inflate the read count without ever being mappable and deflate percent_mapped for reasons that have nothing to do with whether the modbam and assembly actually match. 500 sits well clear of that junk (empirically ~50-90bp) while staying below the ~1000bp floor long-read assemblers use — that threshold is tuned for assembly-graph overlap detection, a different job than this task's: every correctly-placed read here is a vote toward modkit's coverage floor, so a shorter-but-real read is still worth keeping (default = 500)"
        min_mapped_percent: "Fail the task if fewer than this percent of primary reads (AFTER the length filter) map. Self-mapping should exceed 95%; a low rate almost always means the modbam and the assembly belong to different isolates. 0 disables the check (default = 0)"
        cpu:                "Number of CPUs delegated to task (default = 8)"
        mem_gb:             "Amount of memory in GB delegated to task (default = 32)"
        disk_gb:            "Amount of disk space in GB delegated to task (default = 150)"
        docker:             "Container image. Needs minimap2 and samtools together; staphb/minimap2 ships minimap2 alone, which is why this points at nanozoo."
    }

    command <<<
        set -euo pipefail

        # Guard the input before spending an alignment on it. MM/ML are what
        # carry the modification calls; a BAM without them produces an empty
        # pileup several tasks later with no error anywhere in between.
        INPUT_MM=$(samtools view ~{modbam} 2>/dev/null | head -n 10000 | grep -c 'MM:Z:' || true)
        if [ "${INPUT_MM}" -eq 0 ]; then
            echo "ERROR: no MM:Z tags in the first 10000 records of ~{modbam}." >&2
            echo "       The run was probably basecalled without modified-base models," >&2
            echo "       or the tags were stripped by an intermediate samtools step." >&2
            exit 1
        fi
        echo "Input carries MM tags (${INPUT_MM} of first 10000 records)"

        # Stage the reference locally so the .fai lands somewhere writable and
        # can be emitted as an output; downstream tasks need it and re-indexing
        # per task is wasted work. Bakta does not read gzipped FASTA either, so
        # decompress here and keep one canonical copy.
        case "~{reference_fasta}" in
            *.gz) gunzip -c ~{reference_fasta} > ref.fa ;;
            *)    cp ~{reference_fasta} ref.fa ;;
        esac
        samtools faidx ref.fa

        # Three things here are load-bearing:
        #
        #   -T MM,ML,MN  samtools fastq drops every auxiliary tag unless asked
        #                for them by name. Omitting this is the single most
        #                common way a methylation pipeline silently yields
        #                nothing. MN is included because modkit uses it to
        #                validate tag length against the read.
        #
        #   -F 0x900     Drop secondary and supplementary records. If the input
        #                was previously aligned, supplementary records carry
        #                hard-clipped SEQ and would corrupt the read set.
        #                Primary records are restored to original orientation
        #                by samtools, which keeps MM/ML — stored in original
        #                read coordinates — consistent.
        #
        #   -y           minimap2 parks the tags in the FASTQ comment field;
        #                -y is what copies that comment back out as SAM tags.
        #                Without it the tags reach minimap2 and die there.
        samtools fastq -@ ~{cpu} -T MM,ML,MN -F 0x900 ~{modbam} > reads.fastq
        N_RAW=$(($(wc -l < reads.fastq) / 4))

        # Length filter, applied here rather than left implicit in the mapping
        # rate: a read too short to place confidently was never going to map
        # regardless of whether the modbam and assembly correspond, so leaving
        # it in the denominator makes percent_mapped measure read-length
        # composition as much as it measures isolate/assembly agreement.
        # samtools fastq has no native length filter, hence the awk pass.
        awk -v minlen=~{min_read_length} '
            NR%4==1 { h=$0 }
            NR%4==2 { s=$0 }
            NR%4==3 { p=$0 }
            NR%4==0 { if (length(s) >= minlen) print h "\n" s "\n" p "\n" $0 }
        ' reads.fastq > filtered.fastq
        N_FILTERED=$(($(wc -l < filtered.fastq) / 4))

        PCT_RETAINED=$(awk -v r="${N_RAW}" -v f="${N_FILTERED}" \
            'BEGIN {if (r>0) printf "%.1f", 100*f/r; else print "NA"}')
        echo "Length filter (>=~{min_read_length} bp): ${N_RAW} -> ${N_FILTERED} reads (${PCT_RETAINED}% retained)"
        echo "${N_RAW}"      > N_READS_RAW
        echo "${N_FILTERED}" > N_READS_FILTERED

        minimap2 -y -ax ~{minimap2_preset} --MD -t ~{cpu} ref.fa filtered.fastq \
            | samtools sort -@ ~{cpu} -o ~{sample_name}.aligned.bam -

        samtools index -@ ~{cpu} ~{sample_name}.aligned.bam

        # The tags survived the input check; confirm they survived the round
        # trip too, rather than discovering it from an empty bedMethyl.
        OUTPUT_MM=$(samtools view ~{sample_name}.aligned.bam | head -n 10000 | grep -c 'MM:Z:' || true)
        if [ "${OUTPUT_MM}" -eq 0 ]; then
            echo "ERROR: MM tags were present in the input but absent after alignment." >&2
            echo "       Check that minimap2 was invoked with -y and samtools fastq with -T." >&2
            exit 1
        fi

        samtools flagstat -@ ~{cpu} ~{sample_name}.aligned.bam > ~{sample_name}_flagstat.txt

        # Mapping rate against a foreign reference is the number that tells you
        # whether this isolate is close enough to the reference to interpret.
        awk '/primary mapped \(/ {gsub(/[(%]/,"",$6); print $6; found=1}
             END {if (!found) print "NA"}' ~{sample_name}_flagstat.txt > MAPPED_PCT

        samtools coverage ~{sample_name}.aligned.bam > ~{sample_name}_coverage.txt

        # Depth-weighted mean across contigs. Per-site modification fractions
        # get noisy fast below ~20x, so this is a gating number, not a nicety.
        awk 'NR>1 {len=$3-$2; total+=$7*len; bases+=len}
             END {if (bases>0) printf "%.2f\n", total/bases; else print "NA"}' \
            ~{sample_name}_coverage.txt > MEAN_DEPTH

        # modbams and assemblies are matched positionally by the caller, which
        # is easy to get wrong and produces no error on its own — just a quietly
        # terrible alignment. A self-mapping run that maps poorly is almost
        # always a mismatched pair.
        MIN_PCT="~{min_mapped_percent}"
        MAPPED="$(cat MAPPED_PCT)"
        if [ "${MIN_PCT}" != "0" ] && [ "${MIN_PCT}" != "0.0" ] && [ "${MAPPED}" != "NA" ]; then
            if awk -v m="${MAPPED}" -v t="${MIN_PCT}" 'BEGIN {exit !(m < t)}'; then
                echo "ERROR: only ${MAPPED}% of primary reads mapped, below the ${MIN_PCT}% floor." >&2
                echo "       For self-mapping this normally means the modbam and the" >&2
                echo "       assembly are from different isolates — check that the input" >&2
                echo "       arrays are in the same order." >&2
                exit 1
            fi
        fi

        minimap2 --version > VERSION
    >>>

    output {
        File    aligned_bam         = "~{sample_name}.aligned.bam"
        File    aligned_bam_index   = "~{sample_name}.aligned.bam.bai"
        File    reference_fai       = "ref.fa.fai"
        File    flagstat            = "~{sample_name}_flagstat.txt"
        File    coverage_txt        = "~{sample_name}_coverage.txt"
        Int     n_reads_raw         = read_int("N_READS_RAW")
        Int     n_reads_length_filtered = read_int("N_READS_FILTERED")
        String  percent_mapped      = read_string("MAPPED_PCT")
        String  mean_depth          = read_string("MEAN_DEPTH")
        String  minimap2_version    = read_string("VERSION")
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
