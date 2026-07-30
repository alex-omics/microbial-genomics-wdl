version 1.0

task align_modbam {

    input {
        File    modbam
        String  sample_name
        File    reference_fasta
        String  minimap2_preset = "map-ont"
        Int     cpu             = 8
        Int     mem_gb          = 32
        Int     disk_gb         = 150
        String  docker          = "nanozoo/minimap2:2.31--c2b4c91"
    }

    parameter_meta {
        modbam:          "Modified-basecalled BAM from dorado / ont_basecall_client, carrying MM and ML tags. Aligned or unaligned both work; an aligned input is reduced back to primary records in original orientation before realignment."
        sample_name:     "Some identifier for naming outputs"
        reference_fasta: "Reference to align against. For comparative methylomics across isolates this must be the SAME reference for every sample, or the resulting bedMethyl files share no coordinate space."
        minimap2_preset: "minimap2 -x preset (default = map-ont)"
        cpu:             "Number of CPUs delegated to task (default = 8)"
        mem_gb:          "Amount of memory in GB delegated to task (default = 32)"
        disk_gb:         "Amount of disk space in GB delegated to task (default = 150)"
        docker:          "Container image. Needs minimap2 and samtools together; staphb/minimap2 ships minimap2 alone, which is why this points at nanozoo."
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

        samtools faidx ~{reference_fasta} 2>/dev/null || true

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
        samtools fastq -@ ~{cpu} -T MM,ML,MN -F 0x900 ~{modbam} \
            | minimap2 -y -ax ~{minimap2_preset} --MD -t ~{cpu} ~{reference_fasta} - \
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

        minimap2 --version > VERSION
    >>>

    output {
        File    aligned_bam         = "~{sample_name}.aligned.bam"
        File    aligned_bam_index   = "~{sample_name}.aligned.bam.bai"
        File    flagstat            = "~{sample_name}_flagstat.txt"
        File    coverage_txt        = "~{sample_name}_coverage.txt"
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
