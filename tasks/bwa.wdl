version 1.0

# bwa 0.7.17 + samtools 1.16.1 in one image. staphb/bwa ships bwa only, and
# bwa_mem pipes straight into samtools sort, so it cannot be used here.

task bwa_index {

    input {
        File    reference_fasta
        String  reference_name
        File?   annotation
        Int     cpu        = 2
        Int     mem_gb     = 8
        Int     disk_gb    = 20
        String  docker     = "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40@sha256:9548dc56bdc0b734cd3767f9eea0f9d0ea1c44b35ef5fb35b0f746807cacbeea"
    }

    parameter_meta {
        reference_fasta: "Reference or assembly FASTA (optionally gzipped)"
        reference_name:  "Identifier used to name the output bundle"
        annotation:      "Optional GFF3 for this reference. Only used to check its sequence names against the FASTA before any alignment is billed."
        cpu:             "Number of CPUs delegated to task (default = 2)"
        mem_gb:          "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:         "Amount of disk space in GB delegated to task (default = 20)"
        docker:          "Container image"
    }

    meta {
        description: "Build a BWA index for one reference and return it as a single bundle (FASTA + index files), so a per-sample alignment can never be given an index that belongs to a different genome."
    }

    command <<<
        set -euo pipefail
        # `bwa` with no arguments prints its usage and exits 1, which pipefail would treat as fatal.
        (bwa 2>&1 || true) | grep -i '^Version' | sed 's/^/BWA /' | tee BWA_VERSION

        # Read a file that may or may not be gzipped. `zcat -f` is not portable
        # across the images used here (some zcats reject plain text outright).
        plain() {
            if [ "$(head -c2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then gzip -dc "$1"; else cat "$1"; fi
        }

        mkdir ref
        plain ~{reference_fasta} > ref/reference.fasta

        # With a per-isolate assembly the GFF and FASTA must agree on contig
        # names. If they do not, featureCounts reports zero for everything and
        # nothing errors, so check here, before alignment.
        if [ -n "~{select_first([annotation, ''])}" ]; then
            grep '^>' ref/reference.fasta | sed 's/^>//; s/[[:space:]].*//' | sort -u > fasta_contigs.txt
            plain ~{select_first([annotation, ''])} \
                | awk '/^##FASTA/{skip=1} skip||/^#/{next} NF{print $1}' | sort -u > gff_contigs.txt
            # (comm is not in this image; awk set-difference instead.)
            MISSING="$(awk 'NR==FNR{seen[$0]; next} !($0 in seen)' fasta_contigs.txt gff_contigs.txt)"
            if [ -n "${MISSING}" ]; then
                echo "ERROR: annotation for ~{reference_name} refers to sequences that are not in its FASTA:" >&2
                echo "${MISSING}" | head -n 10 >&2
                echo "The GFF must have been made from this exact assembly with contig names kept" >&2
                echo "(Bakta needs --keep-contig-headers)." >&2
                exit 1
            fi
        fi

        bwa index ref/reference.fasta
        # BusyBox tar in this image has no -z, so compress through a pipe.
        tar -cf - -C ref . | gzip > "~{reference_name}.bwa_ref.tar.gz"
    >>>

    output {
        File   reference_bundle = "~{reference_name}.bwa_ref.tar.gz"
        String bwa_version      = read_string("BWA_VERSION")
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


task bwa_mem {

    input {
        File     read1
        File     read2
        String   sample_name
        File     reference_bundle
        Int      mapq_min              = 20
        Boolean  mark_secondary        = true
        String   samtools_filter_flags = "0x904"
        String   extra_args            = ""
        Int      cpu                   = 8
        Int      mem_gb                = 16
        Int      disk_gb               = 100
        String   docker                = "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40@sha256:9548dc56bdc0b734cd3767f9eea0f9d0ea1c44b35ef5fb35b0f746807cacbeea"
    }

    parameter_meta {
        read1:                 "Forward reads (FASTQ, optionally gzipped)"
        read2:                 "Reverse reads (FASTQ, optionally gzipped)"
        sample_name:           "Sample identifier; names outputs and is written into the read group"
        reference_bundle:      "Reference FASTA + BWA index from bwa_index"
        mapq_min:              "Minimum MAPQ to retain. 0 = no filter (default = 20)"
        mark_secondary:        "bwa -M: mark split hits as secondary, for Picard compatibility (default = true)"
        samtools_filter_flags: "samtools view -F. Default 0x904 drops unmapped, secondary and supplementary. Set to 0 to disable."
        extra_args:            "Additional BWA-MEM flags, passed through"
        cpu:                   "Number of CPUs delegated to task (default = 8)"
        mem_gb:                "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:               "Amount of disk space in GB delegated to task (default = 100)"
        docker:                "Container image"
    }

    meta {
        description: "Align paired-end reads with BWA-MEM, coordinate-sort, filter, and index."
    }

    command <<<
        set -euo pipefail
        (bwa 2>&1 || true) | grep -i '^Version' | sed 's/^/BWA /' | tee BWA_VERSION
        samtools --version | head -n1 | tee SAMTOOLS_VERSION

        mkdir ref
        gzip -dc ~{reference_bundle} | tar -xf - -C ref

        # Optional flags go in bash arrays so an absent flag expands to nothing
        # rather than an empty positional argument.
        bwa_optional_flags=()
        if [ "~{mark_secondary}" == "true" ]; then
            bwa_optional_flags+=("-M")
        fi

        extra_args_array=()
        if [ -n "~{extra_args}" ]; then
            read -r -a extra_args_array <<< "~{extra_args}"
        fi

        samtools_filter_flags=()
        if [ "~{samtools_filter_flags}" != "0" ]; then
            samtools_filter_flags+=("-F" "~{samtools_filter_flags}")
        fi

        bwa mem \
            -t ~{cpu} \
            "${bwa_optional_flags[@]+"${bwa_optional_flags[@]}"}" \
            -R "@RG\tID:~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA\tLB:~{sample_name}\tPU:~{sample_name}" \
            "${extra_args_array[@]+"${extra_args_array[@]}"}" \
            ref/reference.fasta \
            "~{read1}" \
            "~{read2}" | \
        samtools sort \
            -@ ~{cpu} \
            -o "~{sample_name}.sorted.bam" -

        samtools view \
            -@ ~{cpu} \
            "${samtools_filter_flags[@]+"${samtools_filter_flags[@]}"}" \
            -q ~{mapq_min} \
            -b \
            -o "~{sample_name}.sorted.filtered.bam" \
            "~{sample_name}.sorted.bam"

        samtools index "~{sample_name}.sorted.filtered.bam"
        samtools flagstat "~{sample_name}.sorted.filtered.bam" | tee "~{sample_name}.flagstat.txt"

        # Fraction of reads that align at all is the first thing to look at when
        # a sample looks wrong, and it is what tells you a sample was aligned to
        # the wrong isolate. Emit it as a value rather than making people open
        # the flagstat file.
        samtools flagstat "~{sample_name}.sorted.bam" | awk '/ mapped \(/{gsub(/[(%]/,"",$5); print $5; exit}' > PCT_MAPPED
        [ -s PCT_MAPPED ] || echo "NA" > PCT_MAPPED
    >>>

    output {
        File   bam              = "~{sample_name}.sorted.filtered.bam"
        File   bai              = "~{sample_name}.sorted.filtered.bam.bai"
        File   flagstat         = "~{sample_name}.flagstat.txt"
        String pct_mapped       = read_string("PCT_MAPPED")
        String bwa_version      = read_string("BWA_VERSION")
        String samtools_version = read_string("SAMTOOLS_VERSION")
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
