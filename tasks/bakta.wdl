version 1.0

task bakta {

    input {
        File     assembly
        String   sample_name
        File?    bakta_db
        File?    proteins
        String   genus                = "Pseudomonas"
        String   species              = "aeruginosa"
        String?  strain
        String?  locus_tag
        Boolean  skip_plot            = true
        Int      cpu                  = 8
        Int      mem_gb               = 32
        Int      disk_gb              = 100
        String   docker               = "staphb/bakta:1.12.0-6.0-light@sha256:9e6684870e21b5195addc5d17589d3a50aacb11441edba94b97f26bb464fff2d"
    }

    parameter_meta {
        assembly:            "Assembly FASTA to annotate"
        sample_name:         "Some identifier for naming outputs"
        bakta_db:            "Optional .tar.gz of the FULL Bakta database. Omit to use the light database baked into the staphb image at /db/db-light, which is smaller and faster but annotates fewer CDS confidently."
        proteins:            "FASTA of trusted protein sequences for first-pass CDS assignment. Supplying the PAO1 proteome here transfers PAO1 gene names and products onto this isolate's genes, which is what makes per-isolate annotations comparable to each other and to reference-keyed omics data."
        genus:               "Genus name for the annotation (default = Pseudomonas)"
        species:             "Species name for the annotation (default = aeruginosa)"
        strain:              "Optional strain name; defaults to sample_name"
        locus_tag:           "Locus tag prefix. Defaults to a sanitised sample_name so every isolate gets deterministic, mutually distinguishable tags rather than Bakta's random prefix."
        skip_plot:           "Skip circular genome plots, which are slow and not used downstream (default = true)"
        cpu:                 "Number of CPUs delegated to task (default = 8)"
        mem_gb:              "Amount of memory in GB delegated to task (default = 32)"
        disk_gb:             "Amount of disk space in GB delegated to task (default = 100)"
        docker:              "Container image"
    }

    command <<<
        set -euo pipefail

        # The staphb light image ships its database at /db/db-light. A full
        # database has to be supplied; it is far too large to bake in, so it
        # follows the same stage-in-GCS pattern as checkm2_db.
        DB_TARBALL="~{select_first([bakta_db, ''])}"
        if [ -n "${DB_TARBALL}" ]; then
            echo "Extracting supplied Bakta database"
            mkdir -p bakta_db
            tar -xzf "${DB_TARBALL}" -C bakta_db
            DB_PATH="$(dirname "$(find bakta_db -name 'version.json' | head -n1)")"
            if [ -z "${DB_PATH}" ]; then
                echo "ERROR: no version.json found in the supplied database tarball" >&2
                exit 1
            fi
        else
            echo "WARNING: no bakta_db supplied; using the light database in the image." >&2
            DB_PATH="/db/db-light"
        fi
        echo "Bakta database: ${DB_PATH}"

        # Bakta wants an alphanumeric locus tag prefix. Derive a deterministic
        # one from the sample name so tags stay stable across reruns and remain
        # traceable to their isolate once tables from several isolates are
        # concatenated.
        LT="~{select_first([locus_tag, ''])}"
        if [ -z "${LT}" ]; then
            LT="$(echo "~{sample_name}" | tr -cd '[:alnum:]' | tr '[:lower:]' '[:upper:]' | cut -c1-12)"
            [ -z "${LT}" ] && LT="ISOLATE"
        fi
        echo "Locus tag prefix: ${LT}"

        # Bakta does not read gzipped FASTA.
        case "~{assembly}" in
            *.gz) gunzip -c ~{assembly} > input.fasta ;;
            *)    cp ~{assembly} input.fasta ;;
        esac

        # --keep-contig-headers is NOT optional and is deliberately not exposed
        # as an input. Without it Bakta renames contigs to contig_1, contig_2...,
        # the GFF then shares no sequence names with the BAM, and every
        # downstream intersect returns an empty table that looks exactly like
        # "no methylation found" rather than an error.
        EXTRA_ARGS=(--keep-contig-headers)
        ~{if skip_plot then "EXTRA_ARGS+=(--skip-plot)" else ""}
        ~{if defined(proteins) then "EXTRA_ARGS+=(--proteins " + proteins + ")" else ""}
        ~{if defined(strain) then "EXTRA_ARGS+=(--strain '" + strain + "')" else ""}

        bakta \
            --db "${DB_PATH}" \
            --output bakta_out \
            --prefix ~{sample_name} \
            --genus "~{genus}" \
            --species "~{species}" \
            --locus-tag "${LT}" \
            --threads ~{cpu} \
            --force \
            "${EXTRA_ARGS[@]}" \
            input.fasta

        for ext in gff3 gbff faa ffn fna tsv txt json; do
            if [ -f "bakta_out/~{sample_name}.${ext}" ]; then
                cp "bakta_out/~{sample_name}.${ext}" "~{sample_name}.${ext}"
            fi
        done

        # Guard the thing that silently breaks everything downstream: the GFF's
        # sequence names must match the assembly's contig names, or the
        # methylation intersect returns an empty table that looks exactly like
        # "no methylation found".
        grep -v '^#' ~{sample_name}.gff3 | cut -f1 | sort -u > gff_contigs.txt
        grep '^>' input.fasta | sed 's/^>//' | awk '{print $1}' | sort -u > fasta_contigs.txt
        SHARED=$(comm -12 gff_contigs.txt fasta_contigs.txt | wc -l)
        if [ "${SHARED}" -eq 0 ]; then
            echo "ERROR: Bakta GFF sequence names do not match the assembly contig names." >&2
            echo "  GFF:      $(head -3 gff_contigs.txt   | tr '\n' ' ')" >&2
            echo "  assembly: $(head -3 fasta_contigs.txt | tr '\n' ' ')" >&2
            echo "  Is --keep-contig-headers set?" >&2
            exit 1
        fi
        echo "Contig names agree (${SHARED} shared)"

        pull() {
            awk -F': ' -v key="$1" '$1 ~ key {gsub(/^[ \t]+/, "", $2); print $2; found=1; exit}
                 END {if (!found) print "NA"}' ~{sample_name}.txt
        }

        pull "^CDSs"          > N_CDS
        pull "^tRNAs"         > N_TRNA
        pull "^rRNAs"         > N_RRNA
        pull "^Contigs"       > N_CONTIGS
        pull "^Genome size"   > GENOME_SIZE

        bakta --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1 > VERSION
    >>>

    output {
        File    gff3            = "~{sample_name}.gff3"
        File    gbff            = "~{sample_name}.gbff"
        File    faa             = "~{sample_name}.faa"
        File    ffn             = "~{sample_name}.ffn"
        File    fna             = "~{sample_name}.fna"
        File    annotation_tsv  = "~{sample_name}.tsv"
        File    summary_txt     = "~{sample_name}.txt"
        String  n_cds           = read_string("N_CDS")
        String  n_trna          = read_string("N_TRNA")
        String  n_rrna          = read_string("N_RRNA")
        String  n_contigs       = read_string("N_CONTIGS")
        String  genome_size     = read_string("GENOME_SIZE")
        String  bakta_version   = read_string("VERSION")
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
