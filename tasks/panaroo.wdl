version 1.0

task panaroo {

    input {
        Array[File]  gff3s
        Array[File]  fnas
        Array[String]? sample_names
        String       clean_mode           = "strict"
        Float        core_threshold       = 0.95
        Float        seq_id               = 0.95
        Float        family_threshold     = 0.7
        String       refind_mode          = "default"
        Boolean      merge_paralogs       = false
        Boolean      remove_invalid_genes = true
        String?      alignment
        String?      aligner
        Int          cpu                  = 16
        Int          mem_gb               = 64
        Int          disk_gb              = 200
        Int          preemptible          = 0
        String       docker               = "staphb/panaroo:1.7.0@sha256:0e163bed1b5228bb0cc2bc07b406e9061f2fdaf3ce7b83a574fc0a49ad55ecf1"
    }

    parameter_meta {
        gff3s:                "Per-isolate Bakta GFF3s. Panaroo needs each GFF to carry its own sequence; if the ##FASTA block is absent it is appended from the matching entry in fnas."
        fnas:                 "Per-isolate assembly FASTAs from Bakta, positionally matched to gff3s"
        sample_names:         "Optional isolate names, positionally matched to gff3s. Panaroo labels each isolate's column in gene_presence_absence.csv after its GFF filename; supplying names makes the label independent of the filename. Omitted, the filename (minus .gff3) is used."
        clean_mode:           "Panaroo's error-correction aggressiveness: strict, moderate, or sensitive. Complete ONT assemblies justify 'strict'; loosen it only for fragmented or contaminated input (default = strict)"
        core_threshold:       "Fraction of isolates a gene must appear in to count as core (default = 0.95)"
        seq_id:               "Sequence identity threshold for initial clustering, panaroo -c (default = 0.95)"
        family_threshold:     "Family-level sequence identity threshold, panaroo -f (default = 0.7)"
        refind_mode:          "Gene refinding aggressiveness: default, strict, or off (default = default)"
        merge_paralogs:       "Collapse paralogous families into one group. Set per run. Against: real paralogues (e.g. efflux systems) are fused into one group, averaging away per-gene differences. For: Panaroo can over-split, scattering one gene across several groups so that isolates look like they lack it, which appears as NA in an ortholog matrix. Over-splitting is worst on large, mixed-quality panels. Default matches Panaroo's own (false)."
        remove_invalid_genes: "Drop gene calls failing basic validity checks (default = true)"
        alignment:            "Optionally 'core' or 'pan' to also emit gene alignments. Unset by default: gene_presence_absence.csv does not need the aligner."
        aligner:              "Aligner to use when alignment is set: mafft, prank, or clustal. Ignored otherwise."
        cpu:                  "Number of CPUs delegated to task (default = 16)"
        mem_gb:               "Memory in GB. Sized for a few dozen complete assemblies. Scale up for genus-level panels: graph construction is memory-hungry, and running out mid-build wastes the whole run (default = 64)"
        disk_gb:              "Amount of disk space in GB delegated to task (default = 200)"
        preemptible:          "Preemptible attempts. Defaults to 0 because this is the panel-wide aggregation step and a preemption late in a long graph build wastes the whole run (default = 0)"
        docker:               "Container image"
    }

    meta {
        description: "Build a pangenome across the isolate panel and the ortholog groups that let independently assembled genomes be compared. Independent assemblies have independent locus tags, so orthology, not coordinates, is the cross-isolate join key."
    }

    command <<<
        set -euo pipefail

        mkdir -p gffs
        GFFS=(~{sep=' ' gff3s})
        FNAS=(~{sep=' ' fnas})
        NAMES=(~{sep=' ' select_first([sample_names, []])})

        if [ "${#NAMES[@]}" -ne 0 ] && [ "${#NAMES[@]}" -ne "${#GFFS[@]}" ]; then
            echo "ERROR: ${#GFFS[@]} GFFs but ${#NAMES[@]} sample_names; these must be matched." >&2
            exit 1
        fi

        if [ "${#GFFS[@]}" -ne "${#FNAS[@]}" ]; then
            echo "ERROR: ${#GFFS[@]} GFFs but ${#FNAS[@]} FASTAs; these must be matched." >&2
            exit 1
        fi

        # Panaroo reads the nucleotide sequence out of the GFF itself. Bakta may
        # or may not append a ##FASTA block depending on version and flags, so
        # normalise rather than assume: a GFF without sequence yields a cryptic
        # downstream failure rather than a clear one.
        for i in "${!GFFS[@]}"; do
            g="${GFFS[$i]}"
            f="${FNAS[$i]}"
            if [ "${#NAMES[@]}" -ne 0 ]; then
                base="${NAMES[$i]}"
            else
                base="$(basename "${g}" .gff3)"
            fi
            if grep -q '^##FASTA' "${g}"; then
                cp "${g}" "gffs/${base}.gff"
            else
                echo "Appending sequence to ${base}"
                { cat "${g}"; echo "##FASTA"; cat "${f}"; } > "gffs/${base}.gff"
            fi
        done

        echo "Prepared $(find gffs -name '*.gff' | wc -l) annotated genomes"

        # Panaroo accepts a file of paths, which sidesteps shell glob expansion
        # and argument-length limits once the panel grows.
        find "$(pwd)/gffs" -name '*.gff' | sort > local_gffs.txt

        EXTRA_ARGS=()
        ~{if defined(alignment) then "EXTRA_ARGS+=(--alignment " + alignment + ")" else ""}
        ~{if defined(aligner)   then "EXTRA_ARGS+=(--aligner "   + aligner   + ")" else ""}
        ~{if merge_paralogs       then "EXTRA_ARGS+=(--merge_paralogs)"       else ""}
        ~{if remove_invalid_genes then "EXTRA_ARGS+=(--remove-invalid-genes)" else ""}

        panaroo \
            -i local_gffs.txt \
            -o panaroo_out \
            --clean-mode ~{clean_mode} \
            --core_threshold ~{core_threshold} \
            -c ~{seq_id} \
            -f ~{family_threshold} \
            --refind-mode ~{refind_mode} \
            -t ~{cpu} \
            "${EXTRA_ARGS[@]}"

        cp panaroo_out/gene_presence_absence.csv       ./gene_presence_absence.csv
        cp panaroo_out/gene_presence_absence.Rtab      ./gene_presence_absence.Rtab
        cp panaroo_out/summary_statistics.txt          ./summary_statistics.txt
        [ -f panaroo_out/pan_genome_reference.fa ] && cp panaroo_out/pan_genome_reference.fa ./pan_genome_reference.fa || true

        # summary_statistics.txt is "label\tdefinition\tcount" per row.
        pull() {
            awk -F'\t' -v key="$1" '$1==key {print $NF; found=1; exit}
                 END {if (!found) print "NA"}' summary_statistics.txt
        }
        pull "Core genes"       > N_CORE
        pull "Total genes"      > N_TOTAL

        panaroo --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1 > VERSION
    >>>

    output {
        File    gene_presence_absence      = "gene_presence_absence.csv"
        File    gene_presence_absence_rtab = "gene_presence_absence.Rtab"
        File    summary_statistics         = "summary_statistics.txt"
        File?   pan_genome_reference       = "pan_genome_reference.fa"
        String  n_core_genes               = read_string("N_CORE")
        String  n_total_genes              = read_string("N_TOTAL")
        String  panaroo_version            = read_string("VERSION")
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    preemptible
        maxRetries:     1
    }
}
