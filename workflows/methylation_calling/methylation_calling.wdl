version 1.0

import "../../tasks/align_modbam.wdl" as align_task
import "../../tasks/modkit.wdl" as modkit_task
import "../../tasks/annotate_methylation.wdl" as annotate_task

workflow methylation_calling {

    meta {
        description: "Per-sample bacterial methylation calling from ONT modified-basecalled BAMs. Aligns modBAMs to a common reference, pileups 6mA/4mC/5mC with modkit, joins the calls to the reference annotation, and collapses everything into one long-format table keyed by locus tag for comparative methylomics and multi-omics integration."
        author: "Alex Arvanitis"
    }

    input {
        Array[File]     modbams
        Array[String]?  sample_names

        # Every sample must go against the SAME reference. Cross-sample
        # comparison is coordinate-based, so a per-isolate reference would give
        # each sample its own incompatible coordinate space.
        File            reference_fasta
        File            reference_gff

        Boolean         run_find_motifs   = true

        # modkit
        Boolean         no_filtering      = false
        Float?          filter_threshold
        Int             min_coverage      = 10
        Float           min_percent       = 50.0

        # Annotation
        Int             flank_upstream     = 300
        Boolean         trim_to_intergenic = true
        String          feature_type       = "CDS"

        String          summary_basename   = "methylation_summary"
        String          merged_basename    = "methylation_calls_long"
    }

    parameter_meta {
        modbams:            "Modified-basecalled BAMs, one per sample, carrying MM/ML tags. Unaligned is the expected input; already-aligned files are reduced to primary records and realigned."
        sample_names:       "Optional labels, positionally matched to modbams. If omitted, names are derived from each filename."
        reference_fasta:    "Common reference for all samples (e.g. PAO1). Comparative methylomics requires one shared coordinate space."
        reference_gff:      "Annotation matching reference_fasta. Use the reference's own curated GFF so locus tags join cleanly to transcriptomics and proteomics."
        run_find_motifs:    "Run de novo motif discovery per sample. Strongly recommended for multi-isolate studies: the motif inventory reports which restriction-modification systems each isolate carries, which is the control needed before calling any cross-isolate difference regulatory (default = true)"
        min_coverage:       "Minimum Nvalid_cov for a site to be counted or annotated (default = 10)"
        min_percent:        "Minimum percent-modified for a site to count as methylated in summary stats (default = 50.0)"
        flank_upstream:     "Bases upstream of each CDS treated as putative promoter region (default = 300)"
        trim_to_intergenic: "Trim upstream windows that run into neighbouring genes (default = true)"
    }

    # Indexing once and passing the .fai around beats re-indexing the reference
    # inside every scattered task.
    call index_reference {
        input:
            reference_fasta = reference_fasta
    }

    scatter (i in range(length(modbams))) {

        String resolved_name = if defined(sample_names)
            then select_first([sample_names])[i]
            else sub(basename(modbams[i]), "\\.(bam|modbam)$", "")

        call align_task.align_modbam {
            input:
                modbam          = modbams[i],
                sample_name     = resolved_name,
                reference_fasta = reference_fasta
        }

        call modkit_task.modkit_pileup {
            input:
                aligned_bam       = align_modbam.aligned_bam,
                aligned_bam_index = align_modbam.aligned_bam_index,
                sample_name       = resolved_name,
                reference_fasta   = reference_fasta,
                no_filtering      = no_filtering,
                filter_threshold  = filter_threshold,
                min_coverage      = min_coverage,
                min_percent       = min_percent
        }

        call annotate_task.annotate_methylation {
            input:
                bedmethyl          = modkit_pileup.bedmethyl,
                sample_name        = resolved_name,
                reference_gff      = reference_gff,
                reference_fai      = index_reference.fai,
                flank_upstream     = flank_upstream,
                trim_to_intergenic = trim_to_intergenic,
                min_coverage       = min_coverage,
                feature_type       = feature_type
        }

        if (run_find_motifs) {
            call modkit_task.modkit_find_motifs {
                input:
                    bedmethyl       = modkit_pileup.bedmethyl,
                    sample_name     = resolved_name,
                    reference_fasta = reference_fasta
            }
        }

        # One row per sample. Everything stringified and NA-defaulted so the
        # table stays rectangular whether or not motif discovery ran.
        Array[String] summary_row = [
            resolved_name,
            align_modbam.percent_mapped,
            align_modbam.mean_depth,
            "~{modkit_pileup.n_positions_covered}",
            "~{modkit_pileup.n_sites_6ma}",
            "~{modkit_pileup.n_sites_4mc}",
            "~{modkit_pileup.n_sites_5mc}",
            modkit_pileup.mean_percent_6ma,
            modkit_pileup.mean_percent_4mc,
            modkit_pileup.mean_percent_5mc,
            "~{annotate_methylation.n_genic}",
            "~{annotate_methylation.n_upstream}",
            "~{annotate_methylation.n_loci_with_methylation}",
            if defined(modkit_find_motifs.n_motifs) then "~{modkit_find_motifs.n_motifs}" else "NA",
            if defined(modkit_find_motifs.motifs)   then "~{modkit_find_motifs.motifs}"   else "NA",
            modkit_pileup.modkit_version
        ]
    }

    Array[String] summary_header = [
        "sample", "percent_mapped", "mean_depth", "positions_covered",
        "n_sites_6mA", "n_sites_4mC", "n_sites_5mC",
        "mean_percent_6mA", "mean_percent_4mC", "mean_percent_5mC",
        "n_genic_sites", "n_upstream_sites", "n_loci_with_methylation",
        "n_motifs", "motifs", "modkit_version"
    ]

    call name_summary {
        input:
            rows     = flatten([[summary_header], summary_row]),
            basename = summary_basename
    }

    # The long-format table is the actual deliverable for comparative work:
    # one row per (sample, site, feature), keyed by locus tag so it joins
    # directly against expression or abundance matrices.
    call merge_annotated {
        input:
            annotated_tables = annotate_methylation.annotated_tsv,
            basename         = merged_basename
    }

    output {
        File summary_tsv        = name_summary.summary
        File methylation_long   = merge_annotated.merged

        Array[File] aligned_bams        = align_modbam.aligned_bam
        Array[File] aligned_bam_indexes = align_modbam.aligned_bam_index
        Array[File] flagstats           = align_modbam.flagstat

        Array[File] bedmethyl_gz        = modkit_pileup.bedmethyl_gz
        Array[File] bed_6ma             = modkit_pileup.bed_6ma
        Array[File] bed_4mc             = modkit_pileup.bed_4mc
        Array[File] bed_5mc             = modkit_pileup.bed_5mc

        Array[File] annotated_tables    = annotate_methylation.annotated_tsv

        Array[File?] motif_tables       = modkit_find_motifs.motifs_tsv
    }
}

# Same helper as assembly_qc.wdl: write_tsv alone emits a temp-named file that
# lands in the bucket as something like "tmpzsaunsvb", so route it through a
# task to give the deliverable a recognisable name. Duplicated rather than
# imported because it is defined inside assembly_qc.wdl; worth lifting into a
# shared tasks/utils.wdl if a third workflow needs it.
task name_summary {
    input {
        Array[Array[String]] rows
        String               basename
    }

    command <<<
        set -euo pipefail
        cp ~{write_tsv(rows)} "~{basename}.tsv"
    >>>

    output {
        File summary = "~{basename}.tsv"
    }

    runtime {
        docker:         "ubuntu:22.04"
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     2
    }
}

task index_reference {
    input {
        File reference_fasta
    }

    command <<<
        set -euo pipefail
        cp ~{reference_fasta} ref.fa
        samtools faidx ref.fa
    >>>

    output {
        File fai = "ref.fa.fai"
    }

    runtime {
        docker:         "staphb/samtools:1.24"
        memory:         "4 GB"
        cpu:            1
        disks:          "local-disk 20 SSD"
        preemptible:    1
        maxRetries:     2
    }
}

task merge_annotated {
    input {
        Array[File] annotated_tables
        String      basename
    }

    command <<<
        set -euo pipefail

        # Header from the first table, then every table's body. Each row already
        # carries its sample name, so the result is tidy long format.
        head -n1 ~{annotated_tables[0]} > "~{basename}.tsv"
        for f in ~{sep=' ' annotated_tables}; do
            awk 'NR>1' "${f}" >> "~{basename}.tsv"
        done

        gzip -c "~{basename}.tsv" > "~{basename}.tsv.gz"
        awk 'NR>1' "~{basename}.tsv" | wc -l > N_ROWS
    >>>

    output {
        File merged     = "~{basename}.tsv"
        File merged_gz  = "~{basename}.tsv.gz"
        Int  n_rows     = read_int("N_ROWS")
    }

    runtime {
        docker:         "ubuntu:22.04"
        memory:         "8 GB"
        cpu:            1
        disks:          "local-disk 50 SSD"
        preemptible:    1
        maxRetries:     2
    }
}
