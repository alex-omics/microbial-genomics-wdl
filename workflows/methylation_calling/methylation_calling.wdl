version 1.0

import "../../tasks/align_modbam.wdl" as align_task
import "../../tasks/modkit.wdl" as modkit_task
import "../../tasks/bakta.wdl" as bakta_task
import "../../tasks/annotate_methylation.wdl" as annotate_task
import "../../tasks/panaroo.wdl" as panaroo_task
import "../../tasks/methylation_orthologs.wdl" as ortholog_task

workflow methylation_calling {

    meta {
        description: "Per-isolate bacterial methylation calling from ONT modified-basecalled BAMs. Each isolate's reads are mapped to its OWN assembly, pileup'd with modkit, scanned for methylated motifs, and joined to a Bakta annotation of that same assembly. Self-mapping is deliberate: for an organism with substantial accessory genome and frequent rearrangement, no single reference is adequate, and motif discovery against a foreign reference reads sequence context that the isolate does not actually have."
        author: "Alex Arvanitis"
    }

    input {
        Array[File]     modbams
        Array[File]     assemblies
        Array[String]?  sample_names

        # Trusted proteins for Bakta's first-pass CDS assignment. Supplying the
        # PAO1 proteome transfers PAO1 gene names and products onto each
        # isolate's own genes, which is what makes independently-annotated
        # isolates comparable without forcing them into shared coordinates.
        File?           proteins
        File?           bakta_db

        Boolean         run_bakta         = true
        Boolean         run_find_motifs   = true
        Boolean         run_pangenome     = true

        String          panaroo_clean_mode = "strict"
        Float           core_threshold     = 0.95

        # modkit
        Boolean         no_filtering      = false
        Float?          filter_threshold
        Int             min_coverage      = 10
        Float           min_percent       = 50.0

        # Self-mapping should comfortably exceed this; the floor exists to catch
        # modbams and assemblies passed in mismatched order.
        Float           min_mapped_percent = 85.0

        # Annotation
        Int             flank_upstream     = 300
        Boolean         trim_to_intergenic = true
        String          feature_type       = "CDS"

        String          genus              = "Pseudomonas"
        String          species            = "aeruginosa"

        String          summary_basename   = "methylation_summary"
        String          merged_basename    = "methylation_calls_long"
    }

    parameter_meta {
        modbams:            "Modified-basecalled BAMs, one per isolate, carrying MM/ML tags"
        assemblies:         "Each isolate's own assembly, positionally matched to modbams. Produced upstream (TheiaProk ONT, Autocycler, etc.) — this workflow does not assemble."
        sample_names:       "Optional labels, positionally matched to modbams. If omitted, names are derived from each filename."
        proteins:           "FASTA of trusted proteins for Bakta --proteins, e.g. the PAO1 proteome"
        bakta_db:           "Optional .tar.gz of the full Bakta database. Omitted, the light database baked into the staphb image is used."
        run_bakta:          "Annotate each assembly with Bakta. Turn off only if supplying annotations another way (default = true)"
        run_find_motifs:    "Per-isolate de novo motif discovery. The motif inventory is a primary characterisation axis here, effectively reporting which restriction-modification systems each isolate carries (default = true)"
        run_pangenome:      "Build a pangenome across the panel and collapse methylation onto ortholog groups. This is what makes independently-assembled isolates comparable; without it the outputs are a per-isolate catalogue only (default = true)"
        panaroo_clean_mode: "Panaroo error-correction mode: strict, moderate, or sensitive. Complete ONT assemblies justify 'strict' (default = strict)"
        core_threshold:     "Fraction of isolates a gene must appear in to be called core (default = 0.95)"
        min_coverage:       "Minimum Nvalid_cov for a site to be counted or annotated (default = 10)"
        min_percent:        "Minimum percent-modified for a site to count as methylated in summary stats (default = 50.0)"
        min_mapped_percent: "Mapping-rate floor, as a guard against mismatched modbam/assembly pairs (default = 85.0)"
        flank_upstream:     "Bases upstream of each CDS treated as putative promoter region (default = 300)"
        trim_to_intergenic: "Trim upstream windows that run into neighbouring genes (default = true)"
    }

    scatter (i in range(length(modbams))) {

        String resolved_name = if defined(sample_names)
            then select_first([sample_names])[i]
            else sub(basename(modbams[i]), "\\.(bam|modbam)$", "")

        # Everything below is in this isolate's own coordinate space.
        call align_task.align_modbam {
            input:
                modbam             = modbams[i],
                sample_name        = resolved_name,
                reference_fasta    = assemblies[i],
                min_mapped_percent = min_mapped_percent
        }

        call modkit_task.modkit_pileup {
            input:
                aligned_bam       = align_modbam.aligned_bam,
                aligned_bam_index = align_modbam.aligned_bam_index,
                sample_name       = resolved_name,
                reference_fasta   = assemblies[i],
                no_filtering      = no_filtering,
                filter_threshold  = filter_threshold,
                min_coverage      = min_coverage,
                min_percent       = min_percent
        }

        if (run_find_motifs) {
            call modkit_task.modkit_find_motifs {
                input:
                    bedmethyl       = modkit_pileup.bedmethyl,
                    sample_name     = resolved_name,
                    reference_fasta = assemblies[i]
            }
        }

        if (run_bakta) {
            call bakta_task.bakta {
                input:
                    assembly    = assemblies[i],
                    sample_name = resolved_name,
                    bakta_db    = bakta_db,
                    proteins    = proteins,
                    genus       = genus,
                    species     = species,
                    strain      = resolved_name
            }

            call annotate_task.annotate_methylation {
                input:
                    bedmethyl          = modkit_pileup.bedmethyl,
                    sample_name        = resolved_name,
                    reference_gff      = bakta.gff3,
                    reference_fai      = align_modbam.reference_fai,
                    flank_upstream     = flank_upstream,
                    trim_to_intergenic = trim_to_intergenic,
                    min_coverage       = min_coverage,
                    feature_type       = feature_type
            }
        }

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
            if defined(modkit_find_motifs.n_motifs) then "~{modkit_find_motifs.n_motifs}" else "NA",
            if defined(modkit_find_motifs.motifs)   then "~{modkit_find_motifs.motifs}"   else "NA",
            if defined(bakta.n_contigs)             then "~{bakta.n_contigs}"             else "NA",
            if defined(bakta.genome_size)           then "~{bakta.genome_size}"           else "NA",
            if defined(bakta.n_cds)                 then "~{bakta.n_cds}"                 else "NA",
            if defined(annotate_methylation.n_genic)    then "~{annotate_methylation.n_genic}"    else "NA",
            if defined(annotate_methylation.n_upstream) then "~{annotate_methylation.n_upstream}" else "NA",
            if defined(annotate_methylation.n_loci_with_methylation)
                then "~{annotate_methylation.n_loci_with_methylation}" else "NA",
            modkit_pileup.modkit_version
        ]
    }

    Array[String] summary_header = [
        "sample", "percent_mapped", "mean_depth", "positions_covered",
        "n_sites_6mA", "n_sites_4mC", "n_sites_5mC",
        "mean_percent_6mA", "mean_percent_4mC", "mean_percent_5mC",
        "n_motifs", "motifs",
        "n_contigs", "genome_size", "n_cds",
        "n_genic_sites", "n_upstream_sites", "n_loci_with_methylation",
        "modkit_version"
    ]

    call name_summary {
        input:
            rows     = flatten([[summary_header], summary_row]),
            basename = summary_basename
    }

    # Rows carry each isolate's own locus tags, so this concatenation is a
    # per-isolate catalogue rather than a cross-isolate matrix. Comparing
    # isolates to each other needs an orthology layer to map those tags onto
    # shared gene groups; coordinates cannot do that job across independent
    # assemblies.
    if (run_bakta) {
        call merge_annotated {
            input:
                annotated_tables = select_all(annotate_methylation.annotated_tsv),
                basename         = merged_basename
        }
    }

    # Orthology is the cross-isolate join key. Coordinates cannot serve that
    # role once every isolate has been assembled and annotated independently,
    # so the pangenome is what turns a stack of per-isolate catalogues into a
    # gene-by-isolate matrix.
    if (run_bakta && run_pangenome) {
        call panaroo_task.panaroo {
            input:
                gff3s          = select_all(bakta.gff3),
                fnas           = select_all(bakta.fna),
                clean_mode     = panaroo_clean_mode,
                core_threshold = core_threshold
        }

        call ortholog_task.methylation_orthologs {
            input:
                gene_presence_absence = panaroo.gene_presence_absence,
                annotated_tables      = select_all(annotate_methylation.annotated_tsv),
                sample_names          = resolved_name
        }
    }

    output {
        File  summary_tsv          = name_summary.summary
        File? methylation_long     = merge_annotated.merged

        # The comparative deliverables: ortholog groups x isolates, with gene
        # absence held as NA rather than collapsed to zero.
        File? ortholog_matrix      = methylation_orthologs.matrix
        File? ortholog_long        = methylation_orthologs.long_table
        File? pangenome_presence   = panaroo.gene_presence_absence
        File? pangenome_summary    = panaroo.summary_statistics

        Array[File] aligned_bams         = align_modbam.aligned_bam
        Array[File] aligned_bam_indexes  = align_modbam.aligned_bam_index
        Array[File] flagstats            = align_modbam.flagstat
        Array[File] coverage_reports     = align_modbam.coverage_txt

        Array[File] bedmethyl_gz         = modkit_pileup.bedmethyl_gz
        Array[File] bed_6ma              = modkit_pileup.bed_6ma
        Array[File] bed_4mc              = modkit_pileup.bed_4mc
        Array[File] bed_5mc              = modkit_pileup.bed_5mc

        Array[File?] motif_tables        = modkit_find_motifs.motifs_tsv
        Array[File?] bakta_gff3          = bakta.gff3
        Array[File?] bakta_faa           = bakta.faa
        Array[File?] bakta_tsv           = bakta.annotation_tsv
        Array[File?] annotated_tables    = annotate_methylation.annotated_tsv
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

task merge_annotated {
    input {
        Array[File] annotated_tables
        String      basename
    }

    command <<<
        set -euo pipefail

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
