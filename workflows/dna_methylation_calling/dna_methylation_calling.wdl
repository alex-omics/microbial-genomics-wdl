version 1.0

import "../../tasks/align_modbam.wdl" as align_task
import "../../tasks/modkit.wdl" as modkit_task
import "../../tasks/validate_annotation.wdl" as validate_task
import "../../tasks/annotate_methylation.wdl" as annotate_task
import "../../tasks/panaroo.wdl" as panaroo_task
import "../../tasks/methylation_orthologs.wdl" as ortholog_task
import "../../tasks/rebase_mtase_search.wdl" as rebase_task
import "../../tasks/motif_landscape.wdl" as landscape_task
import "../../tasks/motif_landscape_summary.wdl" as landscape_summary_task
import "../../tasks/annotate_motifs.wdl" as annotate_motifs_task
import "../../tasks/utils.wdl" as utils

workflow dna_methylation_calling {

    meta {
        description: "Per-isolate bacterial methylation calling from ONT modified-basecalled BAMs. Each isolate's reads are mapped to its own assembly, piled up with modkit, scanned for methylated motifs, and joined to a supplied annotation of the same assembly. Works for any bacterium for which the assemblies and their annotations are available; nothing in it is specific to a species. Self-mapping is deliberate: with a large accessory genome and frequent rearrangement no single reference is adequate, and motif discovery against a foreign reference reads sequence context the isolate does not have. This workflow does not assemble or annotate."
        author: "Alex Arvanitis"
    }

    input {
        Array[File]     modbams
        Array[File]     assemblies

        # Each isolate's annotation of its own assembly, from whatever annotator
        # produced it (e.g. the Bakta output of TheiaProk). All three are
        # positionally matched to modbams/assemblies, and checked against the
        # assembly before anything expensive runs.
        Array[File]     annotation_gff3s
        Array[File]     annotation_faas
        Array[File]     annotation_fnas

        Array[String]?  sample_names

        Boolean         run_find_motifs   = true
        Boolean         run_pangenome     = true

        # REBASE homology identification. Both files must be supplied to run it.
        # They are inputs rather than part of the image or repo because REBASE
        # carries its own redistribution terms.
        File?           rebase_goldset_fasta
        File?           rebase_motif_tsv
        String          rebase_evalue      = "1e-25"

        # Tier 1/2/3 landscape characterisation (see motif_landscape.wdl and
        # motif_landscape_summary.wdl). Runs off find-motifs and/or REBASE
        # output, whichever is available; needs neither MICs nor phenotype
        # groups, since it is descriptive across the panel, not a case/control
        # test.
        Boolean         run_motif_landscape       = true
        Float           heterogeneous_low_cutoff  = 50.0
        Array[String]?  highlight_genes

        # Reuse a pangenome already built by workflows/pangenome rather than
        # recomputing it here. Its isolate column names must match this run's
        # sample names.
        File?           gene_presence_absence

        String          panaroo_clean_mode = "strict"
        Float           core_threshold     = 0.95
        Boolean         merge_paralogs     = false

        # modkit
        Boolean         no_filtering      = false
        Float?          filter_threshold
        Int             min_coverage      = 10
        Float           min_percent       = 50.0
        Int             min_mod_reads     = 3

        # Reads shorter than this are dropped before alignment, not after --
        # a raw modBAM commonly carries a tail of short fragments mechanically
        # incapable of a confident placement, and leaving them in makes
        # percent_mapped measure read-length composition as much as it
        # measures whether the modbam and assembly actually correspond.
        Int             min_read_length    = 500

        # Self-mapping should comfortably exceed this; the floor exists to catch
        # modbams and assemblies passed in mismatched order.
        Float           min_mapped_percent = 85.0

        # Annotation
        Int             flank_upstream     = 300
        Boolean         trim_to_intergenic = true
        String          feature_type       = "CDS"

        String          summary_basename   = "methylation_summary"
        String          merged_basename    = "methylation_calls_long"
    }

    parameter_meta {
        modbams:            "Modified-basecalled BAMs, one per isolate, carrying MM/ML tags"
        assemblies:         "Each isolate's own assembly, positionally matched to modbams. This workflow does not assemble."
        annotation_gff3s:   "Each isolate's GFF3 annotation of its own assembly, positionally matched to modbams/assemblies. Its sequence names must be the assembly's contig names. Locus tags are what the pangenome later maps onto shared ortholog groups."
        annotation_faas:    "Each isolate's protein FASTA from the same annotation, positionally matched. Used for REBASE homology."
        annotation_fnas:    "Each isolate's nucleotide FASTA from the same annotation (the assembly as the annotator wrote it), positionally matched. Used for the pangenome and to check the annotation against the assembly."
        sample_names:       "Optional labels, positionally matched to modbams. If omitted, names are derived from each filename."
        run_find_motifs:    "Per-isolate de novo motif discovery. The motif inventory is a primary characterisation axis here, effectively reporting which restriction-modification systems each isolate carries (default = true)"
        run_pangenome:      "Build a pangenome across the panel and collapse methylation onto ortholog groups. This is what makes independently assembled isolates comparable; without it the outputs are a per-isolate catalogue. The panel should be a single species: ortholog groups across species are not meaningful (default = true)"
        panaroo_clean_mode: "Panaroo error-correction mode: strict, moderate, or sensitive. Complete ONT assemblies justify 'strict' (default = strict)"
        core_threshold:     "Fraction of isolates a gene must appear in to be called core (default = 0.95)"
        min_coverage:       "Minimum Nvalid_cov for a site to be counted or annotated (default = 10)"
        min_percent:        "Minimum percent-modified for a site to count as methylated. Re-thresholding reruns only annotation and the ortholog join (default = 50.0)"
        min_mod_reads:      "Minimum reads actually carrying the modification. Decouples the percent floor from depth, so min_percent can be lowered to catch partial methylation without admitting two-read calls (default = 3)"
        min_read_length:    "Reads shorter than this are dropped before alignment, removing the short-fragment noise common in a raw modBAM while keeping shorter reads that still contribute coverage (default = 500)"
        min_mapped_percent: "Mapping-rate floor, evaluated AFTER the length filter, as a guard against mismatched modbam/assembly pairs (default = 85.0)"
        flank_upstream:     "Bases upstream of each CDS treated as putative promoter region (default = 300)"
        trim_to_intergenic: "Trim upstream windows that run into neighbouring genes (default = true)"
        rebase_goldset_fasta:      "REBASE Gold Standard protein set for MTase homology identification. Both this and rebase_motif_tsv must be supplied to run it; omitting either skips REBASE cleanly rather than failing."
        rebase_motif_tsv:          "REBASE enzyme -> recognition motif -> modification type table, keyed on rebase_goldset_fasta's headers"
        rebase_evalue:             "BLASTP e-value cutoff for a REBASE homology call, as a String -- a Float this small renders as the literal text 0.000000 in WDL's interpolation, which blastp rejects (default = \"1e-25\")"
        run_motif_landscape:       "Test de novo and REBASE candidate motifs against this isolate's own data (enrichment and within-genome heterogeneity), then summarise variability across the panel by motif and by gene. Descriptive only; needs no phenotype data (default = true)"
        heterogeneous_low_cutoff:  "Below this percent-modified, a motif occurrence counts as 'low' in the heterogeneity summary -- meaningful relative to the panel's typical housekeeping level, usually near 100 (default = 50.0)"
        highlight_genes:           "Case-insensitive substrings (e.g. ['hsd','mod','res']) to flag in the gene-level tier-3 table. A sort/flag convenience, not a filter -- every gene is still reported."
    }

    # basename() evaluates against the path string without localising anything,
    # so this scatter costs nothing and runs no tasks.
    scatter (m in modbams) {
        String derived_name = sub(basename(m), "\\.(bam|modbam)$", "")
    }

    # Gating the main scatter on this is the point: a mismatched panel or an
    # unusable sample name otherwise surfaces only after alignment and Panaroo
    # have already been paid for.
    call utils.validate_panel {
        input:
            derived_names    = derived_name,
            n_primary        = length(modbams),
            companion_counts = [length(assemblies), length(annotation_gff3s),
                                length(annotation_faas), length(annotation_fnas)],
            sample_names     = sample_names
    }

    scatter (i in range(length(modbams))) {

        String resolved_name = validate_panel.sample_ids[i]

        # Cheap, and everything below uses its checked_name, so nothing
        # expensive starts for an isolate whose assembly and annotation do not
        # belong together.
        call validate_task.validate_annotation {
            input:
                assembly    = assemblies[i],
                gff3        = annotation_gff3s[i],
                fna         = annotation_fnas[i],
                faa         = annotation_faas[i],
                sample_name = resolved_name
        }

        String isolate = validate_annotation.checked_name

        # Everything below is in this isolate's own coordinate space.
        call align_task.align_modbam {
            input:
                modbam             = modbams[i],
                sample_name        = isolate,
                reference_fasta    = assemblies[i],
                min_read_length    = min_read_length,
                min_mapped_percent = min_mapped_percent
        }

        call modkit_task.modkit_pileup {
            input:
                aligned_bam       = align_modbam.aligned_bam,
                aligned_bam_index = align_modbam.aligned_bam_index,
                sample_name       = isolate,
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
                    sample_name     = isolate,
                    reference_fasta = assemblies[i]
            }
        }

        call annotate_task.annotate_methylation {
            input:
                bedmethyl          = modkit_pileup.bedmethyl,
                sample_name        = isolate,
                reference_gff      = annotation_gff3s[i],
                reference_fai      = align_modbam.reference_fai,
                flank_upstream     = flank_upstream,
                trim_to_intergenic = trim_to_intergenic,
                min_coverage       = min_coverage,
                min_percent        = min_percent,
                min_mod_reads      = min_mod_reads,
                feature_type       = feature_type
        }

        if (defined(rebase_goldset_fasta) && defined(rebase_motif_tsv)) {
            call rebase_task.rebase_blastp {
                input:
                    faa                   = annotation_faas[i],
                    sample_name           = isolate,
                    rebase_goldset_fasta  = select_first([rebase_goldset_fasta]),
                    evalue                = rebase_evalue
            }

            call rebase_task.rebase_join_motifs {
                input:
                    blast_hits       = rebase_blastp.blast_hits,
                    sample_name      = isolate,
                    rebase_motif_tsv = select_first([rebase_motif_tsv])
            }
        }

        if (run_motif_landscape) {
            # Neither input is required -- an isolate with find-motifs off and
            # no REBASE hits just tests zero motifs rather than failing, which
            # keeps this composable with every other toggle in the workflow.
            call landscape_task.build_motif_list {
                input:
                    find_motifs_tsv   = modkit_find_motifs.motifs_tsv,
                    rebase_mtases_tsv = rebase_join_motifs.rebase_mtases,
                    basename          = "~{isolate}_motif_list"
            }

            call landscape_task.motif_landscape {
                input:
                    bedmethyl                = modkit_pileup.bedmethyl,
                    motif_list               = build_motif_list.motif_list,
                    sample_name              = isolate,
                    reference_fasta          = assemblies[i],
                    min_coverage             = min_coverage,
                    heterogeneous_low_cutoff = heterogeneous_low_cutoff
            }

            # The single per-isolate deliverable: gene/product/region annotation,
            # promoter windows, and motif context alongside each methylation
            # call, in one row per site.
            call annotate_motifs_task.annotate_motif_membership {
                input:
                    annotated_tsv   = annotate_methylation.annotated_tsv,
                    motif_list      = build_motif_list.motif_list,
                    sample_name     = isolate,
                    reference_fasta = assemblies[i]
            }
        }

        Array[String] summary_row = [
            isolate,
            "~{align_modbam.n_reads_raw}",
            "~{align_modbam.n_reads_length_filtered}",
            align_modbam.percent_mapped,
            align_modbam.mean_depth,
            "~{modkit_pileup.n_positions_covered}",
            "~{modkit_pileup.n_sites_6ma}",
            "~{modkit_pileup.n_sites_4mc}",
            "~{modkit_pileup.n_sites_5mc}",
            modkit_pileup.mean_percent_6ma,
            modkit_pileup.mean_percent_4mc,
            modkit_pileup.mean_percent_5mc,
            modkit_pileup.mod_codes_seen,
            if defined(modkit_find_motifs.n_motifs) then "~{modkit_find_motifs.n_motifs}" else "NA",
            if defined(modkit_find_motifs.motifs)   then "~{modkit_find_motifs.motifs}"   else "NA",
            "~{validate_annotation.n_contigs}",
            "~{validate_annotation.genome_size}",
            validate_annotation.gc_percent,
            "~{validate_annotation.n_cds}",
            "~{annotate_methylation.n_genic}",
            "~{annotate_methylation.n_upstream}",
            "~{annotate_methylation.n_loci_with_methylation}",
            modkit_pileup.modkit_version
        ]
    }

    Array[String] summary_header = [
        "sample", "n_reads_raw", "n_reads_length_filtered", "percent_mapped", "mean_depth", "positions_covered",
        "n_sites_6mA", "n_sites_4mC", "n_sites_5mC",
        "mean_percent_6mA", "mean_percent_4mC", "mean_percent_5mC", "mod_codes_seen",
        "n_motifs", "motifs",
        "n_contigs", "genome_size", "gc_percent", "n_cds",
        "n_genic_sites", "n_upstream_sites", "n_loci_with_methylation",
        "modkit_version"
    ]

    call utils.name_summary {
        input:
            rows     = flatten([[summary_header], summary_row]),
            basename = summary_basename
    }

    # Rows carry each isolate's own locus tags, so this concatenation is a
    # per-isolate catalogue rather than a cross-isolate matrix. Comparing
    # isolates to each other needs an orthology layer to map those tags onto
    # shared gene groups; coordinates cannot do that job across independent
    # assemblies.
    call utils.concat_tables {
        input:
            tables   = annotate_methylation.annotated_tsv,
            basename = merged_basename
    }

    # Orthology is the cross-isolate join key. Coordinates cannot serve that
    # role once every isolate has been assembled and annotated independently,
    # so the pangenome is what turns a stack of per-isolate catalogues into a
    # gene-by-isolate matrix.
    #
    # sample_names is passed explicitly: Panaroo otherwise labels each isolate's
    # column after its GFF filename, which for a supplied annotation is whatever
    # the upstream pipeline called it, and the ortholog join below looks the
    # columns up by sample name.
    if (run_pangenome && !defined(gene_presence_absence)) {
        call panaroo_task.panaroo {
            input:
                gff3s          = annotation_gff3s,
                fnas           = annotation_fnas,
                sample_names   = isolate,
                clean_mode     = panaroo_clean_mode,
                core_threshold = core_threshold,
                merge_paralogs = merge_paralogs
        }
    }

    if (run_pangenome || defined(gene_presence_absence)) {
        call ortholog_task.methylation_orthologs {
            input:
                gene_presence_absence = select_first([gene_presence_absence,
                                                      panaroo.gene_presence_absence]),
                annotated_tables      = annotate_methylation.annotated_tsv,
                sample_names          = isolate
        }
    }

    # Tier 3: cross-isolate variability, both by motif (needs nothing but the
    # per-isolate landscape tables) and by gene (reuses the ortholog long
    # table above -- if pangenome/ortholog work did not run, the gene-level
    # half is simply omitted, since the motif-level half does not depend on it).
    if (run_motif_landscape) {
        call landscape_summary_task.motif_landscape_summary {
            input:
                per_isolate_landscape = select_all(motif_landscape.landscape_tsv),
                ortholog_long_table   = methylation_orthologs.long_table,
                highlight_genes       = highlight_genes
        }
    }

    output {
        File  summary_tsv          = name_summary.summary
        File  methylation_long     = concat_tables.merged

        # The comparative deliverables: ortholog groups x isolates, with gene
        # absence held as NA rather than collapsed to zero.
        File? ortholog_matrix      = methylation_orthologs.matrix
        File? ortholog_long        = methylation_orthologs.long_table
        File? pangenome_presence   = panaroo.gene_presence_absence
        File? pangenome_summary    = panaroo.summary_statistics

        # Tier 1/2/3 landscape characterisation.
        Array[File?] rebase_mtases       = rebase_join_motifs.rebase_mtases
        Array[File?] motif_lists         = build_motif_list.motif_list
        Array[File?] motif_landscapes    = motif_landscape.landscape_tsv
        File?        motif_summary       = motif_landscape_summary.motif_summary
        File?        gene_landscape_summary = motif_landscape_summary.gene_summary

        Array[File] aligned_bams         = align_modbam.aligned_bam
        Array[File] aligned_bam_indexes  = align_modbam.aligned_bam_index
        Array[File] flagstats            = align_modbam.flagstat
        Array[File] coverage_reports     = align_modbam.coverage_txt

        Array[File] bedmethyl_gz         = modkit_pileup.bedmethyl_gz
        Array[File] bed_6ma              = modkit_pileup.bed_6ma
        Array[File] bed_4mc              = modkit_pileup.bed_4mc
        Array[File] bed_5mc              = modkit_pileup.bed_5mc

        Array[File?] motif_tables        = modkit_find_motifs.motifs_tsv
        Array[File]  annotated_tables    = annotate_methylation.annotated_tsv

        # annotate_methylation's table plus a motif column: the single
        # per-isolate file combining gene/product/region annotation, promoter
        # windows, and motif context alongside each methylation call.
        Array[File?] annotated_tables_with_motifs = annotate_motif_membership.combined_tsv
    }
}
