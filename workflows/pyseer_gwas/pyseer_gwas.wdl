version 1.0

import "../../tasks/pyseer.wdl" as pyseer_tasks

# workflow pyseer_gwas
#
# Standalone microbial pangenome-wide association study using pyseer's LMM
# (likelihood-ratio test), independently callable in Terra. Accepts any
# block_id x isolate presence/absence Rtab - a Panaroo gene_presence_absence
# matrix or a pangenome-network module matrix both work unchanged - and a
# core-genome phylogeny for population-structure correction.
#
# Deliberately separate from the multimodal_gwas workflow (feature/
# multimodal-gwas), which re-implements the much larger microGWAS Snakemake
# DAG (unitigs, structural variants, panfeed, whole-genome elastic net) as
# one component of a bigger pipeline. This workflow exists to run pyseer on
# its own; multimodal_gwas can later import it once that pipeline needs the
# same tasks.

workflow pyseer_gwas {

    meta {
        description: "Bacterial GWAS using pyseer's LMM (likelihood-ratio test) over a presence/absence Rtab, with optional lineage effects, per-covariate-column scanning, and a gene/module-annotation join. Accepts Panaroo gene-family Rtabs or pangenome-network module Rtabs unchanged."
        author: "Alex Arvanitis"
    }

    input {
        File            phylogeny_newick
        File            presence_absence_rtab
        File            phenotype_tsv

        File?           variant_vcf

        # Joint covariate adjustment, applied to every result. Also see
        # covariate_combinations below for testing chosen groups of
        # covariates on their own, without loading all of them at once.
        File?           covariates_file
        String?         use_covariates
        Array[String]   covariate_combinations = []

        Boolean         run_lineage_effects = false
        File?           lineage_clusters

        # Gene/module-name annotation join, applied to gene_results and
        # gene_significant. Left unset, the workflow still runs and simply
        # skips the join (see pyseer_annotate_results).
        File?           annotation_table
        String          annotation_id_column         = "Gene"
        String          annotation_name_column        = "Non-unique Gene name"
        String          annotation_description_column = "Annotation"

        Float           min_af = 0.05
        Float           max_af = 0.95
        Int             association_cpu    = 4
        Int             association_mem_gb = 8

        String          docker = "aarvani1/pyseer:1.4.2@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c"
    }

    parameter_meta {
        phylogeny_newick:       "Midpoint-rooted core-genome phylogeny (Newick), used only to build the kinship matrix"
        presence_absence_rtab:  "block_id x isolate 0/1 matrix - Panaroo's gene_presence_absence.Rtab, or a module Rtab from a pangenome-network decomposition"
        phenotype_tsv:          "Two columns: sample_id\\tphenotype_value"
        covariates_file:        "Tab-separated: sample_id, then one named column per covariate (e.g. RST, OspC, MLST, BAPS)"
        use_covariates:         "pyseer --use-covariates value, applied jointly across every covariates_file column named here"
        covariate_combinations: "Groups of covariates_file column names to test together, one group per pyseer run, each producing its own result file under covariate_scan_results. A bare name is a single-covariate run (e.g. 'BAPS'); '+' joins several into one run (e.g. 'RST+OspC+BAPS'). Deliberately not an all-combinations search - see README for why loading every column jointly is a real degrees-of-freedom risk, not just a style choice."
        run_lineage_effects:    "Report per-lineage effects as a separate, minimal-input pyseer call, not folded into the main association (default = false). The distance matrix pyseer's --lineage requires is derived automatically from phylogeny_newick - nothing extra to supply."
        lineage_clusters:       "Two columns sample_id\\tcluster_id (e.g. BAPS) for --lineage-clusters. Omit to use pyseer's MDS-derived lineages."
        annotation_table:       "Panaroo/Roary gene_presence_absence.csv, or an equivalent module-annotation table, for the name/annotation join"
        association_mem_gb:     "pyseer association memory. 8 GB covers a Panaroo gene Rtab for a few dozen isolates; scale up for genus-level panels (default = 8)"
        docker:                 "Container image for every task in this workflow"
    }

    call pyseer_tasks.pyseer_similarity_from_phylogeny {
        input:
            phylogeny_newick = phylogeny_newick,
            docker           = docker
    }

    call pyseer_tasks.pyseer_association {
        input:
            phenotype_tsv          = phenotype_tsv,
            presence_absence_rtab  = presence_absence_rtab,
            kinship_matrix         = pyseer_similarity_from_phylogeny.kinship_matrix,
            variant_vcf            = variant_vcf,
            covariates_file        = covariates_file,
            use_covariates         = use_covariates,
            covariate_combinations = covariate_combinations,
            min_af                 = min_af,
            max_af                 = max_af,
            cpu                    = association_cpu,
            mem_gb                 = association_mem_gb,
            docker                 = docker
    }

    # Lineage effects are deliberately their own call, against a minimal
    # slice of presence_absence_rtab, not folded into pyseer_association's
    # full-scale --lmm run - see pyseer_lineage_effects's meta.description.
    if (run_lineage_effects) {
        call pyseer_tasks.pyseer_lineage_effects {
            input:
                phenotype_tsv         = phenotype_tsv,
                presence_absence_rtab = presence_absence_rtab,
                distance_matrix       = pyseer_similarity_from_phylogeny.distance_matrix,
                lineage_clusters      = lineage_clusters,
                covariates_file       = covariates_file,
                use_covariates        = use_covariates,
                docker                = docker
        }
    }

    if (defined(annotation_table)) {
        call pyseer_tasks.pyseer_annotate_results as annotate_all {
            input:
                pyseer_results    = pyseer_association.gene_results,
                annotation_table  = annotation_table,
                id_column         = annotation_id_column,
                name_column       = annotation_name_column,
                annotation_column = annotation_description_column,
                basename          = "pyseer_gene_results",
                docker            = docker
        }

        call pyseer_tasks.pyseer_annotate_results as annotate_significant {
            input:
                pyseer_results    = pyseer_association.gene_significant,
                annotation_table  = annotation_table,
                id_column         = annotation_id_column,
                name_column       = annotation_name_column,
                annotation_column = annotation_description_column,
                basename          = "pyseer_gene_significant",
                docker            = docker
        }
    }

    output {
        File         kinship_matrix          = pyseer_similarity_from_phylogeny.kinship_matrix
        File         distance_matrix         = pyseer_similarity_from_phylogeny.distance_matrix

        File         gene_results            = pyseer_association.gene_results
        File         gene_significant        = pyseer_association.gene_significant
        File         gene_patterns           = pyseer_association.gene_patterns
        File         significance_threshold  = pyseer_association.significance_threshold
        File         gene_log                = pyseer_association.gene_log

        File?        snp_results             = pyseer_association.snp_results
        File?        snp_patterns            = pyseer_association.snp_patterns
        File?        snp_log                 = pyseer_association.snp_log

        File?        lineage_effects         = pyseer_lineage_effects.lineage_effects
        File?        lineage_pass_log        = pyseer_lineage_effects.lineage_pass_log
        Array[File]  covariate_scan_results  = pyseer_association.covariate_scan_results
        File         covariate_scan_combined = pyseer_association.covariate_scan_combined

        File?        gene_results_annotated       = annotate_all.annotated
        File?        gene_results_readable        = annotate_all.readable
        File?        gene_significant_annotated   = annotate_significant.annotated
        File?        gene_significant_readable    = annotate_significant.readable
    }
}
