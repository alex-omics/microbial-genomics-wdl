version 1.0

import "../../tasks/gwas/task_setup.wdl" as setup
import "../../tasks/gwas/task_variants.wdl" as variants
import "../../tasks/gwas/task_population_structure.wdl" as popstruct
import "../../tasks/gwas/task_rare_variants.wdl" as rare
import "../../tasks/gwas/task_pyseer.wdl" as pyseer
import "../../tasks/gwas/task_mapback.wdl" as mapback
import "../../tasks/gwas/task_heritability.wdl" as heritability
import "../../tasks/gwas/task_annotation.wdl" as annotation
import "../../tasks/gwas/task_summaries.wdl" as summaries
import "../../tasks/gwas/task_enrichment.wdl" as enrichment
import "../../tasks/gwas/task_amr.wdl" as amr

# wf_multimodal_gwas.wdl
#
# Terra-native bacterial GWAS across multiple classes of genetic variation.
#
# This is a WDL translation of microGWAS
# (https://github.com/microbial-pangenomes-lab/microGWAS), the Snakemake
# pipeline described in Burgaya et al. 2025, Microbial Genomics 11:001349.
# MIT-licensed, Copyright (c) 2022, Marco Galardini. Pinned at commit 1307250
# (VERSION 0.9.1-dev). Shell commands and script invocations are preserved
# verbatim wherever possible; see the individual task files for the places
# where they are not, and why.
#
# ---------------------------------------------------------------------------
# The call graph
# ---------------------------------------------------------------------------
# Snakemake infers dependencies by matching each rule's output paths against
# other rules' input paths. WDL requires them stated. This block is that
# statement, reconstructed by reading each rule's input:/output: config-key
# references rather than inferring from rule names. Broadly:
#
#   assemblies ─┬─> ggcaller ──> GFFs ─┬─> panaroo ─┬─> core alignment ─┬─> FastTree
#               │                      │            │                   ├─> snp-sites ─> VCF ─> similarity (kinship)
#               │                      │            │                   └─> variant sites
#               │                      │            ├─> gene presence/absence
#               │                      │            └─> structural variants
#               │                      └─> panfeed ──> k-mer patterns
#               ├─> unitig-counter ──> unitigs
#               ├─> mash sketch ──> mash dist ──> distance matrix
#               ├─> mlst ──> lineages
#               └─> snippy ──> per-sample VCFs ──> bcftools ──> common/rare SNPs
#
#   then, scattered per phenotype:
#     prepare_pyseer ──> run_pyseer / run_pyseer_vcf / run_panfeed
#                   └──> heritability
#     prepare_wg ──> run_wg (ridge, lasso) ──> metrics, map-back, summaries
#     hits ──> map back ──> summaries ──> eggNOG annotation ──> enrichment ──> plots
#
# The ~12 pure "collector" rules upstream (manhattan_plots, wg, qq_plots,
# map_back, map_summary, annotate_summary, enrichment, enrichment_plots,
# pyseer, pyseer_vcf, heritability, panfeed) have no WDL equivalent by design:
# they exist only to force Snakemake to materialise every scattered output via
# expand(), and a task consuming an Array[File] from a scatter already means
# "wait for all of these".

workflow wf_multimodal_gwas {

  meta {
    version: "0.1.0"
    author: "Alex Arvanitis"
    description: "Bacterial GWAS across unitigs, gene presence/absence, gene-cluster k-mers, common SNPs and a whole-genome elastic net, with population-structure correction via pyseer's LMM. WDL translation of the microGWAS Snakemake pipeline."
  }

  parameter_meta {
    samples_tsv: "Tab-separated sample table. First column is the strain ID; must also contain one column per phenotype named in `targets`. Covariate columns are optional."
    sample_fastas: "Assembly FASTAs, one per sample, named SAMPLE.fasta to match the strain IDs in samples_tsv."
    sample_gffs: "Optional pre-computed annotations, named SAMPLE.gff. Supply these to skip gene calling; omit to run ggCaller."
    targets: "Phenotype column names in samples_tsv to run associations for. Each becomes an independent scatter."
    covariates: "Optional per-phenotype pyseer covariate flag strings, parallel to `targets`, e.g. '--use-covariates 6q 7' (1-based column indices; q suffix = quantitative). Leave empty for no covariates."
    assembly_accessions: "NCBI RefSeq accessions for the reference genomes, e.g. GCF_000013305.1. These are downloaded and included in the pangenome and unitig counting, as upstream's bootstrap.sh does."
    reference_strain: "Strain name of the primary reference used for rare-variant calling, Manhattan plot coordinates and enrichment background, e.g. IAI39."
    eggnog_db_tarball: "Pre-downloaded eggNOG database tarball. Strongly recommended: leaving this unset triggers a tens-of-GB download on every run."
    deleteriousness_dir_tarball: "Sequence UNET deleteriousness predictions. Omit to skip the rare-variant gene burden test, which is out of scope for Milestone 1."
  }

  input {
    # --- Samples and phenotypes ---------------------------------------------
    File samples_tsv
    Array[File] sample_fastas
    Array[File]? sample_gffs

    Array[String] targets
    # Parallel to `targets`. Must be the same length when supplied.
    Array[String]? covariates

    # --- References -----------------------------------------------------------
    Array[String] assembly_accessions = []
    String genus
    String species
    String reference_strain

    Array[File] local_reference_fastas = []
    Array[File] local_reference_gffs   = []
    Array[File] local_reference_gbks   = []

    # config["summary_references"] / ["annotation_references"] / ["enrichment_reference"]
    String summary_references
    String annotation_references
    String enrichment_reference

    # --- Analysis parameters -------------------------------------------------
    String mlst_scheme
    String species_amr = ""

    # Optional bring-your-own lineage assignments, replacing the MLST-derived
    # ones. Mirrors upstream's `lineages_file` config escape hatch.
    File? lineages_file

    # Spurious-hit filters (config: length / min_hits / max_genes)
    Int length    = 10
    Int min_hits  = 1
    Int max_genes = 25

    # --- Optional modules ----------------------------------------------------
    File? eggnog_db_tarball
    String eggnog_taxid = "2"
    File? deleteriousness_dir_tarball
    Boolean run_amr_screen        = false
    Boolean run_whole_genome_model = true
    Boolean run_panfeed           = true
    Boolean run_heritability      = true
  }

  # The rare-variant gene burden test is gated on its one hard prerequisite
  # rather than on a separate flag someone could set inconsistently:
  # vcf2deleterious.py cannot run without deleteriousness predictions, so
  # supplying them IS enabling the burden test. Out of scope for Milestone 1.
  Boolean run_rare_variants = defined(deleteriousness_dir_tarball)

  # ==========================================================================
  # Stage 1 - reference genomes
  # ==========================================================================
  # Upstream does this in bootstrap.sh, outside Snakemake entirely. The
  # references it produces are appended to the sample lists, so they take part
  # in unitig counting and the pangenome, and one of them becomes the
  # coordinate system for the Manhattan plot.
  call setup.bootstrap_references {
    input:
      assembly_accessions     = assembly_accessions,
      genus                   = genus,
      species                 = species,
      reference_strain        = reference_strain,
      local_reference_fastas  = local_reference_fastas,
      local_reference_gffs    = local_reference_gffs,
      local_reference_gbks    = local_reference_gbks
  }

  # ==========================================================================
  # Stage 2 - gene calling
  # ==========================================================================
  # When no GFFs are supplied, upstream's aid_bootstrap.py points
  # panaroo_input.txt at out/ggcaller/GFF/{sample}.gff and the ggcaller rule
  # fills them in. The tutorial path (and therefore Tier 2) takes this branch,
  # since the 370-strain E. coli download is assemblies only.
  if (!defined(sample_gffs)) {
    call variants.ggcaller {
      input:
        sample_fastas = sample_fastas
    }
  }

  Array[File] resolved_sample_gffs = select_first([sample_gffs, ggcaller.sample_gffs])

  # ==========================================================================
  # Stage 3 - variant and feature generation
  # ==========================================================================
  call variants.unitigs {
    input:
      sample_fastas    = sample_fastas,
      reference_fastas = bootstrap_references.reference_fastas
  }

  call variants.pangenome {
    input:
      sample_gffs    = resolved_sample_gffs,
      reference_gffs = bootstrap_references.reference_gffs
  }

  if (run_panfeed) {
    call variants.panfeed_kmers {
      input:
        sample_gffs      = resolved_sample_gffs,
        reference_gffs   = bootstrap_references.reference_gffs,
        sample_fastas    = sample_fastas,
        reference_fastas = bootstrap_references.reference_fastas,
        gene_presence_absence_csv = pangenome.gene_presence_absence_csv
    }
  }

  # ==========================================================================
  # Stage 4 - population structure
  # ==========================================================================
  call popstruct.mash_sketch {
    input:
      sample_fastas    = sample_fastas,
      reference_fastas = bootstrap_references.reference_fastas
  }

  call popstruct.distance {
    input:
      sketches = mash_sketch.sketches
  }

  # Skipped entirely when lineages are supplied - MLST is slow and, on very
  # small or very clonal genome sets, can fail to type at all
  if (!defined(lineages_file)) {
    call popstruct.lineage_st {
      input:
        sample_fastas    = sample_fastas,
        reference_fastas = bootstrap_references.reference_fastas,
        mlst_scheme      = mlst_scheme
    }
  }

  File resolved_lineages = select_first([lineages_file, lineage_st.lineages])

  call popstruct.tree {
    input:
      core_genome_aln = pangenome.core_genome_aln
  }

  call popstruct.variant_sites {
    input:
      core_genome_aln = pangenome.core_genome_aln
  }

  call popstruct.aln2vcf {
    input:
      core_genome_aln = pangenome.core_genome_aln
  }

  call popstruct.similarity {
    input:
      core_genome_aln        = pangenome.core_genome_aln,
      core_genome_vcf        = aln2vcf.core_genome_vcf,
      core_genome_vcf_index  = aln2vcf.core_genome_vcf_index
  }

  # ==========================================================================
  # Stage 5 - reference-based variant calling
  # ==========================================================================
  # snippy is scattered per sample; bcftools then merges the per-sample VCFs
  # and splits them into the common and (optionally) rare sets.
  scatter (fasta in sample_fastas) {
    String snippy_samplename = basename(fasta, ".fasta")
    call rare.get_snps {
      input:
        samplename     = snippy_samplename,
        assembly_fasta = fasta,
        reference_gbk  = bootstrap_references.snps_reference_gbk
    }
  }

  call rare.prepare_vcf_variants {
    input:
      snps_vcfs        = get_snps.snps_vcf,
      snps_vcf_indexes = get_snps.snps_vcf_index,
      deleteriousness_dir_tarball = deleteriousness_dir_tarball
  }

  if (run_rare_variants) {
    call setup.prepare_regions {
      input:
        reference_gff = bootstrap_references.snps_reference_gff
    }
  }

  # ==========================================================================
  # Stage 6 - pangenome annotation (shared across all phenotypes)
  # ==========================================================================
  call annotation.sample_whole_pangenome {
    input:
      gene_presence_absence_csv = pangenome.gene_presence_absence_csv,
      gene_data                 = pangenome.gene_data,
      annotation_references     = annotation_references
  }

  if (!defined(eggnog_db_tarball)) {
    call annotation.download_eggnog {
      input:
        eggnog_taxid = eggnog_taxid
    }
  }

  call annotation.annotate_pangenome {
    input:
      pangenome_faa     = sample_whole_pangenome.pangenome_faa,
      eggnog_db_tarball = select_first([eggnog_db_tarball, download_eggnog.eggnog_db_tarball])
  }

  call annotation.annotate_reference {
    input:
      gene_presence_absence_csv = pangenome.gene_presence_absence_csv,
      annotations               = annotate_pangenome.annotations,
      enrichment_reference      = enrichment_reference
  }

  call setup.download_obo {}

  # ==========================================================================
  # Stage 7 - optional AMR/virulence side branch
  # ==========================================================================
  if (run_amr_screen) {
    call amr.find_amr_vag {
      input:
        sample_fastas = sample_fastas,
        species       = species_amr
    }
  }

  # ==========================================================================
  # Stage 8 - per-phenotype associations
  # ==========================================================================
  scatter (i in range(length(targets))) {

    String target = targets[i]
    # Covariate flag string for this phenotype, or "" when none were supplied
    String target_covariates = if defined(covariates)
                               then select_first([covariates])[i]
                               else ""

    call pyseer.prepare_pyseer {
      input:
        phenotype        = target,
        sample_fastas    = sample_fastas,
        reference_fastas = bootstrap_references.reference_fastas,
        phenotypes_tsv   = samples_tsv,
        all_similarities = similarity.similarities,
        all_distances    = distance.distances,
        all_lineages     = resolved_lineages
    }

    # --- Unitigs, gene presence/absence, structural variants ----------------
    call pyseer.run_pyseer {
      input:
        phenotype = target,
        unitigs   = unitigs.unitigs,
        gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab,
        struct_presence_absence    = pangenome.struct_presence_absence,
        phenotypes  = prepare_pyseer.phenotypes,
        similarity  = prepare_pyseer.similarity,
        distances   = prepare_pyseer.distances,
        lineages    = prepare_pyseer.lineages,
        covariates  = target_covariates
    }

    # --- Common SNPs, and rare variants when enabled -------------------------
    call pyseer.run_pyseer_vcf {
      input:
        phenotype         = target,
        common_snps       = prepare_vcf_variants.common_snps,
        common_snps_index = prepare_vcf_variants.common_snps_index,
        rare_snps         = prepare_vcf_variants.rare_snps,
        rare_snps_index   = prepare_vcf_variants.rare_snps_index,
        regions           = prepare_regions.regions,
        phenotypes        = prepare_pyseer.phenotypes,
        similarity        = prepare_pyseer.similarity,
        unitigs_patterns  = run_pyseer.unitigs_patterns,
        covariates        = target_covariates
    }

    # --- Gene-cluster k-mers -------------------------------------------------
    if (run_panfeed) {
      call pyseer.run_panfeed as panfeed_association {
        input:
          phenotype        = target,
          panfeed_patterns = select_first([panfeed_kmers.panfeed_patterns]),
          phenotypes       = prepare_pyseer.phenotypes,
          similarity       = prepare_pyseer.similarity,
          unitigs_patterns = run_pyseer.unitigs_patterns,
          covariates       = target_covariates
      }

      call summaries.annotate_panfeed_small {
        input:
          phenotype             = target,
          panfeed_associations  = panfeed_association.panfeed_results,
          panfeed_conversion    = select_first([panfeed_kmers.panfeed_conversion]),
          panfeed_patterns      = select_first([panfeed_kmers.panfeed_patterns]),
          unitigs_patterns      = run_pyseer.unitigs_patterns
      }

      call summaries.panfeed_downstream {
        input:
          phenotype            = target,
          phenotypes           = prepare_pyseer.phenotypes,
          panfeed_associations = panfeed_association.panfeed_results,
          panfeed_conversion   = select_first([panfeed_kmers.panfeed_conversion]),
          unitigs_patterns     = run_pyseer.unitigs_patterns,
          gene_presence_absence_csv = pangenome.gene_presence_absence_csv,
          sample_gffs      = resolved_sample_gffs,
          reference_gffs   = bootstrap_references.reference_gffs,
          sample_fastas    = sample_fastas,
          reference_fastas = bootstrap_references.reference_fastas
      }

      call summaries.map_summary_panfeed {
        input:
          phenotype     = target,
          panfeed_kmers = annotate_panfeed_small.panfeed_kmers,
          phenotypes    = prepare_pyseer.phenotypes,
          gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab,
          gene_presence_absence_csv  = pangenome.gene_presence_absence_csv,
          reference_gffs      = bootstrap_references.reference_gffs,
          summary_references  = summary_references,
          min_hits            = min_hits,
          max_genes           = max_genes
      }

      call summaries.annotate_summary as annotate_panfeed_summary {
        input:
          phenotype   = target,
          label       = "panfeed",
          summary     = map_summary_panfeed.summary,
          annotations = annotate_pangenome.annotations
      }

      call enrichment.run_enrich as enrich_panfeed {
        input:
          phenotype           = target,
          label               = "panfeed",
          annotated_summary   = annotate_panfeed_summary.annotated_summary,
          annotated_reference = annotate_reference.annotated_reference,
          go_obo              = download_obo.go_obo
      }

      call enrichment.run_enrichment_plots as enrich_plots_panfeed {
        input:
          phenotype = target,
          label     = "panfeed",
          cog       = enrich_panfeed.cog,
          go        = enrich_panfeed.go,
          kegg      = enrich_panfeed.kegg
      }
    }

    # --- Heritability --------------------------------------------------------
    if (run_heritability) {
      call heritability.lineages2covariance {
        input:
          phenotype = target,
          lineages  = prepare_pyseer.lineages
      }

      call heritability.run_heritability as estimate_heritability {
        input:
          phenotype           = target,
          phenotypes          = prepare_pyseer.phenotypes,
          similarity          = prepare_pyseer.similarity,
          lineages_covariance = lineages2covariance.lineages_covariance
      }

      call heritability.combine_heritability {
        input:
          phenotype              = target,
          heritability           = estimate_heritability.heritability,
          heritability_lineages  = estimate_heritability.heritability_lineages,
          heritability_ci        = estimate_heritability.heritability_ci
      }
    }

    # --- Map significant unitigs back onto genomes ---------------------------
    call mapback.run_map_back {
      input:
        phenotype         = target,
        filtered_variants = run_pyseer.unitigs_filtered,
        sample_fastas     = sample_fastas,
        reference_fastas  = bootstrap_references.reference_fastas,
        sample_gffs       = resolved_sample_gffs,
        reference_gffs    = bootstrap_references.reference_gffs,
        gene_presence_absence_csv = pangenome.gene_presence_absence_csv
    }

    call mapback.run_map_back_all {
      input:
        phenotype        = target,
        all_variants     = run_pyseer.unitigs_results,
        reference_fastas = bootstrap_references.reference_fastas,
        reference_gffs   = bootstrap_references.reference_gffs
    }

    call mapback.run_manhattan_plot {
      input:
        phenotype        = target,
        mapped_all       = run_map_back_all.mapped_all,
        unitigs_patterns = run_pyseer.unitigs_patterns,
        reference_strain = enrichment_reference
    }

    # --- Summaries and annotation for the unitig and GPA hits ----------------
    call summaries.map_summary {
      input:
        phenotype         = target,
        mapped            = run_map_back.mapped,
        phenotypes        = prepare_pyseer.phenotypes,
        filtered_variants = run_pyseer.unitigs_filtered,
        gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab,
        gene_presence_absence_csv  = pangenome.gene_presence_absence_csv,
        reference_gffs     = bootstrap_references.reference_gffs,
        summary_references = summary_references,
        length             = length,
        min_hits           = min_hits,
        max_genes          = max_genes
    }

    call summaries.annotate_summary as annotate_unitigs_summary {
      input:
        phenotype   = target,
        label       = "unitigs",
        summary     = map_summary.summary,
        annotations = annotate_pangenome.annotations
    }

    call summaries.gpa_summary {
      input:
        phenotype    = target,
        gpa_filtered = run_pyseer.gpa_filtered,
        gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab,
        gene_presence_absence_csv  = pangenome.gene_presence_absence_csv,
        reference_gffs     = bootstrap_references.reference_gffs,
        summary_references = summary_references
    }

    call summaries.annotate_summary as annotate_gpa_summary {
      input:
        phenotype   = target,
        label       = "gpa",
        summary     = gpa_summary.summary,
        annotations = annotate_pangenome.annotations
    }

    # --- Rare variant summaries (only when the burden test ran) --------------
    if (run_rare_variants) {
      call summaries.rare_summary {
        input:
          phenotype     = target,
          rare_filtered = select_first([run_pyseer_vcf.rare_filtered]),
          gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab,
          gene_presence_absence_csv  = pangenome.gene_presence_absence_csv,
          reference_gffs       = bootstrap_references.reference_gffs,
          summary_references   = summary_references,
          enrichment_reference = enrichment_reference
      }

      call summaries.annotate_summary as annotate_rare_summary {
        input:
          phenotype   = target,
          label       = "rare",
          summary     = rare_summary.summary,
          annotations = annotate_pangenome.annotations
      }

      call enrichment.run_enrich as enrich_rare {
        input:
          phenotype           = target,
          label               = "rare",
          annotated_summary   = annotate_rare_summary.annotated_summary,
          annotated_reference = annotate_reference.annotated_reference,
          go_obo              = download_obo.go_obo
      }
    }

    # --- Enrichment on the main and GPA summaries ----------------------------
    call enrichment.run_enrich as enrich_unitigs {
      input:
        phenotype           = target,
        label               = "unitigs",
        annotated_summary   = annotate_unitigs_summary.annotated_summary,
        annotated_reference = annotate_reference.annotated_reference,
        go_obo              = download_obo.go_obo
    }

    call enrichment.run_enrichment_plots as enrich_plots_unitigs {
      input:
        phenotype = target,
        label     = "unitigs",
        cog       = enrich_unitigs.cog,
        go        = enrich_unitigs.go,
        kegg      = enrich_unitigs.kegg
    }

    call enrichment.run_enrich as enrich_gpa {
      input:
        phenotype           = target,
        label               = "gpa",
        annotated_summary   = annotate_gpa_summary.annotated_summary,
        annotated_reference = annotate_reference.annotated_reference,
        go_obo              = download_obo.go_obo
    }

    call enrichment.run_enrichment_plots as enrich_plots_gpa {
      input:
        phenotype = target,
        label     = "gpa",
        cog       = enrich_gpa.cog,
        go        = enrich_gpa.go,
        kegg      = enrich_gpa.kegg
    }

    # --- QQ plots ------------------------------------------------------------
    call pyseer.run_qq_plot as qq_unitigs {
      input:
        phenotype = target, variant_class = "unitigs",
        association_results = run_pyseer.unitigs_results
    }

    call pyseer.run_qq_plot as qq_gpa {
      input:
        phenotype = target, variant_class = "gpa",
        association_results = run_pyseer.gpa_results
    }

    call pyseer.run_qq_plot as qq_vcf {
      input:
        phenotype = target, variant_class = "vcf",
        association_results = run_pyseer_vcf.common_results
    }

    if (run_panfeed) {
      call pyseer.run_qq_plot as qq_panfeed {
        input:
          phenotype = target, variant_class = "panfeed",
          association_results = select_first([panfeed_association.panfeed_results])
      }
    }

    if (run_rare_variants) {
      call pyseer.run_qq_plot as qq_rare {
        input:
          phenotype = target, variant_class = "rare",
          association_results = select_first([run_pyseer_vcf.rare_results])
      }
    }

    # --- Whole-genome elastic net -------------------------------------------
    if (run_whole_genome_model) {
      call pyseer.prepare_wg {
        input:
          phenotype        = target,
          sample_fastas    = sample_fastas,
          reference_fastas = bootstrap_references.reference_fastas,
          unitigs          = unitigs.unitigs,
          phenotypes_tsv   = samples_tsv,
          all_similarities = similarity.similarities,
          all_distances    = distance.distances,
          all_lineages     = resolved_lineages,
          covariates       = target_covariates
      }

      call pyseer.run_wg {
        input:
          phenotype    = target,
          unitigs      = unitigs.unitigs,
          variants_pkl = prepare_wg.variants_pkl,
          phenotypes   = prepare_wg.phenotypes,
          distances    = prepare_wg.distances,
          lineages     = prepare_wg.lineages,
          covariates   = target_covariates
      }

      # Both elastic net models take the same downstream path, so they are
      # scattered rather than spelled out twice. This mirrors upstream's
      # {model} wildcard over ["ridge", "lasso"].
      # Indexed parallel arrays rather than tuples: WDL 1.0's Pair holds only
      # two values, and these three have to travel together.
      Array[String] wg_models      = ["ridge", "lasso"]
      Array[File] wg_variant_files = [run_wg.ridge, run_wg.lasso]
      Array[File] wg_prediction_files = [run_wg.ridge_predictions, run_wg.lasso_predictions]

      scatter (j in range(length(wg_models))) {

        String wg_model     = wg_models[j]
        File wg_variants    = wg_variant_files[j]
        File wg_predictions = wg_prediction_files[j]

        call pyseer.run_wg_metrics {
          input:
            phenotype   = target,
            model       = wg_model,
            predictions = wg_predictions
        }

        call mapback.run_map_back_wg {
          input:
            phenotype        = target,
            model            = wg_model,
            model_variants   = wg_variants,
            sample_fastas    = sample_fastas,
            reference_fastas = bootstrap_references.reference_fastas,
            sample_gffs      = resolved_sample_gffs,
            reference_gffs   = bootstrap_references.reference_gffs,
            gene_presence_absence_csv = pangenome.gene_presence_absence_csv
        }

        call summaries.map_summary_wg {
          input:
            phenotype      = target,
            model          = wg_model,
            mapped         = run_map_back_wg.mapped,
            phenotypes     = prepare_wg.phenotypes,
            model_variants = wg_variants,
            gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab,
            gene_presence_absence_csv  = pangenome.gene_presence_absence_csv,
            reference_gffs     = bootstrap_references.reference_gffs,
            summary_references = summary_references,
            length             = length,
            min_hits           = min_hits,
            max_genes          = max_genes
        }

        call summaries.annotate_summary as annotate_wg_summary {
          input:
            phenotype   = target,
            label       = wg_model,
            summary     = map_summary_wg.summary,
            annotations = annotate_pangenome.annotations
        }

        call enrichment.run_enrich as enrich_wg {
          input:
            phenotype           = target,
            label               = wg_model,
            annotated_summary   = annotate_wg_summary.annotated_summary,
            annotated_reference = annotate_reference.annotated_reference,
            go_obo              = download_obo.go_obo
        }

        call enrichment.run_enrichment_plots as enrich_plots_wg {
          input:
            phenotype = target,
            label     = wg_model,
            cog       = enrich_wg.cog,
            go        = enrich_wg.go,
            kegg      = enrich_wg.kegg
        }
      }
    }
  }

  output {
    # --- Provenance -----------------------------------------------------------
    String pyseer_version  = similarity.pyseer_version
    String panaroo_version = pangenome.panaroo_version
    String mash_version    = mash_sketch.mash_version
    Int n_references       = bootstrap_references.n_references
    Int n_unitigs          = unitigs.n_unitigs

    # --- Shared intermediates worth keeping ----------------------------------
    File gene_presence_absence_rtab = pangenome.gene_presence_absence_rtab
    File gene_presence_absence_csv  = pangenome.gene_presence_absence_csv
    File core_genome_alignment      = pangenome.core_genome_aln
    # Upstream's `variant_sites` rule is a leaf - nothing downstream reads it -
    # but the variable-sites alignment is useful on its own, so it is surfaced
    # rather than dropped.
    File core_genome_variable_sites = variant_sites.variant_alignment
    File core_genome_tree           = tree.treefile
    File kinship_matrix             = similarity.similarities
    File distance_matrix            = distance.distances
    File lineages                   = resolved_lineages
    File pangenome_annotations      = annotate_pangenome.annotations
    File? panaroo_summary_statistics = pangenome.summary_statistics

    # --- Per-phenotype association results -----------------------------------
    Array[String] significance_thresholds = run_pyseer.significance_threshold
    Array[Int] n_unitig_hits = run_pyseer.n_unitig_hits
    Array[Int] n_gpa_hits    = run_pyseer.n_gpa_hits

    Array[File] unitigs_results  = run_pyseer.unitigs_results
    Array[File] unitigs_filtered = run_pyseer.unitigs_filtered
    Array[File] gpa_results      = run_pyseer.gpa_results
    Array[File] gpa_filtered     = run_pyseer.gpa_filtered
    Array[File] struct_filtered  = run_pyseer.struct_filtered
    Array[File] common_snp_results   = run_pyseer_vcf.common_results
    Array[File] common_snp_annotated = run_pyseer_vcf.common_annotated
    Array[File?] rare_results        = run_pyseer_vcf.rare_filtered
    Array[File?] panfeed_results     = panfeed_association.panfeed_filtered

    # --- Summaries: the tables to read for the Tier 2 loci check -------------
    Array[File] annotated_summary_unitigs = annotate_unitigs_summary.annotated_summary
    Array[File] annotated_summary_gpa     = annotate_gpa_summary.annotated_summary
    Array[File?] annotated_summary_rare   = annotate_rare_summary.annotated_summary
    Array[File?] annotated_summary_panfeed = annotate_panfeed_summary.annotated_summary
    Array[Array[File]?] annotated_summary_wg = annotate_wg_summary.annotated_summary

    # --- Plots ----------------------------------------------------------------
    Array[File] manhattan_png = run_manhattan_plot.manhattan_png
    Array[File] manhattan_svg = run_manhattan_plot.manhattan_svg
    Array[File] qq_unitigs_png = qq_unitigs.qq_plot
    Array[File] qq_gpa_png     = qq_gpa.qq_plot
    Array[File] qq_vcf_png     = qq_vcf.qq_plot
    Array[File?] qq_panfeed_png = qq_panfeed.qq_plot
    Array[File?] qq_rare_png    = qq_rare.qq_plot

    Array[File] enrichment_cog_unitigs = enrich_unitigs.cog
    Array[File] enrichment_go_unitigs  = enrich_unitigs.go
    Array[File] enrichment_kegg_unitigs = enrich_unitigs.kegg
    Array[File] enrichment_kegg_plot_unitigs = enrich_plots_unitigs.kegg_png

    Array[Array[File]?] panfeed_plots = panfeed_downstream.panfeed_plots

    # --- Whole-genome model ---------------------------------------------------
    Array[Array[File]?] wg_metrics = run_wg_metrics.metrics
    Array[File?] wg_ridge = run_wg.ridge
    Array[File?] wg_lasso = run_wg.lasso

    # --- Heritability ---------------------------------------------------------
    Array[File?] heritability_all = combine_heritability.heritability_all

    # --- Optional AMR screen --------------------------------------------------
    File? amr_summary_matches   = find_amr_vag.summary_matches
    File? amr_summary_virulence = find_amr_vag.summary_virulence
  }
}
