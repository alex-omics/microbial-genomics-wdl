version 1.0

import "../../tasks/sample_matching.wdl" as matching
import "../../tasks/bwa.wdl" as bwa_task
import "../../tasks/picard.wdl" as picard_task
import "../../tasks/bakta.wdl" as bakta_task
import "../../tasks/panaroo.wdl" as panaroo_task
import "../../tasks/featurecounts.wdl" as counts_task

workflow rna_seq_counts {

    meta {
        description: "Align paired-end RNA-seq reads with BWA-MEM and count reads per feature, with each sample aligned to its own isolate's assembly rather than one shared reference. Replicates are matched to their parent assembly by name. Isolates aligned to their own assemblies carry independent locus tags, so an ortholog table (Panaroo) is built to make their counts comparable, after annotating with Bakta; either step is skipped if its output is supplied. Also works against a single shared reference, where neither is needed. Emits per-isolate count matrices, a long-format table, and one ortholog-level matrix."
        author: "Alex Arvanitis"
    }

    input {
        # Per-sample, positionally matched
        Array[File]    read1_trimmed
        Array[File]    read2_trimmed
        Array[String]  sample_ids

        # Per-reference, positionally matched. For per-isolate alignment these
        # are the isolates' assemblies; for a shared reference, supply one and
        # set match_mode = "all_to_one".
        Array[File]    reference_fastas
        Array[String]  reference_ids

        # How samples are matched to references
        String         match_mode           = "match"
        File?          sample_reference_map
        Int            max_strip_depth      = 1
        String?        replicate_regex

        # Alignment
        Int            mapq_min             = 20
        Boolean        mark_secondary       = true
        String         samtools_filter_flags = "0x904"
        String         bwa_extra_args       = ""

        # Duplicates
        Boolean        skip_markdup         = false
        Boolean        remove_duplicates    = false

        # Counting. run_featurecounts = false stops after alignment.
        Boolean        run_featurecounts    = true
        String         strandness           = "2"
        String         feature_type         = "CDS"
        String         attribute_type       = "locus_tag"
        String         annotation_format    = "GFF"
        Boolean        paired_end_counting  = true
        Boolean        count_read_pairs     = true
        Boolean        require_both_mates   = true
        Boolean        count_chimeric       = false
        Boolean        ignore_duplicates    = false
        Boolean        fraction_counting    = false
        Int            min_overlap          = 1
        String         output_prefix        = "counts"

        # Annotation. Supplied GFFs skip Bakta entirely.
        Array[File]?   reference_annotations
        File?          bakta_db
        File?          proteins
        String         genus                = "Pseudomonas"
        String         species              = "aeruginosa"

        # Ortholog table. Supplied, it is used as is; otherwise built with Panaroo.
        File?          ortholog_table
        Boolean        run_pangenome        = true
        String         panaroo_clean_mode   = "sensitive"
        Boolean        merge_paralogs       = false

        # Compute
        Int            index_cpu            = 2
        Int            index_mem_gb         = 8
        Int            align_cpu            = 8
        Int            align_mem_gb         = 16
        Int            align_disk_gb        = 100
        Int            markdup_mem_gb       = 16
        Int            count_cpu            = 8
        Int            count_mem_gb         = 16
    }

    parameter_meta {
        read1_trimmed:          "Trimmed forward reads, one per sample"
        read2_trimmed:          "Trimmed reverse reads, one per sample"
        sample_ids:             "Sample names. Replicates are expected to be the isolate name plus a suffix (see match_mode)."
        reference_fastas:       "Reference/assembly FASTAs, one per isolate"
        reference_ids:          "Isolate names, one per FASTA. Samples are matched against these, and they become the isolate column headers of the ortholog table."
        match_mode:             "'match' aligns each sample to the reference its name resolves to; 'all_to_one' aligns every sample to the single supplied reference (default = match)"
        sample_reference_map:   "Optional TSV (sample_id, reference_id) overriding name matching for the samples it lists"
        max_strip_depth:        "Replicate suffixes that may be peeled off a sample name while looking for its reference (default = 1)"
        replicate_regex:        "Regex for the replicate suffix, replacing the built-in grammar (1/2/3, a/b/c, rep2, ...)"
        run_featurecounts:      "Count reads per feature. Off, the workflow stops after alignment and needs no annotation (default = true)"
        strandness:             "0 = unstranded, 1 = stranded, 2 = reverse-stranded (default = 2)"
        count_read_pairs:       "Count fragments not reads (default = true). Without it subread >= 2.0.2 counts each mate."
        ignore_duplicates:      "Ignore duplicate-flagged reads when counting (default = false). Duplicates are still marked, and PERCENT_DUPLICATION still reported, as a QC metric."
        reference_annotations:  "Pre-computed GFF3 per reference (e.g. Bakta), positionally matched to reference_fastas. Supplying them skips Bakta. Contig names must match the assembly (Bakta: --keep-contig-headers), which is checked before alignment. To build the ortholog table they should include Bakta's ##FASTA block (its default), or the assembly must be uncompressed."
        bakta_db:               "Optional full Bakta database .tar.gz. Omitted, the light database in the image is used. Ignored when reference_annotations is supplied."
        proteins:               "Trusted proteins for Bakta --proteins, e.g. the PAO1 proteome. Ignored when reference_annotations is supplied."
        genus:                  "Genus for Bakta (default = Pseudomonas). Ignored when reference_annotations is supplied."
        species:                "Species for Bakta (default = aeruginosa). Ignored when reference_annotations is supplied."
        ortholog_table:         "Panaroo gene_presence_absence.csv, whose isolate columns equal reference_ids. Supplying it skips Panaroo. Checked before any alignment."
        run_pangenome:          "Build the ortholog table with Panaroo when several isolates are used and none is supplied. Off, you get the per-isolate matrices and the long table but no cross-isolate matrix (default = true)."
        panaroo_clean_mode:     "Panaroo error-correction: strict, moderate or sensitive (default = sensitive). strict and moderate discard genes found in few genomes as likely annotation error, which removes the accessory genes this workflow exists to keep."
        merge_paralogs:         "Panaroo --merge_paralogs (default = false). See the panaroo task for the tradeoff."
    }

    # With no annotations supplied there is nothing to length-check; comparing
    # against the FASTA count keeps that case from reading as a mismatch.
    Boolean annotations_supplied = defined(reference_annotations)
    Int annotation_count = if annotations_supplied
        then length(select_first([reference_annotations]))
        else length(reference_fastas)

    # Fail fast, and work out which references are actually used.
    call matching.match_samples_to_references as match {
        input:
            sample_ids              = sample_ids,
            reference_ids           = reference_ids,
            sample_array_lengths    = [length(read1_trimmed), length(read2_trimmed)],
            reference_array_lengths = [length(reference_fastas), annotation_count],
            mode                    = match_mode,
            sample_reference_map    = sample_reference_map,
            max_strip_depth         = max_strip_depth,
            replicate_regex         = replicate_regex,
            ortholog_table          = ortholog_table
    }

    # Per used reference (not per sample): index it, and get it an annotation.
    scatter (a in match.used_reference_idx) {
        String used_reference_name = reference_ids[a]

        # Only a SUPPLIED GFF is checked against the FASTA here. Bakta runs with
        # --keep-contig-headers hardcoded, so its own GFF cannot disagree, and
        # not depending on it lets indexing and alignment proceed while Bakta runs.
        if (annotations_supplied) {
            File annotation_for_check = select_first([reference_annotations])[a]
        }
        call bwa_task.bwa_index {
            input:
                reference_fasta = reference_fastas[a],
                reference_name  = reference_ids[a],
                annotation      = annotation_for_check,
                cpu             = index_cpu,
                mem_gb          = index_mem_gb
        }

        if (run_featurecounts) {
            if (!annotations_supplied) {
                call bakta_task.bakta {
                    input:
                        assembly    = reference_fastas[a],
                        sample_name = reference_ids[a],
                        bakta_db    = bakta_db,
                        proteins    = proteins,
                        genus       = genus,
                        species     = species,
                        strain      = reference_ids[a]
                }
            }
            # Whichever branch ran, this isolate's annotation now lives under one name.
            File resolved_gff3 = if annotations_supplied
                then select_first([reference_annotations])[a] else select_first([bakta.gff3])
            File resolved_fna  = if annotations_supplied
                then reference_fastas[a] else select_first([bakta.fna])
        }
    }

    scatter (i in range(length(sample_ids))) {

        call bwa_task.bwa_mem {
            input:
                read1                 = read1_trimmed[i],
                read2                 = read2_trimmed[i],
                sample_name           = sample_ids[i],
                reference_bundle      = bwa_index.reference_bundle[match.sample_to_used_pos[i]],
                mapq_min              = mapq_min,
                mark_secondary        = mark_secondary,
                samtools_filter_flags = samtools_filter_flags,
                extra_args            = bwa_extra_args,
                cpu                   = align_cpu,
                mem_gb                = align_mem_gb,
                disk_gb               = align_disk_gb
        }

        if (!skip_markdup) {
            call picard_task.mark_duplicates {
                input:
                    input_bam         = bwa_mem.bam,
                    input_bai         = bwa_mem.bai,
                    sample_name       = sample_ids[i],
                    remove_duplicates = remove_duplicates,
                    mem_gb            = markdup_mem_gb
            }
        }

        # Count from the duplicate-marked BAM when there is one.
        File bam_for_counting = select_first([mark_duplicates.markdup_bam, bwa_mem.bam])
    }

    if (run_featurecounts) {

        # Isolates aligned to their own assemblies carry independent locus tags, so
        # orthology, not coordinates, is what makes their counts comparable. Built
        # only when it is needed (several isolates), wanted, and not already supplied.
        if (run_pangenome && !defined(ortholog_table) && length(match.used_reference_idx) > 1) {
            call panaroo_task.panaroo {
                input:
                    gff3s          = select_all(resolved_gff3),
                    fnas           = select_all(resolved_fna),
                    sample_names   = used_reference_name,
                    clean_mode     = panaroo_clean_mode,
                    merge_paralogs = merge_paralogs
            }
        }
        File? resolved_ortholog_table = if defined(ortholog_table)
            then ortholog_table else panaroo.gene_presence_absence

        # Group each isolate's samples together and count them against that
        # isolate's own annotation.
        scatter (k in range(length(match.used_reference_idx))) {

            scatter (i in range(length(sample_ids))) {
                if (match.sample_to_used_pos[i] == k) {
                    File   group_bam    = bam_for_counting[i]
                    String group_sample = sample_ids[i]
                }
            }

            call counts_task.featurecounts {
                input:
                    input_bams         = select_all(group_bam),
                    sample_ids         = select_all(group_sample),
                    annotation         = select_first([resolved_gff3[k]]),
                    reference_name     = used_reference_name[k],
                    strandness         = strandness,
                    feature_type       = feature_type,
                    attribute_type     = attribute_type,
                    annotation_format  = annotation_format,
                    paired_end         = paired_end_counting,
                    count_read_pairs   = count_read_pairs,
                    require_both_mates = require_both_mates,
                    count_chimeric     = count_chimeric,
                    ignore_duplicates  = ignore_duplicates,
                    fraction_counting  = fraction_counting,
                    min_overlap        = min_overlap,
                    cpu                = count_cpu,
                    mem_gb             = count_mem_gb
            }
        }

        call counts_task.merge_count_matrices as merge {
            input:
                matrices       = featurecounts.count_matrix,
                reference_ids  = used_reference_name,
                ortholog_table = resolved_ortholog_table,
                basename       = output_prefix
        }
    }

    output {
        File           sample_reference_mapping = match.mapping

        Array[File]    aligned_bams             = bwa_mem.bam
        Array[File]    aligned_bais             = bwa_mem.bai
        Array[File]    flagstats                = bwa_mem.flagstat
        Array[String]  pct_mapped               = bwa_mem.pct_mapped
        Array[File?]   markdup_bams             = mark_duplicates.markdup_bam
        Array[File?]   markdup_bais             = mark_duplicates.markdup_bai
        Array[File?]   dup_metrics              = mark_duplicates.metrics_file

        Array[File?]   annotations_used         = resolved_gff3
        Array[File]?   count_matrices           = featurecounts.count_matrix
        Array[File]?   count_summaries          = featurecounts.summary
        File?          ortholog_table_used      = resolved_ortholog_table
        File?          pangenome_summary        = panaroo.summary_statistics
        File?          counts_long              = merge.counts_long
        File?          counts_matrix            = merge.counts_matrix
    }
}
