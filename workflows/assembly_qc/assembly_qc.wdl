version 1.0

import "../../tasks/quast.wdl" as quast_task
import "../../tasks/busco.wdl" as busco_task
import "../../tasks/checkm2.wdl" as checkm2_task

workflow assembly_qc {

    meta {
        description: "Assess the quality and completeness of one or more assembled microbial genomes with QUAST, BUSCO, and CheckM2, and collapse the results into a single summary table."
        author: "Alex Arvanitis"
    }

    input {
        Array[File]     assemblies
        Array[String]?  sample_names

        Boolean         run_quast       = true
        Boolean         run_busco       = true
        Boolean         run_checkm2     = true

        # BUSCO
        String          busco_lineage   = "auto"
        File?           busco_lineage_tarball

        # CheckM2 — stage the DIAMOND DB in GCS once and point at it here.
        File?           checkm2_db

        # QUAST comparative mode
        File?           reference_genome
        File?           reference_annotation
        String          genome_type     = "prokaryote"
        Int             min_contig      = 500

        String          summary_basename = "assembly_qc_summary"
    }

    parameter_meta {
        assemblies:            "Assembly FASTAs to assess. This is the whole point of the workflow: point it at a pile of assemblies gathered from anywhere and get one table back."
        sample_names:          "Optional labels, positionally matched to assemblies. If omitted, names are derived from each filename."
        busco_lineage:         "BUSCO dataset name, or 'auto' to place each assembly in the prokaryote tree independently (default = auto)"
        busco_lineage_tarball: "Optional BUSCO lineage .tar.gz for a dataset not baked into the image"
        checkm2_db:            "CheckM2 DIAMOND database (~3 GB). Omitting it makes every scattered task download its own copy."
        genome_type:           "QUAST domain hint: prokaryote, eukaryote, or fungus (default = prokaryote)"
    }

    # Deriving names from filenames keeps the common case a one-input run.
    scatter (i in range(length(assemblies))) {

        String resolved_name = if defined(sample_names)
            then select_first([sample_names])[i]
            else sub(basename(assemblies[i]), "\\.(fasta|fa|fna|fsa|contigs)(\\.gz)?$", "")

        if (run_quast) {
            call quast_task.quast {
                input:
                    assembly             = assemblies[i],
                    sample_name          = resolved_name,
                    reference_genome     = reference_genome,
                    reference_annotation = reference_annotation,
                    genome_type          = genome_type,
                    min_contig           = min_contig
            }
        }

        if (run_busco) {
            call busco_task.busco {
                input:
                    assembly        = assemblies[i],
                    sample_name     = resolved_name,
                    busco_lineage   = busco_lineage,
                    lineage_tarball = busco_lineage_tarball
            }
        }

        if (run_checkm2) {
            call checkm2_task.checkm2 {
                input:
                    assembly    = assemblies[i],
                    sample_name = resolved_name,
                    checkm2_db  = checkm2_db
            }
        }

        # One row per assembly. Every field is stringified and defaulted to NA so
        # the table stays rectangular regardless of which tools were switched off.
        Array[String] summary_row = [
            resolved_name,
            if defined(quast.n_contigs)         then "~{quast.n_contigs}"         else "NA",
            if defined(quast.total_length)      then "~{quast.total_length}"      else "NA",
            if defined(quast.largest_contig)    then "~{quast.largest_contig}"    else "NA",
            if defined(quast.n50)               then "~{quast.n50}"               else "NA",
            if defined(quast.l50)               then "~{quast.l50}"               else "NA",
            if defined(quast.gc_percent)        then "~{quast.gc_percent}"        else "NA",
            if defined(busco.complete_pct)      then "~{busco.complete_pct}"      else "NA",
            if defined(busco.single_copy_pct)   then "~{busco.single_copy_pct}"   else "NA",
            if defined(busco.duplicated_pct)    then "~{busco.duplicated_pct}"    else "NA",
            if defined(busco.fragmented_pct)    then "~{busco.fragmented_pct}"    else "NA",
            if defined(busco.missing_pct)       then "~{busco.missing_pct}"       else "NA",
            if defined(busco.n_markers)         then "~{busco.n_markers}"         else "NA",
            if defined(busco.lineage_used)      then "~{busco.lineage_used}"      else "NA",
            if defined(checkm2.completeness)    then "~{checkm2.completeness}"    else "NA",
            if defined(checkm2.contamination)   then "~{checkm2.contamination}"   else "NA",
            if defined(checkm2.model_used)      then "~{checkm2.model_used}"      else "NA",
            if defined(checkm2.coding_density)  then "~{checkm2.coding_density}"  else "NA"
        ]
    }

    Array[String] summary_header = [
        "sample", "quast_contigs", "quast_total_length", "quast_largest_contig",
        "quast_n50", "quast_l50", "quast_gc_percent",
        "busco_complete_pct", "busco_single_copy_pct", "busco_duplicated_pct",
        "busco_fragmented_pct", "busco_missing_pct", "busco_n_markers", "busco_lineage",
        "checkm2_completeness", "checkm2_contamination", "checkm2_model", "checkm2_coding_density"
    ]

    # write_tsv alone emits a temp-named file, which lands in the bucket as
    # something like "tmpzsaunsvb". Route it through a task so the deliverable
    # arrives with a name you can recognize.
    call name_summary {
        input:
            rows     = flatten([[summary_header], summary_row]),
            basename = summary_basename
    }

    output {
        File summary_tsv = name_summary.summary

        Array[File?] quast_report_txt   = quast.report_txt
        Array[File?] quast_report_tsv   = quast.report_tsv
        Array[File?] quast_report_html  = quast.report_html

        Array[File?] busco_summary_txt  = busco.summary_txt
        Array[File?] busco_full_table   = busco.full_table_tsv

        Array[File?] checkm2_report_tsv = checkm2.report_tsv

        Array[String?] busco_lineages_used = busco.lineage_used
    }
}

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
