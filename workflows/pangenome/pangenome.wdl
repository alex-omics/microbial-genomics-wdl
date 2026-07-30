version 1.0

import "../../tasks/bakta.wdl" as bakta_task
import "../../tasks/panaroo.wdl" as panaroo_task

workflow pangenome {

    meta {
        description: "Build a pangenome across a panel of assemblies with Bakta and Panaroo, and stop there. Deliberately produces no alignment and no tree: gene presence/absence and ortholog groups are useful on their own, and coupling them to phylogenetics forces every pangenome run to pay for a tree search that may not converge on fragmented input."
        author: "Alex Arvanitis"
    }

    input {
        Array[File]     assemblies
        Array[String]?  sample_names

        # Skip Bakta by supplying annotations directly. GFF3s must carry their
        # own sequence or be accompanied by matching FASTAs.
        Array[File]?    gff3s
        Array[File]?    fnas

        File?           proteins
        File?           bakta_db
        String          genus                = "Pseudomonas"
        String          species              = "aeruginosa"

        String          clean_mode           = "strict"
        Float           core_threshold       = 0.95
        Float           seq_id               = 0.95
        Float           family_threshold     = 0.7
        String          refind_mode          = "default"
        Boolean         merge_paralogs       = false
        Boolean         remove_invalid_genes = true

        # Alignment is off by default and should usually stay off. Turning it on
        # is what makes a pangenome run expensive, and a core alignment is only
        # worth building if something downstream actually consumes it.
        String?         alignment
        String?         aligner

        Int             panaroo_cpu          = 16
        Int             panaroo_mem_gb       = 64
        Int             panaroo_disk_gb      = 200
    }

    parameter_meta {
        assemblies:           "Assembly FASTAs making up the panel. Ignored if gff3s is supplied."
        sample_names:         "Optional labels, positionally matched to assemblies. If omitted, names are derived from each filename. These become the isolate column headers in gene_presence_absence.csv."
        gff3s:                "Pre-existing annotations, to skip Bakta. Must be positionally matched to fnas."
        fnas:                 "Assembly FASTAs matching gff3s, used to append sequence to any GFF lacking a ##FASTA block"
        proteins:             "FASTA of trusted proteins for Bakta --proteins, e.g. the PAO1 proteome. Passing one keeps gene names consistent across independently annotated isolates."
        bakta_db:             "Optional .tar.gz of the full Bakta database. Omitted, the light database in the image is used."
        merge_paralogs:       "Collapse paralogous families. Off by default; paralogues are usually biologically meaningful (default = false)"
        alignment:            "Optionally 'core' or 'pan'. Off by default — see the workflow description."
        panaroo_mem_gb:       "Panaroo memory. 64 GB suits a panel of a few dozen complete assemblies; genus-scale panels of ~1000+ genomes need substantially more (default = 64)"
    }

    Boolean annotate_here = !defined(gff3s)

    if (annotate_here) {
        scatter (i in range(length(assemblies))) {

            String resolved_name = if defined(sample_names)
                then select_first([sample_names])[i]
                else sub(basename(assemblies[i]), "\\.(fasta|fa|fna|fsa|contigs)(\\.gz)?$", "")

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
        }
    }

    call panaroo_task.panaroo {
        input:
            gff3s                = select_first([gff3s, bakta.gff3]),
            fnas                 = select_first([fnas,  bakta.fna]),
            clean_mode           = clean_mode,
            core_threshold       = core_threshold,
            seq_id               = seq_id,
            family_threshold     = family_threshold,
            refind_mode          = refind_mode,
            merge_paralogs       = merge_paralogs,
            remove_invalid_genes = remove_invalid_genes,
            alignment            = alignment,
            aligner              = aligner,
            cpu                  = panaroo_cpu,
            mem_gb               = panaroo_mem_gb,
            disk_gb              = panaroo_disk_gb
    }

    output {
        File   gene_presence_absence      = panaroo.gene_presence_absence
        File   gene_presence_absence_rtab = panaroo.gene_presence_absence_rtab
        File   summary_statistics         = panaroo.summary_statistics
        File?  pan_genome_reference       = panaroo.pan_genome_reference
        String n_core_genes               = panaroo.n_core_genes
        String n_total_genes              = panaroo.n_total_genes
        String panaroo_version            = panaroo.panaroo_version

        Array[File]? bakta_gff3           = bakta.gff3
        Array[File]? bakta_faa            = bakta.faa
        Array[File]? bakta_tsv            = bakta.annotation_tsv
    }
}
