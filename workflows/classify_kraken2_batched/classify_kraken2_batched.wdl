version 1.0

import "../../tasks/kraken2.wdl" as kraken2_task
import "../../tasks/utils.wdl" as utils

workflow classify_kraken2_batched {

    meta {
        description: "Classify many samples (unaligned BAMs, or single- or paired-end FASTQ) with kraken2 without loading the database once per sample. Samples are split into a few batches; each batch runs on one VM that extracts the database once and classifies several samples at once against a shared, memory-mapped copy. Emits one unaligned BAM per sample holding all of its reads, or, if taxa are given, only the reads within (or outside) those clades. Submit as one run over a set of samples, not one run per sample."
        author: "Alex Arvanitis"
    }

    input {
        # Reads: either BAMs, or FASTQ (read 2 only if paired)
        Array[File]?    reads_bams
        Array[File]?    reads_fastq_r1
        Array[File]?    reads_fastq_r2
        Array[String]?  sample_names

        File            kraken2_db_tgz

        # Filtering. With no taxa, every read is kept.
        Array[Int]      filter_taxids       = []
        Array[String]   filter_taxon_names  = []
        Boolean         remove_matching     = false
        Boolean         keep_unclassified   = false

        Float           confidence          = 0.0
        String          kraken2_extra_args  = ""

        # Batching and compute
        Int             samples_per_batch   = 100
        Int             cpu                 = 32
        Int             mem_gb              = 128
        Int             concurrent_samples  = 8
        Int?            disk_gb
        Int             preemptible         = 0
        String          docker              = "quay.io/broadinstitute/viral-classify:2.5.21.0"

        String          summary_basename    = "kraken2_summary"
    }

    parameter_meta {
        reads_bams:         "Unaligned BAMs, paired or single-end (detected per sample). Use this or reads_fastq_r1."
        reads_fastq_r1:     "FASTQ files (plain or gzipped). Use this or reads_bams."
        reads_fastq_r2:     "Mate files for reads_fastq_r1, positionally matched. Omit for single-end; a run is all paired or all single-end."
        sample_names:       "Optional labels, positionally matched to the reads. If omitted, derived from each filename."
        kraken2_db_tgz:     "Kraken2 database tarball (hash.k2d, opts.k2d, taxo.k2d)"
        filter_taxids:      "Taxids whose clades (the taxon and everything beneath it) are kept. Check them against your database's own report: taxonomy changes between releases."
        filter_taxon_names: "Taxon names, matched exactly against the names in each sample's kraken2 report. Combined with filter_taxids."
        remove_matching:    "Drop the reads within the listed taxa instead of keeping them (default = false)"
        keep_unclassified:  "When keeping, also retain unclassified reads (default = false)"
        confidence:         "kraken2 --confidence (default = 0.0)"
        samples_per_batch:  "Upper bound on samples per VM. The database is loaded once per batch, so larger batches cost less, but a failed batch reruns entirely (default = 100)"
        mem_gb:             "VM memory in GB; must exceed the uncompressed hash.k2d (default = 128)"
        concurrent_samples: "Samples classified at once per VM; each gets cpu / concurrent_samples threads (default = 8)"
    }

    Array[File] primary = if defined(reads_bams) then select_first([reads_bams, []]) else select_first([reads_fastq_r1, []])
    Array[File] r2_all  = select_first([reads_fastq_r2, []])
    Boolean     has_r2  = defined(reads_fastq_r2)

    call kraken2_task.check_read_inputs {
        input:
            n_bams     = length(select_first([reads_bams, []])),
            n_fastq_r1 = length(select_first([reads_fastq_r1, []])),
            n_fastq_r2 = length(r2_all)
    }
    Int n = check_read_inputs.n_samples

    scatter (f in primary) {
        String derived_name = sub(basename(f), "([._]R?1)?(_001)?\\.(bam|fastq|fq)(\\.gz)?$", "")
    }

    call utils.validate_panel {
        input:
            derived_names    = derived_name,
            n_primary        = n,
            companion_counts = [n],
            sample_names     = sample_names
    }

    Int n_batches   = (n + samples_per_batch - 1) / samples_per_batch
    Int batch_size  = (n + n_batches - 1) / n_batches

    # Contiguous batches keep every output in input order. Indexing by position
    # costs one expression per sample, not one per sample per batch.
    scatter (b in range(n_batches)) {
        scatter (j in range(batch_size)) {
            Int i = b * batch_size + j
            if (i < n) {
                String batch_name = validate_panel.sample_ids[i]
                File   batch_r1   = primary[i]
            }
            if (i < n && has_r2) {
                File   batch_r2   = r2_all[i]
            }
        }

        call kraken2_task.kraken2_batch {
            input:
                names              = select_all(batch_name),
                reads              = select_all(batch_r1),
                reads_r2           = select_all(batch_r2),
                bam_input          = defined(reads_bams),
                db_tgz             = kraken2_db_tgz,
                filter_taxids      = filter_taxids,
                filter_taxon_names = filter_taxon_names,
                remove_matching    = remove_matching,
                keep_unclassified  = keep_unclassified,
                confidence         = confidence,
                extra_args         = kraken2_extra_args,
                cpu                = cpu,
                mem_gb             = mem_gb,
                concurrent         = concurrent_samples,
                disk_gb            = disk_gb,
                preemptible        = preemptible,
                docker             = docker
        }
    }

    call utils.concat_tables {
        input:
            tables   = kraken2_batch.batch_summary,
            basename = summary_basename
    }

    output {
        File          summary_tsv = concat_tables.merged
        Array[String] sample_ids  = validate_panel.sample_ids
        Array[File]   reports     = flatten(kraken2_batch.reports)
        Array[File]   bams        = flatten(kraken2_batch.bams)
    }
}
