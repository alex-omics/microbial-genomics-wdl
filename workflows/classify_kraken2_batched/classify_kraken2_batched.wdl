version 1.0

import "../../tasks/kraken2.wdl" as kraken2_task
import "../../tasks/utils.wdl" as utils

workflow classify_kraken2_batched {

    meta {
        description: "Classify many unaligned BAMs with kraken2 without loading the database once per sample. Samples are split into a few batches; each batch runs on one VM that extracts the database once and classifies several samples at once against a shared, memory-mapped copy. Optionally keeps or removes the reads that fall within given clades. Submit as one run over a set of samples, not one run per sample."
        author: "Alex Arvanitis"
    }

    input {
        Array[File]+    reads_bams
        Array[String]?  sample_names
        File            kraken2_db_tgz

        # Filtering. Empty filter_taxids means classify only.
        Array[Int]      filter_taxids       = []
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
        reads_bams:         "Unaligned BAMs, paired or single-end"
        sample_names:       "Optional labels, positionally matched to reads_bams. If omitted, derived from each filename."
        kraken2_db_tgz:     "Kraken2 database tarball (hash.k2d, opts.k2d, taxo.k2d)"
        filter_taxids:      "Taxids whose clades are kept or removed. Check them against your database's own report: taxonomy changes between releases."
        remove_matching:    "Drop the reads within filter_taxids instead of keeping them (default = false)"
        keep_unclassified:  "When keeping, also retain unclassified reads (default = false)"
        confidence:         "kraken2 --confidence (default = 0.0)"
        samples_per_batch:  "Upper bound on samples per VM. The database is loaded once per batch, so larger batches cost less, but a failed batch reruns entirely (default = 100)"
        mem_gb:             "VM memory in GB; must exceed the uncompressed hash.k2d (default = 128)"
        concurrent_samples: "Samples classified at once per VM; each gets cpu / concurrent_samples threads (default = 8)"
    }

    Int n = length(reads_bams)

    scatter (f in reads_bams) {
        String derived_name = sub(basename(f), "\\.bam$", "")
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
                File   batch_bam  = reads_bams[i]
            }
        }

        call kraken2_task.kraken2_batch {
            input:
                names             = select_all(batch_name),
                bams              = select_all(batch_bam),
                db_tgz            = kraken2_db_tgz,
                filter_taxids     = filter_taxids,
                remove_matching   = remove_matching,
                keep_unclassified = keep_unclassified,
                confidence        = confidence,
                extra_args        = kraken2_extra_args,
                cpu               = cpu,
                mem_gb            = mem_gb,
                concurrent        = concurrent_samples,
                disk_gb           = disk_gb,
                preemptible       = preemptible,
                docker            = docker
        }
    }

    call utils.concat_tables {
        input:
            tables   = kraken2_batch.batch_summary,
            basename = summary_basename
    }

    output {
        File          summary_tsv   = concat_tables.merged
        Array[String] sample_ids    = validate_panel.sample_ids
        Array[File]   reports       = flatten(kraken2_batch.reports)
        Array[File]   filtered_bams = flatten(kraken2_batch.filtered_bams)
    }
}
