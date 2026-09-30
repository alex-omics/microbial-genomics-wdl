# classify_kraken2_batched

Classify many samples with kraken2 and, optionally, keep or remove the reads that fall within
given clades. Loading the kraken2 database is the expensive part of a run, so this workflow
loads it once per VM instead of once per sample.

Reads can be unaligned BAMs or FASTQ, paired or single-end. The output is one unaligned BAM per
sample: every read if no taxa are given, otherwise only the reads kept (or not removed) by the
taxon filter.

## How it works

Samples are split into batches of at most `samples_per_batch`. Each batch runs on one VM that

1. extracts the database once,
2. reads it into the page cache once,
3. classifies `concurrent_samples` samples at a time with `--memory-mapping`, so the concurrent
   kraken2 processes share that one cached copy instead of each holding their own in RAM.

**Submit it once over a set of samples**, not once per sample. Batching only happens within a
run, so one run per sample loads the database once per sample again. In Terra, launch it on a
sample set with `this.samples.<column>` for each reads input.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `reads_bams` | `Array[File]?` | Unaligned BAMs. Paired or single-end is detected per sample from the read flags |
| `reads_fastq_r1` | `Array[File]?` | FASTQ, plain or gzipped. Use instead of `reads_bams` |
| `reads_fastq_r2` | `Array[File]?` | Mates of `reads_fastq_r1`, in the same order. Omit for single-end; a run is all paired or all single-end |
| `kraken2_db_tgz` | `File` | **Required.** Tarball with `hash.k2d`, `opts.k2d`, `taxo.k2d` |
| `sample_names` | `Array[String]?` | Positionally matched. Defaults to the filename without `.bam` / `_R1` / `.fastq.gz`. Letters, digits, `.`, `_`, `-` only; must be unique |
| `filter_taxids` | `Array[Int]` | Taxa whose clades to keep. Empty (default) keeps every read |
| `filter_taxon_names` | `Array[String]` | Taxon names, as spelled in the kraken2 report. Combined with `filter_taxids` |
| `remove_matching` | `Boolean` | Remove the listed clades instead of keeping them. Default `false` |
| `keep_unclassified` | `Boolean` | When keeping, retain unclassified reads too. Default `false` |
| `confidence` | `Float` | kraken2 `--confidence`. Default `0.0` |
| `kraken2_extra_args` | `String` | Passed through, e.g. `--minimum-hit-groups 3` |
| `samples_per_batch` | `Int` | Default `100` |
| `mem_gb` / `cpu` / `concurrent_samples` | `Int` | Default `128` / `32` / `8` |
| `disk_gb` | `Int?` | Sized from the database and reads if omitted |
| `preemptible` | `Int` | Default `0` |

Give exactly one of `reads_bams` or `reads_fastq_r1`; the workflow fails at the start otherwise.

### Taxon filtering

A listed taxon keeps its whole clade: everything classified at or beneath it. Listing a family
therefore keeps reads assigned to any species or genus within it, which matters because reads
from one organism are often assigned to a close relative.

### Choosing settings

- **`mem_gb`** must exceed the uncompressed `hash.k2d`, plus a few GB of headroom (the tarball
  is compressed, so check the extracted size). If it does not, every read pages from disk and
  most of the saving is lost. The task warns when this happens.
- **`samples_per_batch`** trades database loads against blast radius. Larger batches load the
  database fewer times; a failed batch reruns in full, and preemption does too, hence the
  `preemptible` default of 0.
- **Taxa** are matched against your database's taxonomy, which changes between releases (genera
  get split, species move). Take the ids or names from a report produced with your database. A
  taxon with no reads in a sample is noted in the task log and gives that sample an empty BAM.
- **`confidence`** at `0` lets a single k-mer hit assign a read, which leaks reads into a clade.
  For a purity filter, compare `kept_pct` at `0.1` or so against the default on a few samples.

## Outputs

All per-sample outputs are in input order.

| Output | Notes |
| ------ | ----- |
| `summary_tsv` | One row per sample (below) |
| `sample_ids` | Resolved names, in input order |
| `reports` | kraken2 report per sample |
| `bams` | Unaligned BAM per sample. For BAM input the original records and header are kept; for FASTQ input a BAM is built with the sample name as its read group |

`summary_tsv` columns:

```
sample  total_reads  classified_reads  unclassified_reads  top_species_taxid  top_species
top_species_reads  top_species_pct  kept_reads  kept_pct
```

Counts are reads for single-end samples and read pairs for paired-end. `kept_reads` and
`kept_pct` are `NA` when no taxa are given. The kraken2 reports can be passed to Krona
separately; the workflow does not produce plots or per-read assignment files.

## Running

```bash
miniwdl run workflows/classify_kraken2_batched/classify_kraken2_batched.wdl \
    reads_fastq_r1=a_R1.fastq.gz reads_fastq_r2=a_R2.fastq.gz \
    reads_fastq_r1=b_R1.fastq.gz reads_fastq_r2=b_R2.fastq.gz \
    kraken2_db_tgz=k2_db.tar.gz filter_taxids=1234
```
