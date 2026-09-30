# classify_kraken2_batched

Classify many unaligned BAMs with kraken2 and, optionally, keep or remove the reads that fall
within given clades. Loading the kraken2 database is the expensive part of a run, so this
workflow loads it once per VM instead of once per sample.

## How it works

Samples are split into batches of at most `samples_per_batch`. Each batch runs on one VM that

1. extracts the database once,
2. reads it into the page cache once,
3. classifies `concurrent_samples` samples at a time with `--memory-mapping`, so the concurrent
   kraken2 processes share that one cached copy instead of each holding their own in RAM.

Paired-end BAMs are detected from the read flags; anything else is treated as single-end.

**Submit it once over a set of samples**, not once per sample. Batching only happens within a
run, so one run per sample loads the database once per sample again. In Terra, launch it on a
sample set with `this.samples.<bam column>` as `reads_bams`.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `reads_bams` | `Array[File]+` | **Required.** Unaligned BAMs |
| `kraken2_db_tgz` | `File` | **Required.** Tarball with `hash.k2d`, `opts.k2d`, `taxo.k2d` |
| `sample_names` | `Array[String]?` | Positionally matched. Defaults to the BAM filename without `.bam`. Letters, digits, `.`, `_`, `-` only; must be unique |
| `filter_taxids` | `Array[Int]` | Clades to keep or remove. Empty (default) classifies only |
| `remove_matching` | `Boolean` | Remove the listed clades instead of keeping them. Default `false` |
| `keep_unclassified` | `Boolean` | When keeping, retain unclassified reads too. Default `false` |
| `confidence` | `Float` | kraken2 `--confidence`. Default `0.0` |
| `kraken2_extra_args` | `String` | Passed through, e.g. `--minimum-hit-groups 3` |
| `samples_per_batch` | `Int` | Default `100` |
| `mem_gb` / `cpu` / `concurrent_samples` | `Int` | Default `128` / `32` / `8` |
| `disk_gb` | `Int?` | Sized from the database and BAMs if omitted |
| `preemptible` | `Int` | Default `0` |

### Choosing settings

- **`mem_gb`** must exceed the uncompressed `hash.k2d`, plus a few GB of headroom (the tarball is compressed, so check the extracted size). If it does
  not, every read pages from disk and most of the saving is lost. The task warns when this
  happens.
- **`samples_per_batch`** trades database loads against blast radius. Larger batches load the
  database fewer times; a failed batch reruns in full, and preemption does too, hence the
  `preemptible` default of 0.
- **`filter_taxids`** are matched against your database's taxonomy, which changes between
  releases (genera get split, species move). Take the ids from a report produced with your
  database. A taxid with no reads in a sample is noted in the task log and yields an empty
  filtered BAM for that sample.
- **`confidence`** at `0` lets a single k-mer hit assign a read, which leaks reads into a clade.
  For a purity filter, compare `kept_pct` at `0.1` or so against the default on a few samples.

## Outputs

All per-sample outputs are in input order.

| Output | Notes |
| ------ | ----- |
| `summary_tsv` | One row per sample (below) |
| `sample_ids` | Resolved names, in input order |
| `reports` | kraken2 report per sample |
| `filtered_bams` | Filtered unaligned BAM per sample, with the original header and all records of each kept read. Empty when `filter_taxids` is empty |

`summary_tsv` columns:

```
sample  total_reads  classified_reads  unclassified_reads  top_species_taxid  top_species
top_species_reads  top_species_pct  kept_reads  kept_pct
```

The kraken2 reports can be passed to Krona separately; the workflow does not produce plots or per-read assignment files.

Counts are reads for single-end samples and read pairs for paired-end. `kept_reads` and
`kept_pct` are `NA` when not filtering.

## Running

```bash
miniwdl run workflows/classify_kraken2_batched/classify_kraken2_batched.wdl \
    reads_bams=a.bam reads_bams=b.bam kraken2_db_tgz=k2_db.tar.gz filter_taxids=1234
```
