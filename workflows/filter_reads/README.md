# filter_reads

Filter one sample's ONT reads on **length** and **mean read quality**. Takes an unaligned BAM, including a
modBAM, or FASTQ, and returns reads in the same format. A BAM is filtered record for record, so its `MM`/`ML`
modification tags come through untouched. FASTQ can also be cropped at either end. Built for Terra: one sample
per execution, no `Array` inputs.

This filters on properties of the reads themselves. It does **no taxonomic or reference-based filtering**, so it
is safe on metagenomes and on any run where off-target reads are the thing being measured. To remove or keep
reads by taxon, run `classify_kraken2_batched` on the filtered output.

## Terra usage

Run at the data-table level. Map `reads_bam` (or `reads_fastq`) to the column holding the reads, set the
thresholds, select the rows and launch. Outputs write back as `filter_*` columns, plus `filtered_reads`.

## What it does

| Step | Tool | Notes |
| ---- | ---- | ----- |
| Validate | | Rejects a call that would filter nothing, name two read inputs, crop a BAM, or set `min_length` above `max_length` |
| Statistics of the input | `seqkit stats` | Runs alongside the filter |
| Filter a BAM | `samtools view -e` | Length from `SEQ`, quality from the `qs` tag. Falls back to chopper for quality when records lack `qs` (below) |
| Filter FASTQ | `chopper` | Length, mean quality and optional `head_crop` / `tail_crop` |
| Statistics of the output | `seqkit stats` | |
| Summarise | | Retention, and a failure if too little survived |

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `sample_name` | `String` | **Required.** Letters, digits, `.`, `_`, `-` only |
| `reads_bam` | `File?` | Unaligned BAM, such as a dorado modBAM. Use this or `reads_fastq` |
| `reads_fastq` | `File?` | FASTQ, plain or gzipped |
| `min_length`, `max_length` | `Int?` | bp, both inclusive |
| `min_mean_q` | `Float?` | Minimum mean read quality |
| `head_crop`, `tail_crop` | `Int?` | **FASTQ only.** Bases removed from the start or end of every read |
| `qscore_source` | `String` | BAM only. `auto` (default), `qs_tag` or `computed` |
| `min_reads_out`, `min_bases_out` | `Int` | Fail if fewer survive. Defaults `1` and `0`, so an empty result fails |
| `cpu`, `mem_gb`, `disk_gb`, `preemptible` | | Defaults `8`, `16`, 3× the input + 20, `1` |
| `docker` | `String` | Defaults to `aarvani1/read-filter:1.0.0` |

At least one of the filters or crops must be set. No threshold has a default, so every limit is a decision you
made rather than one inherited from a tool.

## Outputs

| Output | Type | Notes |
| ------ | ---- | ----- |
| `filtered_reads` | `File` | `<sample>.filtered.bam` or `<sample>.filtered.fastq.gz`, in the input's format |
| `filter_stats_tsv` | `File` | seqkit statistics (reads, bases, length range and quartiles, N50, Q20/Q30 %, AvgQual, GC) for the `input` and `filtered` rows |
| `filter_reads_in`, `filter_reads_out` | `Int` | Read counts |
| `filter_bases_in`, `filter_bases_out` | `Int` | Base counts |
| `filter_pct_reads_retained`, `filter_pct_bases_retained` | `Float` | |
| `filter_n50_in`, `filter_n50_out` | `Int` | |
| `filter_qscore_source` | `String` | `qs_tag`, `computed` or `none`: where mean quality came from (`none` if `min_mean_q` was not set) |

`AvgQual` in the stats table is seqkit's: the Phred score of the mean per-base error rate over all bases. It is
not the average of the reads' own mean qualities, and long reads weigh more in it.

## How the filters behave

- **Length bounds are inclusive.** A 500 bp read passes `min_length=500` and `max_length=500`.
- **Mean quality is taken over error probabilities, not over the scores.** A read that is Q30 for half its length
  and Q2 for the rest averages Q16 as scores but Q5 as error rates, and is filtered as Q5.
- **A read exactly on `min_mean_q` is a coin toss.** Mean quality is a float, and ties are settled by rounding:
  chopper `-q 15` kept a 1000 bp Q15 read and dropped a 6000 bp one. Set the threshold between values, not on one.
- **Cropping happens before the length filter and after the quality filter.** The length bounds apply to the
  cropped read. Mean quality is judged on the whole read first, so cropping off a bad end does not rescue a read.
  This is how chopper 0.14.1 orders it, and is covered by a test.
- **A BAM is never cropped.** `MM`/`ML` store positions in the read, and cutting a read without rewriting them
  would leave modification calls on the wrong bases. Crop FASTQ reads, which carry no modification calls.
- **An aligned BAM is refused.** Length is read from `SEQ`, and a hard-clipped or supplementary record's `SEQ` is a
  fragment of the read. Filter before aligning.

### Where a BAM's quality comes from

dorado writes each read's mean quality to a `qs` tag. With `qscore_source=auto`:

1. If **every** record has `qs`, the filter is a single `samtools view -e 'length(seq)>=… && [qs]>=…'` pass. The
   count is over the whole file, because a record missing `qs` would fail the test and be dropped without any message.
2. Otherwise quality is **computed from the quality string** with chopper, the same arithmetic as the FASTQ path,
   and samtools then extracts exactly the records chopper kept. The tags still come through, because the BAM
   itself is what gets filtered.

`qs_tag` fails if any record lacks the tag; `computed` always uses chopper. `filter_qscore_source` records which
route ran. The two routes agree on the fixtures but are not guaranteed to match to the last decimal on real data,
because dorado computes `qs` its own way.

## Using the output

- **Methylation:** a filtered modBAM goes straight into `dna_methylation_calling`. That workflow's `align_modbam`
  has its own `min_read_length` (default 500), which is harmless after this filter but redundant.
- **Taxonomic filtering:** `classify_kraken2_batched` takes unaligned BAMs or FASTQ.

## Choosing thresholds

There is no universal floor. The right `min_mean_q` depends on the pore chemistry and the basecalling model, so
look at the `input` row of `filter_stats_tsv` (or the run's own report) before setting it, then check what a
threshold costs: `filter_pct_bases_retained` is usually the number to watch, since a filter can drop most of the
reads and few of the bases, or the reverse. Raise `min_reads_out` / `min_bases_out` to the amount your next step
needs, so a threshold that is too tight fails here rather than downstream.

## Container

`docker/read-filter` holds samtools, chopper and seqkit together; no public image has all three. The tasks pin
`aarvani1/read-filter:1.0.0` by digest. To rebuild it: `docker build --platform linux/amd64 -t aarvani1/read-filter:1.0.0 docker/read-filter`.

## Testing

`tests/run_unit_tests.sh` runs the whole workflow against eight planted reads (`tests/fixtures/filter_reads`),
whose lengths and qualities make each filter's right answer known in advance.
