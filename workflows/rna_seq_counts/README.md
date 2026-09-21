# rna_seq_counts

Align paired-end RNA-seq reads with BWA-MEM and count reads per feature — with each
sample aligned to **its own isolate's assembly** rather than one shared reference.

The reason: a single reference (PA14, PAO1) can only receive reads from genes it
carries. For a diverse panel, transcripts from the accessory genome have nowhere to map
and are silently lost. Aligning each sample to the assembly it was actually sequenced
from keeps them. The same workflow also runs the classic shared-reference way
(`match_mode = "all_to_one"`), where every sample counts the same features and nothing
below about annotation or orthologs applies.

## What you have to supply

Reads, and one assembly per isolate. That is all: the workflow annotates each assembly
with **Bakta** and builds the ortholog table with **Panaroo** itself. Anything you already
have, you can supply instead and that step is skipped:

| You supply | Skipped | Why you might |
| ---------- | ------- | ------------- |
| `reference_annotations` (Bakta GFF3s) | Bakta | You annotate with your own Bakta setup, or already did |
| `ortholog_table` (Panaroo `gene_presence_absence.csv`) | Panaroo | Reuse one table across reruns, or from `pangenome` / `dna_methylation_calling`, rather than paying for Panaroo again |

Both are optional and independent. `pangenome`'s outputs line up directly:
`bakta_gff3` -> `reference_annotations`, `gene_presence_absence` -> `ortholog_table`
(positionally matched, same order as `reference_fastas` / `reference_ids`).

## What comes out

| Output | What it is |
| ------ | ---------- |
| `sample_reference_mapping` | One row per sample: which assembly it was aligned to, which matching rule fired, and any ambiguity note. **Read this first.** |
| `count_matrices` | One featureCounts matrix per isolate (features x that isolate's samples), columns named by sample |
| `counts_long` | Every isolate, sample and feature in one long-format table. Valid however the isolates were annotated |
| `counts_matrix` | One wide ortholog-group x sample matrix — see below |
| `ortholog_table_used`, `pangenome_summary`, `annotations_used` | What the matrix was built on, whether supplied or built here |
| `aligned_bams`, `markdup_bams`, `flagstats`, `pct_mapped`, `dup_metrics` | Per-sample alignment products. `pct_mapped` is the quick check that a sample went to the right assembly |

### The cross-isolate matrix

With a shared reference, every sample counts the same features, so merging is a column
join. With per-isolate assemblies it is not: each isolate's features carry its own locus
tags, so samples from different isolates share no rows. The **ortholog table** is what
puts them on common rows. It is built here with Panaroo (or supplied), and `counts_matrix`
is then one ortholog-group x sample matrix. It is built only when several isolates are in
play; set `run_pangenome = false` to skip it and keep just the per-isolate matrices and the
long table.

- A gene an isolate does not carry is **`NA`, not `0`**. A zero count means "expressed at
  zero"; absent means "cannot be measured". Conflating them manufactures differences that
  are really gene presence/absence. (Same rule as `dna_methylation_calling`.)
- Paralogues in one ortholog group are **summed**.
- **Nothing is dropped for being unclustered.** A counted feature that is in no ortholog
  group gets its own row, `unclustered:<isolate>:<locus_tag>`, with `NA` for the other
  isolates. This matters because Panaroo's stricter clean modes discard genes seen in few
  genomes as likely annotation error — exactly the accessory transcripts this workflow
  exists to keep. On the 2-genome test panel, `strict` and `moderate` kept none of the five
  isolate-specific genes and only `sensitive` kept all of them. So the workflow defaults
  to **`panaroo_clean_mode = "sensitive"`** (methylation uses `strict`, for its own
  reasons), and the `unclustered` rows are the safety net if you supply a stricter table
  such as methylation's. A `NOTE` in the merge log says how many there were.

A **supplied** table is validated before any alignment (the run stops early if an isolate
has no column, or if you passed the `.Rtab`, which has no locus tags). Whatever the source,
the merge step refuses an isolate whose locus tags match none of its counts, which would
otherwise produce an all-`NA` column with no error. Isolate columns must equal
`reference_ids`; a table built here is guaranteed to, since Panaroo is given the names
explicitly rather than inferring them from GFF filenames.

## Matching samples to assemblies

RNA-seq comes in replicates, so a sample is usually the isolate name plus a suffix.
`match_samples_to_references` resolves each sample in this order:

1. `sample_reference_map` (optional TSV: `sample_id<TAB>reference_id`) — explicit wins.
2. **Exact name** match against `reference_ids`.
3. **Strip a replicate suffix** and match again, up to `max_strip_depth` times
   (default **1**).

It only ever matches against names that are really in `reference_ids`, and an unmatched
sample **stops the run before any alignment**, listing every failure with its closest
assembly names.

Built-in replicate suffixes:

| Form | Examples | Stripped |
| ---- | -------- | -------- |
| letter directly after a digit | `PSA_1b`, `PSA_13-2a` | `b`, `a` |
| separator + letter | `PSA_26a_b`, `Iso-C`, `Iso.c` | `_b`, `-C`, `.c` |
| separator or keyword + 1-2 digits | `Iso_2`, `Iso-3`, `Iso.1`, `Iso_rep2`, `Isorep2`, `Iso_replicate3`, `Iso_R1` | `_2`, `-3`, ... |
| keyword + letter | `Iso_replicateA`, `IsoreplicateB` | `_replicateA`, ... |

Deliberate limits, because the failure mode here is a sample silently aligned to the
*wrong existing* genome:

- **A bare digit is never a replicate.** `PSA_11` is not `PSA_1` plus replicate `1`.
  Digits only strip with a separator or keyword (`_1`, `rep1`).
- **A bare letter after a letter is not a replicate** (`Isob`), or `Ecoli` would be
  `Ecol` plus replicate `i`. Nor is `R` + digits without a separator (`IsoR2`).
- Depth defaults to 1. `PSA_13-2a` resolves to `PSA_13-2`; it will *not* fall through
  to `PSA_13` unless you raise `max_strip_depth` to 2.
- If both `PSA_1` and `PSA_1b` are assemblies, the exact match wins, and the mapping
  table notes `also_strips_to:PSA_1`, since that is the one case a replicate could be
  misfiled.

For a convention the defaults do not cover, pass `replicate_regex` (a Python regex for
the suffix; it *replaces* the built-in grammar), or list the awkward samples in
`sample_reference_map`.

Only assemblies that have samples are indexed and counted.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `read1_trimmed`, `read2_trimmed`, `sample_ids` | arrays | Per sample, positionally matched |
| `reference_fastas`, `reference_ids` | arrays | Per assembly, positionally matched |
| `reference_annotations` | `Array[File]?` | GFF3 per assembly, e.g. Bakta's. Omit and Bakta runs. Must come from the same assembly with contig names kept (Bakta: `--keep-contig-headers`), which is checked before alignment. They should carry Bakta's `##FASTA` block (its default) so Panaroo can read them |
| `bakta_db`, `proteins`, `genus`, `species` | | Passed to Bakta when it runs (`genus`/`species` default to *Pseudomonas aeruginosa*; set them for other organisms). Ignored when `reference_annotations` is supplied |
| `match_mode` | `String` | `match` (default) or `all_to_one` for a single shared reference |
| `sample_reference_map`, `max_strip_depth`, `replicate_regex` | | See above |
| `ortholog_table` | `File?` | Panaroo `gene_presence_absence.csv`. Omit and Panaroo runs. Supplying it skips Panaroo |
| `run_pangenome`, `panaroo_clean_mode`, `merge_paralogs` | | Panaroo settings. `panaroo_clean_mode` defaults to `sensitive` (see above) |
| `run_featurecounts` | `Boolean` | `false` stops after alignment; no annotation needed. Default `true` |
| `strandness` | `String` | `2` (reverse-stranded, dUTP) by default |
| `ignore_duplicates` | `Boolean` | `false` by default — see *Duplicates* |
| `feature_type`, `attribute_type` | `String` | `CDS` / `locus_tag` for Bakta GFF3 |

The alignment, duplicate and featureCounts parameters carry over from the original
standalone BWA-MEM workflow; see the workflow's `parameter_meta`.

## Changes from the standalone BWA-MEM workflow

Verified by running the real images, not by reading:

- **Image pins that could not run.** `staphb/bwa:0.7.19` contains bwa only — the task
  pipes into `samtools`, which is not there. `staphb/picard:3.1.0` and
  `biocontainers/subread:*` do not exist on Docker Hub (BioContainers lives on quay.io).
  Now: a bwa 0.7.17 + samtools 1.16.1 image, Picard 3.1.0 and subread 2.0.6 from
  quay.io, all digest-pinned.
- **`-p` without `--countReadPairs` counts each mate separately** (subread >= 2.0.2),
  doubling every fragment. Now `--countReadPairs` by default (`count_read_pairs`).
- featureCounts columns were labelled with BAM paths; they are now sample names, with a
  check that no reordering swapped samples.
- Bakta's trailing `##FASTA` block is stripped from the GFF before counting.
- `ignore_duplicates` now defaults to **false** (see *Duplicates*).
- Default compute reduced (bacterial genomes are small): 8 CPU / 16 GB align.

## Duplicates

Picard `MarkDuplicates` flags read pairs with identical start and end positions as
presumed PCR duplicates. In DNA-seq that is a fair assumption, because fragments start
almost anywhere in a large genome. In RNA-seq it is not: a highly expressed transcript is
short and heavily sampled, so many independent molecules genuinely share coordinates.
Without UMIs there is no way to tell PCR copies from those, and dropping them removes the
most reads from the genes with the most signal, compressing the top of the dynamic range.

On the synthetic test data — every fragment an independent molecule, no PCR at all —
ignoring flagged duplicates still discarded 19% of the most expressed gene and 0-9% of
the weaker ones. (Illustrative rather than a real-data estimate: the fragments are spread
uniformly over short genes.)

So duplicates are still **marked** (the `PERCENT_DUPLICATION` figure is a useful library
quality check, and a very high value on a low-input library points to over-amplification)
but not **excluded from counting**. Set `ignore_duplicates = true` to restore that.

## Other defaults kept from the original

- `mapq_min = 20` drops multi-mapping reads, including reads from duplicated genes
  (rRNA operons, paralogues).
- `feature_type = CDS` counts protein-coding genes only; sRNAs, tRNAs and rRNA are not
  counted.

## Testing

`tests/run_unit_tests.sh` runs the matcher against the real 81-sample demux sheet and 16
other replicate spellings, then the whole workflow on two synthetic isolates with planted
per-gene fragment counts, asserting every count exactly (both strands, paralogue summing,
`NA` vs `0`). It runs the workflow both with a supplied ortholog table and with Panaroo
building one, and checks the `unclustered` safety net against an emulated strict table. It
also pins the `--countReadPairs` behaviour and the FASTA/GFF contig check. Regenerate the
fixtures with `tests/fixtures/rna_seq_counts/make_fixtures.py` (genes are valid ORFs, since
Panaroo discards anything else).

The Bakta branch is not in the unit suite (the database is too heavy). It was run once by
hand with the image's light database on the same synthetic genomes: Bakta recovered the
planted genes at the right coordinates with matching counts, plus one spurious ORF called
in random sequence, which correctly counted 0.

Not yet run on Terra.
