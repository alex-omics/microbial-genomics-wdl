# rna_seq_counts

Align paired-end RNA-seq reads with BWA-MEM and count reads per feature, with each sample
aligned to **its own isolate's assembly** rather than a shared reference.

A single reference can only receive reads from genes it carries, so accessory-genome
transcripts have nowhere to map, and reads from divergent strains map with a bias against
them. Aligning each sample to the assembly it came from avoids both. The workflow also runs
against a single shared reference (`match_mode = "all_to_one"`), in which case every sample
counts the same features and the annotation and ortholog steps below do not apply.

**Status: in development.**

## What you supply

Reads and one assembly per isolate. The workflow annotates each assembly with **Bakta** and
builds an ortholog table with **Panaroo**. Either step is skipped if you supply its output:

| Supply | Skips | Notes |
| ------ | ----- | ----- |
| `reference_annotations` (GFF3 per assembly) | Bakta | Positionally matched to `reference_fastas` |
| `ortholog_table` (Panaroo `gene_presence_absence.csv`) | Panaroo | Isolate columns must equal `reference_ids` |

## Outputs

| Output | Contents |
| ------ | -------- |
| `sample_reference_mapping` | Per sample: the assembly it was aligned to, the matching rule that fired, and any ambiguity note. Check this first |
| `count_matrices` | One featureCounts matrix per isolate (features x that isolate's samples), columns named by sample |
| `counts_long` | Every isolate, sample and feature in one long-format table |
| `counts_matrix` | One ortholog-group x sample matrix (see below) |
| `ortholog_table_used`, `pangenome_summary`, `annotations_used` | The inputs the matrix was built from |
| `aligned_bams`, `markdup_bams`, `flagstats`, `pct_mapped`, `dup_metrics` | Per-sample alignment outputs. `pct_mapped` is a quick check that a sample went to the right assembly |

### The cross-isolate matrix

Isolates aligned to their own assemblies carry independent locus tags, so their samples share
no rows. The ortholog table puts them on common rows. It is built only when several isolates
are used; set `run_pangenome = false` to skip it and keep the per-isolate matrices and the
long table.

- A gene an isolate lacks is **`NA`, not `0`**: zero means expressed at zero, `NA` means it
  cannot be measured.
- Paralogues in one ortholog group are **summed**.
- A counted feature in no ortholog group gets its own row, `unclustered:<isolate>:<locus_tag>`,
  with `NA` for other isolates, so no expression is dropped. Panaroo's `strict` and `moderate`
  clean modes discard genes seen in few genomes, which are the accessory genes this workflow
  is meant to keep, so `panaroo_clean_mode` defaults to `sensitive`.

A supplied table is validated before any alignment: the run stops if an isolate has no
column, or if the `.Rtab` (which has no locus tags) was supplied instead of the `.csv`. The
merge step also refuses an isolate whose locus tags match none of its counted features. A
table built here names its isolate columns from `reference_ids`.

## Matching samples to assemblies

RNA-seq is usually replicated, so a sample name is typically the isolate name plus a suffix.
`match_samples_to_references` resolves each sample in order:

1. `sample_reference_map` (optional TSV: `sample_id<TAB>reference_id`): explicit wins.
2. **Exact name** match against `reference_ids`.
3. **Strip a replicate suffix** and match again, up to `max_strip_depth` times (default 1).

Only names present in `reference_ids` are ever matched, and an unmatched sample **stops the
run before any alignment**, listing every failure with its closest assembly names.

Built-in replicate suffixes:

| Form | Examples | Stripped |
| ---- | -------- | -------- |
| letter directly after a digit | `ISO_1b`, `ISO_13-2a` | `b`, `a` |
| separator + letter | `ISO_26a_b`, `Iso-C`, `Iso.c` | `_b`, `-C`, `.c` |
| separator or keyword + 1-2 digits | `Iso_2`, `Iso-3`, `Iso_rep2`, `Isorep2`, `Iso_replicate3`, `Iso_R1` | `_2`, `-3`, ... |
| keyword + letter | `Iso_replicateA`, `IsoreplicateB` | `_replicateA`, ... |

The limits are deliberate, because the failure to avoid is a sample silently aligned to the
wrong genome:

- **A bare digit is never a replicate.** `ISO_11` is not `ISO_1` plus replicate `1`.
- **A bare letter after a letter is not a replicate** (`Isob`), or `Ecoli` would be `Ecol`
  plus replicate `i`. Nor is `R` + digits without a separator (`IsoR2`).
- `max_strip_depth` defaults to 1: `ISO_13-2a` resolves to `ISO_13-2`, and only reaches `ISO_13`
  if you set it to 2.
- If both `ISO_1` and `ISO_1b` are assemblies, the exact match wins and the mapping table
  notes `also_strips_to:ISO_1`.

For other conventions, pass `replicate_regex` (a Python regex for the suffix, replacing the
built-in grammar) or list the samples in `sample_reference_map`. Only assemblies with samples
are indexed and counted.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `read1_trimmed`, `read2_trimmed`, `sample_ids` | arrays | Per sample, positionally matched |
| `reference_fastas`, `reference_ids` | arrays | Per assembly, positionally matched |
| `reference_annotations` | `Array[File]?` | GFF3 per assembly. Omitted, Bakta runs. Contig names must match the assembly (Bakta: `--keep-contig-headers`), which is checked before alignment. To build the ortholog table they should include Bakta's `##FASTA` block |
| `bakta_db`, `proteins`, `genus`, `species` | | Passed to Bakta when it runs. `genus`/`species` default to *Pseudomonas aeruginosa*; change them for other organisms |
| `ortholog_table` | `File?` | Omitted, Panaroo runs |
| `run_pangenome`, `panaroo_clean_mode`, `merge_paralogs` | | Panaroo settings. `panaroo_clean_mode` defaults to `sensitive` |
| `match_mode` | `String` | `match` (default) or `all_to_one` |
| `sample_reference_map`, `max_strip_depth`, `replicate_regex` | | See above |
| `run_featurecounts` | `Boolean` | `false` stops after alignment and needs no annotation. Default `true` |
| `strandness` | `String` | `2` (reverse-stranded, dUTP) by default |
| `feature_type`, `attribute_type` | `String` | `CDS` / `locus_tag`, for Bakta GFF3 |
| `count_read_pairs` | `Boolean` | Count fragments rather than reads (default `true`). Without it, featureCounts >= 2.0.2 counts each mate separately |
| `ignore_duplicates` | `Boolean` | Default `false`; see below |

The remaining alignment, duplicate and counting parameters are described in the workflow's
`parameter_meta`.

## Duplicates

Picard `MarkDuplicates` flags read pairs with identical start and end positions as PCR
duplicates. That is a fair assumption for DNA-seq, but in RNA-seq a highly expressed
transcript is short and heavily sampled, so many independent molecules share coordinates.
Without UMIs they cannot be told apart from PCR copies, and excluding them removes the most
reads from the most highly expressed genes.

Duplicates are therefore **marked** (`PERCENT_DUPLICATION` is still a useful library-quality
metric) but not **excluded from counting**. Set `ignore_duplicates = true` to exclude them.

## Other defaults

- `mapq_min = 20` drops multi-mapping reads, including reads from duplicated genes such as
  rRNA operons and paralogues.
- `feature_type = CDS` counts protein-coding genes only; sRNAs, tRNAs and rRNA are not counted.

## Testing

`tests/run_unit_tests.sh` runs the sample matcher against a real sample sheet and a range of
replicate spellings, then the whole workflow on two synthetic isolates with planted per-gene
fragment counts, asserting every count exactly. It covers both strands, paralogue summing,
`NA` vs `0`, a supplied and a Panaroo-built ortholog table, and the `unclustered` rows.
Regenerate the fixtures with `tests/fixtures/rna_seq_counts/make_fixtures.py`. The Bakta step
is not part of the automated tests because its database is large.
