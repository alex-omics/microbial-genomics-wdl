# fetch_reads_from_sra

Fetch a single SRA/ENA/DDBJ run and get back ready-to-use, gzip-compressed FASTQ — plus a
record of what was actually extracted, so a mixed batch of accessions (some Illumina PE,
some ONT SE) doesn't need to be sorted out by hand afterward.

## Terra usage

This workflow is built around Terra's row-level scatter, not an `Array` input, because that
is the pattern it's actually run with: point `accession` at `this.SRA_ID` in a data table,
select every row you want fetched, and Terra scatters one workflow execution per row on its
own. Outputs write straight back to the same row — `this.read1`, `this.read2`,
`this.layout`, etc. — so the table ends up with paths and metadata sitting right next to
whatever sample identifiers you already had. This is validated against real use: a table of
300+ accessions run this way, unmodified, with a 100% success rate.

If you're setting this up fresh: build a data table with one row per SRA run accession
(entity name can be your own sample ID), a column holding the accession
(`SRA_ID` by convention), run the workflow with `accession` wired to `this.SRA_ID`, select
every row, and launch. New output columns will prompt you to map them (or auto-create them)
the first time you run a version of this workflow that adds one.

## What it does

1. `prefetch`s the accession's `.sra` container (resumable — a dropped connection costs the
   download, not a redo of the unpack step).
2. `fasterq-dump --split-3`s it into FASTQ.
3. Detects what actually came out — paired (`_1`/`_2`), single (one file), or paired with a
   handful of leftover unpaired reads — rather than assuming paired-end and failing on
   anything else. This matters in practice: single-end output is the norm for ONT/PacBio
   long-read runs, not an edge case.
4. Looks up the run's platform and instrument from ENA's public metadata API (no
   credentials needed) and reports it alongside the layout it detected from the actual
   files, so you can tell at a glance whether a given row is `ILLUMINA_paired`,
   `OXFORD_NANOPORE_single`, etc.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `accession` | `String` | **Required.** A single-run accession: `SRR…`, `ERR…`, or `DRR…` |
| `disk_gb` | `Int` | Default `50` — plenty for MiSeq-scale Illumina. A full ONT/PacBio flowcell run can be much larger; size up accordingly |
| `cpu` | `Int` | Threads for `fasterq-dump` (default `4`) |
| `memory_gb` | `Int` | Default `8` |

### Only single-run accessions

`accession` must be a run accession (`SRR`/`ERR`/`DRR`) — not an experiment (`SRX`/`ERX`),
sample (`SRS`), study (`SRP`), or BioProject (`PRJ*`) accession. Those can each span
multiple runs, and `sra-tools` does not auto-expand them; feeding one in either fails
outright or silently grabs the wrong thing. The task checks the accession's prefix up front
and fails immediately with a clear message rather than partway through a `prefetch` call.

This is deliberate, not a missing feature: resolving a study/experiment down to its
constituent run accessions is a one-time lookup against ENA's `filereport` API, done once
per batch rather than per row, and doesn't need to live inside this workflow. If that
becomes a recurring need, it's a small addition — a separate lookup step feeding a list of
run accessions into this workflow — rather than a reason to change this workflow's shape.

## Outputs

| Output | Type | Notes |
| ------ | ---- | ----- |
| `read1` | `File` | R1 for paired-end, or the only file for single-end |
| `read2` | `File?` | Unset for single-end runs |
| `orphan_reads` | `File?` | Reads `--split-3` couldn't pair, for an otherwise paired-end run. Usually absent |
| `layout` | `String` | `paired` or `single`, determined from which files were actually produced — not from (sometimes stale) SRA/ENA metadata |
| `platform` | `String` | e.g. `ILLUMINA`, `OXFORD_NANOPORE`, `PACBIO_SMRT`, from ENA's metadata. `unknown` if the lookup fails or the run isn't mirrored to ENA yet — this never blocks the fastq output, which is the point of the workflow |
| `instrument_model` | `String` | e.g. `Illumina MiSeq`, `MinION`. Same best-effort caveat as `platform` |
| `read_format` | `String` | Convenience combination of the two above, e.g. `ILLUMINA_paired`, `OXFORD_NANOPORE_single` |
| `sra_tools_version` | `String` | Version of `sra-tools` used, from the pinned image |

## Running

```bash
miniwdl run workflows/fetch_reads_from_sra/fetch_reads_from_sra.wdl -i tests/fetch_reads_from_sra_inputs.json
```

This one input file works for both local (`miniwdl`/Cromwell) and Terra runs — every input
is a primitive (`String`/`Int`), so there's no `gs://` vs. local-path split to make.
