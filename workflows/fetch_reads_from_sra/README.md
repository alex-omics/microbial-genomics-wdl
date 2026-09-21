# fetch_reads_from_sra

Fetch a single SRA/ENA/DDBJ run and get back gzip-compressed FASTQ, plus a record of what was
actually extracted, so a mixed batch of accessions (Illumina paired-end, ONT single-end, ...)
needs no sorting by hand afterwards.

## Terra usage

The workflow takes one accession per execution and relies on Terra's row-level scatter: point
`accession` at `this.SRA_ID` in a data table, select the rows to fetch, and Terra runs one
execution per row. Outputs write back to the same row (`this.read1`, `this.read2`,
`this.layout`, ...), so the table holds read paths and metadata next to your sample
identifiers.

To set it up: build a data table with one row per run accession and a column holding the
accession (`SRA_ID` by convention), run the workflow with `accession` wired to `this.SRA_ID`,
select every row, and launch. Terra will prompt you to map the output columns the first time.

## What it does

1. `prefetch`es the accession's `.sra` container (resumable, so a dropped connection costs
   the download rather than the unpack step).
2. `fasterq-dump --split-3`s it into FASTQ.
3. Detects what came out — paired (`_1`/`_2`), single (one file), or paired with a few
   leftover unpaired reads — instead of assuming paired-end. Single-end output is the norm
   for ONT/PacBio runs.
4. Looks up the run's platform and instrument from ENA's public metadata API (no credentials)
   and reports them as their own outputs, separate from `layout`, so a table can be filtered
   on either independently.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `accession` | `String` | **Required.** A single-run accession: `SRR…`, `ERR…`, or `DRR…` |
| `disk_gb` | `Int` | Default `50`, enough for MiSeq-scale Illumina. A full ONT/PacBio run can need much more |
| `cpu` | `Int` | Threads for `fasterq-dump` (default `4`) |
| `memory_gb` | `Int` | Default `8` |

### Only single-run accessions

`accession` must be a run accession (`SRR`/`ERR`/`DRR`), not an experiment (`SRX`/`ERX`),
sample (`SRS`), study (`SRP`) or BioProject (`PRJ*`) accession. Those can span multiple runs
and `sra-tools` does not expand them. The task checks the prefix up front and fails with a
clear message. To turn a study or experiment into run accessions, use ENA's `filereport` API
once per batch.

## Outputs

| Output | Type | Notes |
| ------ | ---- | ----- |
| `read1` | `File` | R1 for paired-end, or the only file for single-end |
| `read2` | `File?` | Unset for single-end runs |
| `orphan_reads` | `File?` | Reads `--split-3` could not pair, for an otherwise paired-end run. Usually absent |
| `layout` | `String` | `paired` or `single`, from the files actually produced rather than SRA/ENA metadata |
| `platform` | `String` | e.g. `ILLUMINA`, `OXFORD_NANOPORE`, `PACBIO_SMRT`, from ENA. `unknown` if the lookup fails or the run is not yet mirrored to ENA; this never blocks the FASTQ output |
| `instrument_model` | `String` | e.g. `Illumina MiSeq`, `MinION`. Same best-effort caveat as `platform` |
| `sra_tools_version` | `String` | Version of `sra-tools` used |

## Running

```bash
miniwdl run workflows/fetch_reads_from_sra/fetch_reads_from_sra.wdl -i tests/fetch_reads_from_sra_inputs.json
```

The same input file works locally and on Terra, since every input is a primitive.
