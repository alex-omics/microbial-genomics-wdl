# assembly_qc

Assess the quality and completeness of one or more assembled microbial genomes and
collapse everything into a **single summary table**.

The problem this solves: assembly statistics and the assembly FASTA usually end up
decoupled across Terra tables, and re-sequenced isolates multiply the number of places
you have to look. Point this workflow at a pile of assemblies gathered from anywhere and
get one TSV back.

## Tools

| Tool | What it reports | Image |
| ---- | --------------- | ----- |
| QUAST 5.2.0 | Contiguity: contig count, total length, N50/L50, GC% | `staphb/quast:5.2.0` |
| BUSCO 5.7.1 | Completeness against single-copy orthologs | `aarvani1/busco-prokaryota:5.8.0` |
| CheckM2 1.1.0 | Completeness and **contamination** via ML models | `staphb/checkm2:1.1.0` |

BUSCO and CheckM2 both report completeness but disagree usefully: BUSCO counts marker
genes for a specific clade, CheckM2 predicts from a model trained across genomes and is
the one that will tell you an assembly is *contaminated*.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `assemblies` | `Array[File]` | **Required.** Assembly FASTAs (`.fasta`/`.fa`/`.fna`, optionally gzipped) |
| `sample_names` | `Array[String]?` | Optional; positionally matched. Defaults to filenames with the extension stripped |
| `busco_lineage` | `String` | Dataset name, or `auto` (default) to place each assembly independently |
| `busco_lineage_tarball` | `File?` | A lineage `.tar.gz` not baked into the image |
| `checkm2_db` | `File?` | **Recommended.** See below |
| `reference_genome` | `File?` | Enables QUAST comparative mode |
| `run_quast` / `run_busco` / `run_checkm2` | `Boolean` | All default `true` |

### BUSCO lineages

The `busco-prokaryota:5.8.0` image carries the **complete OrthoDB v10 prokaryote set —
99 datasets** — plus the placement files needed for `--auto-lineage-prok`. Everything
runs fully offline; the task passes `--offline` unconditionally.

`busco_lineage = "auto"` lets BUSCO place each assembly in the prokaryote tree and pick
its own dataset, which is what you want for a mixed bag of isolates. The dataset it
actually chose comes back in the `busco_lineage` column, so the number is never
uninterpretable. Set an explicit name (e.g. `spirochaetales_odb10`) when you know the
clade and want it uniform across samples.

> The image pins odb10 deliberately. Upstream now also publishes ~740 OrthoDB v12
> datasets, which BUSCO 5.7.x cannot consume — `busco --download prokaryota` against the
> live manifest pulls odb12 and fails. See `docker/busco-prokaryota/fetch_busco_data.sh`.

### The CheckM2 database

`staphb/checkm2` ships **without** the ~3 GB DIAMOND database. Stage it in GCS once:

```bash
docker run --rm -v "$PWD:/db" staphb/checkm2:1.1.0 checkm2 database --download --path /db
```

Then upload `CheckM2_database/uniref100.KO.1.dmnd` to your bucket and pass that path as
`checkm2_db`. If you omit it, **every scattered task downloads its own copy**, which is
slow and depends on an external host staying up.

## Outputs

`summary_tsv` is the deliverable — one row per assembly:

```
sample  quast_contigs  quast_total_length  quast_largest_contig  quast_n50  quast_l50
quast_gc_percent  busco_complete_pct  busco_single_copy_pct  busco_duplicated_pct
busco_fragmented_pct  busco_missing_pct  busco_n_markers  busco_lineage
checkm2_completeness  checkm2_contamination  checkm2_model  checkm2_coding_density
```

Disabled tools leave `NA` rather than dropping columns, so the table stays rectangular.
Per-sample reports (`quast_report_html`, `busco_summary_txt`, `busco_full_table`,
`checkm2_report_tsv`) come back as arrays for when a row looks wrong and you need detail.

## Running

```bash
miniwdl run workflows/assembly_qc/assembly_qc.wdl -i tests/assembly_qc_inputs.json
```
