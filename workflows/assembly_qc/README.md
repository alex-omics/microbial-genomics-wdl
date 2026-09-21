# assembly_qc

Assess the quality and completeness of assembled microbial genomes and collapse the results
into a **single summary table**, one row per assembly.

## Tools

| Tool | What it reports | Image |
| ---- | --------------- | ----- |
| QUAST 5.2.0 | Contiguity: contig count, total length, N50/L50, GC% | `staphb/quast:5.2.0` |
| BUSCO 5.7.1 | Completeness against single-copy orthologs | `aarvani1/busco-prokaryota:5.8.0` |
| CheckM2 1.1.0 | Completeness and **contamination** via ML models | `staphb/checkm2:1.1.0` |

BUSCO and CheckM2 both report completeness but answer different questions: BUSCO counts marker
genes for a specific clade, while CheckM2 predicts from a model trained across genomes and can
flag a *contaminated* assembly.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `assemblies` | `Array[File]` | **Required.** Assembly FASTAs (`.fasta`/`.fa`/`.fna`, optionally gzipped) |
| `sample_names` | `Array[String]?` | Positionally matched. Defaults to filenames with the extension stripped |
| `busco_lineage` | `String` | Dataset name, or `auto` (default) to place each assembly independently |
| `busco_lineage_tarball` | `File?` | A lineage `.tar.gz` not included in the image |
| `checkm2_db` | `File?` | **Recommended.** See below |
| `reference_genome` | `File?` | Enables QUAST comparative mode |
| `run_quast` / `run_busco` / `run_checkm2` | `Boolean` | All default `true` |

### BUSCO lineages

The `busco-prokaryota:5.8.0` image includes the complete OrthoDB v10 prokaryote set (99
datasets) and the placement files for `--auto-lineage-prok`. It runs fully offline.

With `busco_lineage = "auto"`, BUSCO places each assembly in the prokaryote tree and chooses
its own dataset, which suits a mixed set of isolates. The dataset chosen is reported in the
`busco_lineage` column. Set an explicit name (e.g. `spirochaetales_odb10`) to use one dataset
for every sample.

The image uses odb10 because BUSCO 5.7.x cannot use OrthoDB v12 datasets. See
`docker/busco-prokaryota/fetch_busco_data.sh`.

### The CheckM2 database

`staphb/checkm2` does not include the ~3 GB DIAMOND database. Download it once:

```bash
docker run --rm -v "$PWD:/db" staphb/checkm2:1.1.0 checkm2 database --download --path /db
```

Upload `CheckM2_database/uniref100.KO.1.dmnd` to your bucket and pass its path as `checkm2_db`.
Without it, every scattered task downloads its own copy.

## Outputs

`summary_tsv` has one row per assembly:

```
sample  quast_contigs  quast_total_length  quast_largest_contig  quast_n50  quast_l50
quast_gc_percent  busco_complete_pct  busco_single_copy_pct  busco_duplicated_pct
busco_fragmented_pct  busco_missing_pct  busco_n_markers  busco_lineage
checkm2_completeness  checkm2_contamination  checkm2_model  checkm2_coding_density
```

Columns for disabled tools are `NA`, so the table stays rectangular. Per-sample reports
(`quast_report_html`, `busco_summary_txt`, `busco_full_table`, `checkm2_report_tsv`) are
returned as arrays.

## Running

```bash
miniwdl run workflows/assembly_qc/assembly_qc.wdl -i tests/assembly_qc_inputs.json
```
