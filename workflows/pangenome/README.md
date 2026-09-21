# pangenome

Bakta + Panaroo across a panel of assemblies. **No alignment, no tree.**

**Status: in development.** The Bakta and Panaroo tasks have been run on real data as part of
other workflows; this workflow has not yet been run as its own Terra entry point.

## Scope

Pangenome construction and phylogenetic tree building are kept in separate workflows, because
a tree search on a large or mixed-quality panel can run for a long time without converging.
Panaroo's aligner (`--alignment core`/`pan`) is not invoked by default: `gene_presence_absence.csv`,
which is what an ortholog join needs, does not require it.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `assemblies` | `Array[File]` | Assembly FASTAs making up the panel. Ignored if `gff3s` is supplied |
| `sample_names` | `Array[String]?` | Positionally matched. Becomes the isolate column headers in `gene_presence_absence.csv` |
| `gff3s`, `fnas` | `Array[File]?` | Existing annotations, to skip Bakta. Positionally matched to each other |
| `proteins` | `File?` | Trusted proteins for Bakta `--proteins` (e.g. a reference proteome), to keep gene names consistent across isolates |
| `bakta_db` | `File?` | Full Bakta database `.tar.gz`. Omitted, the light database in the image is used |
| `merge_paralogs` | `Boolean` | Off by default. See the task's `parameter_meta` for the tradeoff |
| `alignment` | `String?` | `"core"` or `"pan"` to also emit gene alignments. Leave unset unless something needs one |
| `panaroo_mem_gb` | `Int` | Default 64 GB, sized for a few dozen complete assemblies. Scale up for panels of hundreds or more; running out of memory mid-build wastes the whole run |

## Outputs

| Output | Contents |
| ------ | -------- |
| `gene_presence_absence` | Panaroo's ortholog group x isolate table (`.csv`, with each isolate's locus tags) |
| `gene_presence_absence_rtab` | Binary presence/absence matrix |
| `summary_statistics`, `n_core_genes`, `n_total_genes` | Pangenome-wide summary counts |
| `bakta_gff3`, `bakta_faa`, `bakta_tsv` | Per-isolate annotation, if Bakta ran here rather than being supplied via `gff3s` |

## Running

```bash
miniwdl run workflows/pangenome/pangenome.wdl -i inputs.json
```

where `inputs.json` sets `pangenome.assemblies` (and optionally `pangenome.sample_names`),
or set the inputs in the Terra console.
