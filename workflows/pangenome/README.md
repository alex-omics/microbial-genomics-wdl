# pangenome

Bakta + Panaroo across a panel of assemblies, and stops there. **No alignment, no tree.**

**Status: 🚧 in development.** Validated against `miniwdl check` and shares its Bakta and
Panaroo tasks with `dna_methylation_calling`, which has been run through real Docker
execution — but this workflow itself has not yet run end to end on real data.

## Why it stops where it does

Pangenome construction and phylogenetic tree building are deliberately kept in separate
workflows. Coupling them means every pangenome run pays for a tree search — on a large,
mixed-quality panel that search can fail to converge and run up real cost with nothing to
show for it (this happened once, on a ~1,500-genome collection, for roughly $200).
Panaroo's own aligner (`--alignment core`/`pan`) is never invoked here by default —
`gene_presence_absence.csv`, which is all a downstream ortholog join needs, comes out
without it.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `assemblies` | `Array[File]` | Assembly FASTAs making up the panel. Ignored if `gff3s` is supplied. |
| `sample_names` | `Array[String]?` | Optional; positionally matched. Becomes the isolate column headers in `gene_presence_absence.csv`. |
| `gff3s`, `fnas` | `Array[File]?` | Pre-existing annotations, to skip Bakta. Positionally matched to each other. |
| `proteins` | `File?` | Trusted proteins for Bakta `--proteins`, e.g. the PAO1 proteome — keeps gene names consistent across independently-annotated isolates. |
| `bakta_db` | `File?` | Full Bakta database `.tar.gz`. Omitted, the light DB baked into the image is used. |
| `merge_paralogs` | `Boolean` | Off by default. Genuinely contested for multi-copy gene families — see the task's own `parameter_meta`. |
| `alignment` | `String?` | `"core"` or `"pan"` to also emit gene alignments. Leave unset unless something downstream actually consumes an alignment — see above. |
| `panaroo_mem_gb` | `Int` | Defaults to 64 GB, sized for a panel of a few dozen complete assemblies. Scale up hard for genus-scale panels (hundreds to low thousands of genomes) — an OOM mid-graph-build wastes the whole run plus the retry, which costs more than the extra RAM would have. |

## Outputs

| Output | Contents |
| ------ | -------- |
| `gene_presence_absence` | Panaroo's core deliverable — the ortholog-group × isolate table. This is the join key `dna_methylation_calling`'s ortholog stage consumes via its own `gene_presence_absence` input, to avoid recomputing the pangenome. |
| `gene_presence_absence_rtab` | Binary presence/absence matrix. |
| `summary_statistics`, `n_core_genes`, `n_total_genes` | Pangenome-wide summary counts. |
| `bakta_gff3`, `bakta_faa`, `bakta_tsv` | Per-isolate annotation, if Bakta ran here rather than being supplied via `gff3s`. |

## Running

```bash
miniwdl run workflows/pangenome/pangenome.wdl -i tests/pangenome_inputs.json
```

`tests/pangenome_inputs.json` is a placeholder for Dockstore's `testParameterFiles`, not
a ready-to-run local test — fill in real assembly paths, or configure inputs directly in
the Terra console.
