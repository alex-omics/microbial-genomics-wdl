# dna_methylation_calling

Per-isolate bacterial **DNA** methylation calling from ONT modified-basecalled BAMs, for
characterising the methylome landscape across a panel of isolates rather than confirming
any one gene's activity.

**Status: available.** Stages are tested against synthetic ground-truth fixtures, and the
full workflow has been run end to end on Terra.

## Design

Each isolate's reads are mapped to **its own assembly**, not a shared reference. No single
reference is adequate for an organism with a large accessory genome and frequent
rearrangement, and de novo motif discovery reads sequence context from whatever the reads
were aligned to, so a foreign reference gives the wrong context wherever the isolate differs.

Independent assemblies share no coordinates, so cross-isolate comparison uses **orthology**
instead, via a Panaroo pangenome built across the panel.

```
                    ┌─ modkit_find_motifs ─┐
modbam ─┬─ align ─┬─ modkit_pileup ────────┼─ build_motif_list ─ motif_landscape ─┐
        │         └─ (per isolate)         │      (tiers 1 + 2)                  │
assembly┘                                  │                                     │
                    bakta ──┬─ annotate ────┘                                     │
                             │                                                    │
                             └─ rebase_blastp → rebase_join_motifs ───────────────┘
                                                                                    │
   panaroo (pangenome) ── methylation_orthologs ── motif_landscape_summary ────────┘
                                                          (tier 3, panel-wide)
```

## Tools

| Stage | Tool | Image |
| ----- | ---- | ----- |
| Alignment | minimap2 2.31 + samtools 1.24 | `nanozoo/minimap2:2.31--c2b4c91` |
| Pileup + de novo motifs | modkit 0.6.4 | `quay.io/biocontainers/ont-modkit:0.6.4--h7f49ad2_0` |
| Annotation | Bakta 1.12.0 (light DB baked in) | `staphb/bakta:1.12.0-6.0-light` |
| Intersect / flank | bedtools 2.31.1 | `staphb/bedtools:2.31.1` |
| Pangenome | Panaroo 1.7.0 | `staphb/panaroo:1.7.0` |
| MTase homology | BLAST+ 2.16.0 | `staphb/blast:2.16.0` |
| Table joins, statistics | Python 3.11 | `python:3.11-slim` |

All images are pinned to a registry digest, not just a tag.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `modbams` | `Array[File]` | **Required.** Modified-basecalled BAMs (MM/ML tags), one per isolate. Unaligned is expected; an aligned input is reduced to primary records and realigned |
| `assemblies` | `Array[File]` | **Required.** Each isolate's own assembly, positionally matched to `modbams`. Not produced by this workflow |
| `sample_names` | `Array[String]?` | Positionally matched. Defaults to filenames with `.bam`/`.modbam` stripped |
| `proteins` | `File?` | Trusted proteins for Bakta `--proteins` (e.g. a reference proteome), to carry gene names across independently annotated isolates |
| `bakta_db` | `File?` | Full Bakta database `.tar.gz`. Omitted, the light database in the image is used |
| `run_bakta` / `run_find_motifs` / `run_pangenome` / `run_motif_landscape` | `Boolean` | All default `true`; each stage can be turned off |
| `rebase_goldset_fasta`, `rebase_motif_tsv` | `File?` | REBASE Gold Standard protein set and its motif table, for MTase homology identification. Both are required to run it; omit either to skip it. Not included in this repo (REBASE carries its own licensing terms) |
| `gene_presence_absence` | `File?` | An existing Panaroo `gene_presence_absence.csv`, to skip building the pangenome |
| `panaroo_clean_mode`, `core_threshold`, `merge_paralogs` | | Pangenome settings; see the task's `parameter_meta` for the `merge_paralogs` tradeoff |
| `min_coverage`, `min_percent`, `min_mod_reads` | | The three floors a site must clear to count as methylated. `min_percent` defaults to 50; tune it to your data. Re-thresholding reruns only annotation and the ortholog join |
| `min_mapped_percent` | `Float` | Self-mapping should exceed ~95%. The default 85% floor catches a `modbams`/`assemblies` pair passed in mismatched order |
| `flank_upstream`, `trim_to_intergenic`, `feature_type` | | Upstream/promoter window definition for annotation |
| `genus`, `species` | `String` | Passed to Bakta (default *Pseudomonas aeruginosa*) |
| `heterogeneous_low_cutoff`, `highlight_genes` | | Tier 2/3 landscape settings, below |

## The three-tier landscape

`run_motif_landscape` tests every candidate motif (de novo from modkit and/or predicted from
REBASE homology) against each isolate's own data, then summarises across the panel. All three
tiers are descriptive and need no phenotype data.

1. **Enrichment** — whether methylation at a motif's occurrences is elevated over the
   genome-wide background for that base (log2 enrichment + chi-square).
2. **Within-genome heterogeneity** — whether a motif's copies in one genome are uniformly
   methylated or split into high and low populations, which suggests phase variation or a
   regulator competing for specific copies.
3. **Cross-isolate variability** — by motif, and by gene via the pangenome. Uniform
   methylation across the panel is likely restriction-modification housekeeping; sharp
   variation is where isolate-specific biology shows up. `highlight_genes` (e.g.
   `["mex","opr","amp","nal"]`) flags matching genes in the gene-level table without
   filtering anything out.

## Outputs

| Output | Contents |
| ------ | -------- |
| `summary_tsv` | One row per isolate: mapping rate, depth, site counts by modification type, motif/gene counts |
| `methylation_long` | Every isolate's annotated methylation sites, concatenated, at nucleotide resolution |
| `bedmethyl_gz`, `bed_6ma`, `bed_4mc`, `bed_5mc` | Unfiltered per-isolate pileups split by modification type, for re-thresholding |
| `annotated_tables` | Per isolate: one row per methylated site, with `locus_tag`/`gene`/`product`/`region` (genic or upstream) |
| `annotated_tables_with_motifs` | `annotated_tables` plus a `motif` column: which of the isolate's candidate motifs (de novo and REBASE), if any, the site falls inside |
| `ortholog_matrix`, `ortholog_long` | Gene x isolate methylation. A gene an isolate lacks is `NA`, distinct from `0` (present, unmethylated) |
| `rebase_mtases` | Per-isolate candidate MTase inventory: homologous gene, predicted motif, modification type |
| `motif_landscapes` | Per-isolate tier 1 + 2 statistics, one row per motif tested |
| `motif_summary`, `gene_landscape_summary` | Tier 3 cross-isolate variability rankings, by motif and by gene |
| `bakta_gff3`, `bakta_faa`, `bakta_tsv` | Per-isolate annotation |

## Running

```bash
miniwdl run workflows/dna_methylation_calling/dna_methylation_calling.wdl -i inputs.json
```

where `inputs.json` sets at least `modbams` and `assemblies`, or set the inputs in the Terra
console.

`tests/run_unit_tests.sh` holds the regression suite: task-level fixture tests against
planted ground truth, run through the real WDL tasks via `miniwdl run`.
