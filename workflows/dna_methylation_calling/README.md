# dna_methylation_calling

Per-isolate bacterial **DNA** methylation calling from ONT modified-basecalled BAMs,
built for characterizing the methylome landscape across a diverse panel of clinical
isolates — not for confirming any one gene's activity. (Named explicitly `dna_` because
"methylation" alone is ambiguous — this is not histone or RNA methylation.)

**Status: 🚧 in development.** Every stage has been validated against synthetic
ground-truth fixtures and real Docker execution, but the pipeline has not yet run on a
real ONT BAM or a real assembly end to end.

## Design

Each isolate's reads are mapped to **its own assembly**, not a shared reference. No
single reference is adequate for an organism with substantial accessory genome and
frequent rearrangement, and de novo motif discovery specifically requires it: motif
context is read from whatever the reads were aligned to, so aligning to a foreign genome
reads the wrong sequence at every position where the isolate differs from it.

Because assemblies are independent, cross-isolate comparison cannot use coordinates —
it uses **orthology** instead, via a Panaroo pangenome built across the panel.

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
| `modbams` | `Array[File]` | **Required.** Modified-basecalled BAMs (MM/ML tags), one per isolate. Unaligned is expected; an aligned input is reduced to primary records and realigned. |
| `assemblies` | `Array[File]` | **Required.** Each isolate's own assembly, positionally matched to `modbams`. Not produced here — bring your own from TheiaProk ONT, Autocycler, etc. |
| `sample_names` | `Array[String]?` | Optional; positionally matched. Defaults to filenames with `.bam`/`.modbam` stripped. |
| `proteins` | `File?` | Trusted proteins for Bakta `--proteins`, e.g. the PAO1 proteome — transfers gene names across independently-annotated isolates. |
| `bakta_db` | `File?` | Full Bakta database `.tar.gz`. Omitted, the light DB baked into the image is used (no staging required). |
| `run_bakta` / `run_find_motifs` / `run_pangenome` / `run_motif_landscape` | `Boolean` | All default `true`. Every stage is independently toggleable. |
| `rebase_goldset_fasta`, `rebase_motif_tsv` | `File?` | REBASE Gold Standard protein set + its motif table, for MTase homology ID. Both required to run it; omitting either skips it cleanly. Stage from MPore's repo (`DiltheyLab/MPore`) or NEB directly — deliberately not vendored in this repo (REBASE and MPore both carry their own licensing terms). |
| `gene_presence_absence` | `File?` | Reuse a pangenome already built by the standalone `pangenome` workflow instead of recomputing it here. |
| `panaroo_clean_mode`, `core_threshold`, `merge_paralogs` | — | Pangenome construction knobs. `merge_paralogs` is genuinely contested — see the task's own `parameter_meta` for the tradeoff. |
| `min_coverage`, `min_percent`, `min_mod_reads` | — | Three floors a site must clear to count as methylated. `min_percent` (default 50) is a placeholder pending literature review, not a settled threshold — re-thresholding only reruns annotation and the ortholog join, not alignment or pileup. |
| `min_mapped_percent` | `Float` | Self-mapping should clear ~95%+; the default 85% floor exists to catch a `modbams`/`assemblies` pair passed in mismatched order. |
| `flank_upstream`, `trim_to_intergenic`, `feature_type` | — | Upstream/promoter window definition for annotation. |
| `genus`, `species` | `String` | Passed to Bakta. Default `Pseudomonas aeruginosa`. |
| `heterogeneous_low_cutoff`, `highlight_genes` | — | Tier 2/3 landscape knobs — see below. |

## The three-tier landscape

`run_motif_landscape` tests every candidate motif (de novo from modkit, and/or
homology-predicted from REBASE) against each isolate's own data, then summarizes across
the panel. All three tiers are descriptive — none require MICs or phenotype groups.

1. **Enrichment** — is methylation at a motif's genomic occurrences elevated over the
   genome-wide background for that base (log2 enrichment + chi-square).
2. **Within-genome heterogeneity** — across one motif's many copies in a single genome,
   uniform (the RM-housekeeping null) or split into high/low populations, the signature
   of phase variation or a regulator competing for specific copies.
3. **Cross-isolate variability** — by motif, and by gene (via the pangenome). A motif or
   gene that is uniform across the whole panel is very likely RM housekeeping; one that
   varies sharply is where isolate-specific biology would show up. `highlight_genes`
   (e.g. `["mex","opr","amp","nal"]`) flags matches in the gene-level table without
   filtering anything out.

## Outputs

| Output | Contents |
| ------ | -------- |
| `summary_tsv` | One row per isolate: mapping rate, depth, site counts by modification type, motif/gene counts. |
| `methylation_long` | Every isolate's annotated methylation sites concatenated, nucleotide-resolution. |
| `bedmethyl_gz`, `bed_6ma`, `bed_4mc`, `bed_5mc` | Unfiltered per-isolate pileups, split by modification type — the re-thresholding fallback if the built-in filters are too strict or loose. |
| `annotated_tables` | Per-isolate: one row per methylated site, with `locus_tag`/`gene`/`product`/`region` (genic vs. upstream). |
| `ortholog_matrix`, `ortholog_long` | Gene × isolate methylation, gene absence held as `NA`, distinct from `0` (present, unmethylated). |
| `rebase_mtases` | Per-isolate candidate MTase inventory: homologous gene, predicted motif, modification type. |
| `motif_landscapes` | Per-isolate tier 1 + 2 statistics, one row per motif tested. |
| `motif_summary`, `gene_landscape_summary` | Tier 3: cross-isolate variability rankings, by motif and by gene. |
| `bakta_gff3`, `bakta_faa`, `bakta_tsv` | Per-isolate annotation, for anything downstream that needs it directly. |

## Running

```bash
miniwdl run workflows/dna_methylation_calling/dna_methylation_calling.wdl -i tests/dna_methylation_calling_inputs.json
```

`tests/dna_methylation_calling_inputs.json` is a placeholder — it exists to satisfy
Dockstore's `testParameterFiles`, not as a ready-to-run local test (there is no small
public ONT modBAM + assembly fixture checked in). Fill in real paths, or configure inputs
directly in the Terra console.

See `tests/run_unit_tests.sh` for the actual regression suite: task-level fixture tests
against planted ground truth, run through the real WDL tasks via `miniwdl run`.
