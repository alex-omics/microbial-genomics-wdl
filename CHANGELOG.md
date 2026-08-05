# Changelog

## [Unreleased]
### Added
- `workflows/dna_methylation_calling/dna_methylation_calling.wdl` - Per-isolate bacterial methylation calling from ONT modified-basecalled BAMs. Each isolate's reads are mapped to its OWN assembly rather than a shared reference — no single reference captures the accessory genome of an organism this prone to rearrangement, and de novo motif discovery specifically requires it, since motif context is read from whatever the reads were aligned to
- `tasks/align_modbam.wdl` - `samtools fastq -T MM,ML,MN | minimap2 -y`, with MM-tag presence checked before and after alignment; losing the tags produces an empty pileup with no error otherwise
- `tasks/modkit.wdl` - `modkit_pileup` (bedMethyl split by 6mA/4mC/5mC) and `modkit_find_motifs` (de novo motif discovery)
- `tasks/bakta.wdl` - Per-isolate annotation with `--keep-contig-headers` hardcoded (not exposed as an input — a renamed contig silently empties every downstream intersect) and `--proteins` support for carrying reference gene names onto independently-annotated isolates
- `tasks/annotate_methylation.wdl` - Joins bedMethyl to the isolate's own Bakta GFF3 (genic/upstream, three floors: coverage, percent-modified, and modified-read count) with `region_length` carried through so site counts can be normalised into a density rather than a proxy for gene length
- `tasks/panaroo.wdl`, `tasks/methylation_orthologs.wdl` - Pangenome-based ortholog join across independently-assembled isolates. Absence is `NA`, not `0`, in the resulting matrix — collapsing the two would manufacture methylation differences that are really gene presence/absence
- `workflows/pangenome/pangenome.wdl` - Bakta + Panaroo, stopping there. No alignment, no tree — coupling pangenome construction to phylogenetics on a large or mixed-quality panel risks paying for a tree search that never converges
- `tasks/rebase_mtase_search.wdl` - MTase homology identification against the REBASE Gold Standard set (BLASTP + motif join, split into two tasks since the BLAST image carries no Python)
- `tasks/motif_landscape.wdl`, `tasks/motif_landscape_summary.wdl` - Three-tier descriptive methylome characterisation: per-isolate motif enrichment vs. genome background, within-genome heterogeneity (uniform vs. bimodal — the phase-variation/regulator-competition signature), and cross-isolate variability by motif and by gene. Deliberately descriptive, not a phenotype association test — the goal is landscape characterisation across a diverse panel, not confirming any one enzyme's activity
- `tasks/utils.wdl` - Shared `name_summary`, `concat_tables`, and `validate_panel` (fails fast on mismatched input arrays or sample names that would break locus tags or pangenome column headers, before the expensive stages run)
- `tests/run_unit_tests.sh` - Fixture-based regression suite for the methylation tasks, run through the real WDL tasks via `miniwdl run` rather than a copy of their logic. Statistics were validated against a synthetic genome with known planted methylation before being fixed as assertions
- `.dockstore.yml` - Registered `dna_methylation_calling` and `pangenome`
- `multimodal_gwas` - Terra WDL translation of microGWAS (Burgaya et al. 2025), targeting faithful reproduction of the published *E. coli* virulence GWAS as Milestone 1
- `workflows/assembly_qc/assembly_qc.wdl` - Scatters QUAST, BUSCO, and CheckM2 over a set of assemblies and emits one merged summary TSV
- `tasks/checkm2.wdl` - CheckM2 completeness and contamination, with the DIAMOND database supplied as a `File` input (Terra-native) or downloaded at runtime
- `busco-prokaryota:5.8.0` - Docker image with the complete OrthoDB v10 prokaryote set (99 datasets) plus placement files for offline auto-lineage
- `tasks/busco.wdl` - Reusable BUSCO task with auto-lineage and explicit lineage support, plus an optional lineage tarball override

### Fixed
- `tasks/annotate_methylation.wdl` - Was filtering on coverage only, so every evaluated base reached the output table; the ortholog matrix built on top of it was effectively reporting gene length rather than methylation. Now filters on coverage, percent-modified, and modified-read count together
- `tasks/annotate_methylation.wdl` - `bedtools intersect -wa -wb` column offset was wrong by one field, writing the feature's end coordinate into the `feature_strand` column and the strand into `product`; caught by hand-verifying the real column layout against bedtools 2.31.1, not by inspection
- `tasks/motif_landscape_summary.wdl` - Sorted motif rows by a formatted string (`"%.3f" % cv`) instead of the underlying float; crashed on first real invocation, caught only by actually running the task
- `tasks/rebase_mtase_search.wdl` - `evalue` was a `Float` default of `1e-25`; WDL's string interpolation renders a float that small as the fixed-decimal text `"0.000000"`, which BLASTP then rejects as a non-positive cutoff. Changed to a `String` — this would have failed identically on Terra, not just in local testing
- `tasks/busco.wdl` - Pointed `--download_path` at the consolidated `/busco_downloads`; the previous path held no manifest, so auto-lineage could not run and only spirochaete datasets resolved
- `tasks/busco.wdl` - `grep` for the lineage-used line had the filename inside the pattern with no file argument, so it read stdin and returned nothing
- `tasks/busco.wdl` - Summary was copied from `<sample>_busco/` while BUSCO wrote to `<sample>/`; auto-lineage also emits several `short_summary*` files, which broke the glob
- `tasks/busco.wdl` - Runtime used `disk:` / `local_disk`; Terra requires `disks:` / `local-disk`, so the attribute was silently ignored
- `tasks/quast.wdl` - Invoked `quast`, which does not exist in `staphb/quast:5.2.0` (the entrypoint is `quast.py`)
- `tasks/quast.wdl` - Quoted domain flag expanded to an empty positional argument for prokaryotes, which QUAST read as an empty input filename

## [0.1.0] — 17 June 2026
- Initial repo structure
