# Changelog

## [Unreleased]
### Added
- `workflows/rna_seq_counts` - Paired-end RNA-seq alignment (BWA-MEM) and featureCounts, with each sample aligned to its own isolate's assembly or to one shared reference. Annotates with Bakta and builds an ortholog table with Panaroo, each skipped if supplied, and emits per-isolate matrices, a long-format table and one ortholog-level matrix
- `tasks/sample_matching.wdl` - Matches samples to assemblies by name, allowing replicate suffixes (`a/b/c`, `1/2/3`, `rep2`, ...), and fails before alignment if a sample cannot be matched. Also validates a supplied ortholog table
- `tasks/bwa.wdl`, `tasks/picard.wdl`, `tasks/featurecounts.wdl` - BWA-MEM alignment, duplicate marking, per-isolate counting, and a merge into a long table and an ortholog-level matrix. Absent genes are `NA`, not `0`; counted features missing from the ortholog table are kept as `unclustered` rows
- `workflows/dna_methylation_calling` - Per-isolate bacterial methylation calling from ONT modified-basecalled BAMs, with each isolate mapped to its own assembly, plus a three-tier descriptive motif landscape across the panel
- `tasks/align_modbam.wdl`, `tasks/modkit.wdl`, `tasks/annotate_methylation.wdl`, `tasks/annotate_motifs.wdl` - Modified-base alignment and pileup, de novo motif discovery, annotation of methylated sites with gene and motif context
- `tasks/bakta.wdl` - Per-isolate Bakta annotation, with `--keep-contig-headers` fixed on and optional `--proteins`
- `tasks/panaroo.wdl`, `tasks/methylation_orthologs.wdl` - Pangenome-based ortholog join across independently assembled isolates
- `tasks/rebase_mtase_search.wdl`, `tasks/motif_landscape.wdl`, `tasks/motif_landscape_summary.wdl` - MTase homology identification against REBASE and the motif landscape statistics
- `tasks/utils.wdl` - Shared `name_summary`, `concat_tables` and `validate_panel`
- `workflows/pangenome` - Bakta + Panaroo, with no alignment or tree
- `workflows/fetch_reads_from_sra` - Fetches one SRA/ENA/DDBJ run as gzip-compressed FASTQ, detecting the layout that was actually produced and reporting platform and instrument
- `workflows/pyseer_gwas` - pyseer LMM association over a presence/absence Rtab, with optional lineage effects, a covariate combination scan and a gene/module annotation join
- `workflows/assembly_qc` - QUAST, BUSCO and CheckM2 over a set of assemblies, merged into one summary table
- `tasks/busco.wdl`, `tasks/checkm2.wdl`, `tasks/quast.wdl` - Assembly QC tasks
- `docker/busco-prokaryota`, `docker/pyseer` - Docker images for BUSCO (OrthoDB v10 prokaryote set, offline) and pyseer
- `tests/run_unit_tests.sh` - Fixture-based regression suite that runs the real WDL tasks through `miniwdl`
- `.dockstore.yml` - Dockstore registration

### Changed
- `tasks/panaroo.wdl` - New optional `sample_names` input, so isolate column names in `gene_presence_absence.csv` do not depend on GFF filenames
- `workflows/pyseer_gwas` - Removed the separate annotation of `gene_significant`, which is a subset of `gene_results`. The covariate scan's long-format table is now annotated too
- `tasks/modkit.wdl` - `modkit_find_motifs` is no longer preemptible, and defaults to 32 CPUs
- `tasks/motif_landscape_summary.wdl` - Gene-level ranking breaks CV ties on mean density

### Fixed
- `rna_seq_counts` - `-p` without `--countReadPairs` counted each mate separately, doubling every fragment. The pinned images for alignment, duplicate marking and counting could not run (`staphb/bwa` has no samtools; `staphb/picard` and `biocontainers/subread` do not exist on Docker Hub). `ignore_duplicates` now defaults to `false`
- `tasks/annotate_methylation.wdl` - Sites were filtered on coverage only; now on coverage, percent-modified and modified-read count. Also fixed a column offset that put the feature end coordinate in `feature_strand`
- `tasks/motif_landscape_summary.wdl` - Motif rows were sorted by a formatted string rather than the numeric value
- `tasks/motif_landscape.wdl` - Replaced a linear scan per site with a binary search; the scan could take up to two hours on some isolates
- `tasks/rebase_mtase_search.wdl` - `evalue` is now a `String`, because WDL renders a very small `Float` as `0.000000`. REBASE entries with several comma-joined motifs are split per motif
- `tasks/modkit.wdl` - The `motifs` summary column read the modification code instead of the motif sequence
- `tasks/pyseer.wdl` - The covariate scan failed on files with Windows line endings. Lineage effects now run as a separate task on a small slice of the Rtab, since combining them with the full `--lmm` call ran out of memory, and they get the distance matrix pyseer requires
- `tasks/busco.wdl` - Corrected the download path, lineage-used lookup, summary file paths, and the runtime disk attribute names Terra requires
- `tasks/quast.wdl` - Corrected the executable name (`quast.py`) and an empty argument passed for prokaryotes

## [0.1.0] — 17 June 2026
- Initial repo structure
