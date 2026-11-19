# Changelog

## [Unreleased]
### Added
- `tasks/validate_annotation.wdl` - `validate_annotation` stops an isolate before alignment if its assembly and annotation do not belong together: the annotator's FASTA must match the assembly's contig names and lengths, the GFF3 may name only those contigs and must agree with its own `##sequence-region` lengths, and the proteins must be features of the GFF3. Catches a GFF from another isolate even where the annotator numbered every isolate's contigs and genes alike, which otherwise labels every methylated site with the wrong isolate's genes. Also reports contigs, genome size, GC% and CDS count
- `tasks/modkit.wdl` - `modkit_pileup` reports every modification code in the pileup (`mod_codes_seen`) and warns about any other than `a`, `21839` and `m`, which were previously dropped from every count without a trace
- `tests/run_unit_tests.sh` - `validate_annotation` section (clean set and each mismatch, including identical contig names and gene numbering across isolates), composition-independence of motif scoring on a planted bipartite motif in 38% and 66% GC genomes, offset scoring leaving GATC unchanged, `scoring` in the tier 3 table, and the mod-code report
- `workflows/filter_reads` - Length and mean-quality filtering of one sample's ONT reads, from an unaligned BAM or FASTQ. A BAM is filtered record for record, so a modBAM keeps its MM/ML tags; FASTQ can also be cropped at either end. Emits the filtered reads in the input's format, a before/after table (reads, bases, N50, quality) and retention figures, and fails if nothing, or less than `min_reads_out` / `min_bases_out`, survives. Does no taxonomic or reference-based filtering
- `tasks/filter_reads.wdl` - `filter_reads_bam` (samtools; mean quality from dorado's `qs` tag, or computed with chopper when records lack it; refuses an aligned BAM), `filter_reads_fastq` (chopper), `read_stats` and `summarize_filtering` (seqkit), and `validate_filter_inputs`, which rejects a call that would filter nothing, crop a BAM or name two read inputs
- `docker/read-filter` - samtools 1.24, chopper 0.14.1 and seqkit 2.14.0 in one image
- `tests/run_unit_tests.sh` - `filter_reads` section: eight planted reads as a modBAM (with and without `qs`) and as FASTQ, covering both quality routes, inclusive length bounds, cropping and every refusal
- `workflows/call_ospA_type` - OspA in silico type (IST) of one *Borrelia burgdorferi* s.l. isolate with Pfizer's LISTT, from its assembly and, optionally, its Illumina paired-end reads. One sample per execution with no `Array` inputs. Emits one call (`IST<n>`, `Unknown` or `Ambig`) and a flag (`OK`, `SM`, `UK`, `DS`, `QF`), plus species, identity, details, a results table, the consensus proteins and LISTT's raw output
- `tasks/listt.wdl` - `bbsl_ospA_listt` runs one LISTT mode and classifies its consensus; `reconcile_ospA_calls` combines assembly and reads; `validate_ospA_inputs` rejects a lone `read1` or `read2`. Enforces the per-IST identity thresholds that the paper describes and LISTT's code never applies, reports ISTs above `max_trusted_ist` (default 17, those in the paper) as `Unknown`, survives LISTT's crash on variant 68, and records ospA locus QC
- `docker/listt` - LISTT 1.2.0 (commit `b9c6569`) with its pinned tool versions, and `audit_listt_tables.py`, which checks LISTT's reference tables and re-derives its IST thresholds
- `tests/run_unit_tests.sh` - `call_ospA_type` section: every flag, the IST cut-off, variant 68, awkward sample names and a lone `read1`
- `classify_kraken2_batched` - Manifest columns after the ID are prefixed `classify_kraken2_`; README acknowledges the Broad Institute's viral-pipelines
- `classify_kraken2_batched` - `manifest_tsv` output with each sample's BAM and report path and summary statistics, uploadable to a Terra data table to attach outputs to sample rows; `removed_reads` summary column
- `workflows/classify_kraken2_batched` - kraken2 classification of many samples (unaligned BAM, or paired or single-end FASTQ) with the database loaded once per VM rather than once per sample, plus optional keeping or removal of reads by clade, by taxid or name. Emits per-sample reports and BAMs and one summary table
- `tasks/kraken2.wdl` - Batch kraken2 task sharing one memory-mapped database across concurrent samples
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

### Removed
- `workflows/bacterial/`, `workflows/eukaryotic/`, `workflows/viral/` - Empty domain-grouping skeleton directories (each held only a blank `README.md`) from the repo's initial scaffold. Every workflow built so far (`assembly_qc`, `pyseer_gwas`, `fetch_reads_from_sra`) already lives flat under `workflows/<name>/`, so the domain subdirectories were unused digital clutter rather than a structure anything actually depends on

### Changed
- `workflows/dna_methylation_calling` - **Breaking.** Takes each isolate's annotation as required inputs, `annotation_gff3s`, `annotation_faas` and `annotation_fnas`, and no longer runs Bakta, so nothing in it is specific to a species. Removed `run_bakta`, `bakta_db`, `proteins`, `genus`, `species` and `bakta_gff3s`/`bakta_faas`/`bakta_fnas`, and the `bakta_gff3`, `bakta_faa` and `bakta_tsv` outputs, which only echoed inputs or other tools' files. `highlight_genes` and the docs no longer use P. aeruginosa examples. `tasks/bakta.wdl` and the other workflows that call it are unchanged
- `workflows/dna_methylation_calling` - The summary's `n_contigs`, `genome_size` and `n_cds` were always `NA` when annotations were supplied, because they came from the Bakta task's outputs. They are now computed from the assembly and GFF, and `gc_percent` and `mod_codes_seen` are new columns
- `tasks/motif_landscape.wdl` - **Tier 1 scoring no longer depends on genome composition.** A motif with a known modified-base offset (from `modkit find-motifs`, carried through `build_motif_list`) is scored on that base alone, on the strand the motif reads 5'→3' on. Previously every A in an occurrence window counted, including the `N` spacer of a split motif, so a fully methylated bipartite 6mA motif read about 31% in a 38% GC genome and about 44% in a 66% GC genome, against about 95% in both when scored by offset. Motifs with no offset (REBASE-only) keep window scoring. The background now excludes every candidate motif's scored sites. `build_motif_list` and `motif_landscape` outputs gain trailing columns (`offset`; `offset` and `scoring`), so existing column positions are unchanged
- `tasks/motif_landscape_summary.wdl` - The motif table gains a `scoring` column (`offset`, `window` or `mixed`), since percentages scored the two ways are on different scales
- `tasks/panaroo.wdl` - New optional `sample_names` input, so isolate column names in `gene_presence_absence.csv` do not depend on GFF filenames
- `workflows/pyseer_gwas` - Removed the separate annotation of `gene_significant`, which is a subset of `gene_results`. The covariate scan's long-format table is now annotated too
- `tasks/modkit.wdl` - `modkit_find_motifs` is no longer preemptible, and defaults to 32 CPUs
- `tasks/motif_landscape_summary.wdl` - Gene-level ranking breaks CV ties on mean density

### Removed
- `tests/` - Placeholder and workspace-specific input files, and the live-network `fetch_reads_from_sra` tests

### Fixed
- `workflows/dna_methylation_calling` - Panaroo was called without `sample_names`, so its columns took their names from the GFF filenames and the ortholog join, which looks columns up by sample name, failed after Panaroo had finished whenever a supplied annotation's filenames differed from the sample names. The annotation arrays' lengths were also unchecked, surfacing as an index error inside the scatter after alignment had run
- `classify_kraken2_batched` - Output BAMs and reports were never copied off the VM on Terra, because their paths were read from a file written by the task. Output paths are now declared from the inputs
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
