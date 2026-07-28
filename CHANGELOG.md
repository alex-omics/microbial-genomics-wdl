# Changelog

## [Unreleased]
### Added
- `multimodal_gwas` - Terra WDL translation of microGWAS (Burgaya et al. 2025), targeting faithful reproduction of the published *E. coli* virulence GWAS as Milestone 1
- `workflows/assembly_qc/assembly_qc.wdl` - Scatters QUAST, BUSCO, and CheckM2 over a set of assemblies and emits one merged summary TSV
- `tasks/checkm2.wdl` - CheckM2 completeness and contamination, with the DIAMOND database supplied as a `File` input (Terra-native) or downloaded at runtime
- `busco-prokaryota:5.8.0` - Docker image with the complete OrthoDB v10 prokaryote set (99 datasets) plus placement files for offline auto-lineage
- `tasks/busco.wdl` - Reusable BUSCO task with auto-lineage and explicit lineage support, plus an optional lineage tarball override

### Fixed
- `tasks/busco.wdl` - Pointed `--download_path` at the consolidated `/busco_downloads`; the previous path held no manifest, so auto-lineage could not run and only spirochaete datasets resolved
- `tasks/busco.wdl` - `grep` for the lineage-used line had the filename inside the pattern with no file argument, so it read stdin and returned nothing
- `tasks/busco.wdl` - Summary was copied from `<sample>_busco/` while BUSCO wrote to `<sample>/`; auto-lineage also emits several `short_summary*` files, which broke the glob
- `tasks/busco.wdl` - Runtime used `disk:` / `local_disk`; Terra requires `disks:` / `local-disk`, so the attribute was silently ignored
- `tasks/quast.wdl` - Invoked `quast`, which does not exist in `staphb/quast:5.2.0` (the entrypoint is `quast.py`)
- `tasks/quast.wdl` - Quoted domain flag expanded to an empty positional argument for prokaryotes, which QUAST read as an empty input filename

## [0.1.0] — 17 June 2026
- Initial repo structure
