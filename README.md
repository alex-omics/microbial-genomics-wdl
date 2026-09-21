# microbial-genomics-wdl
WDL workflows for studying the genomes of bacterial and other microbial pathogens: assembly QC, annotation, transcriptomics, GWAS, and more. Designed to be implemented in Terra.

## Workflows
| Workflow | Description |  Status   |
| -------- | ----------- | --------- |
| annotate_assembly | Annotate an assembled genome with optional reference input  | 🔜 Planned |
| [assembly_qc](workflows/assembly_qc) | QUAST + BUSCO + CheckM2 over a set of assemblies, collapsed into one summary table | ✅ Available |
| [dna_methylation_calling](workflows/dna_methylation_calling) | Per-isolate ONT methylation calling (self-mapped), REBASE MTase homology ID, and a three-tier descriptive landscape (enrichment, within-genome heterogeneity, cross-isolate variability) across a pangenome | ✅ Available |
| [fetch_reads_from_sra](workflows/fetch_reads_from_sra) | Fetches a single SRA/ENA/DDBJ run and emits gzip-compressed FASTQ, with the platform and layout it actually extracted | ✅ Available |
| multimodal_gwas | Use multiple statistical approaches to correlate genes or variants to phenotypes | 🚧 In development |
| [pangenome](workflows/pangenome) | Bakta + Panaroo across a panel of assemblies — pangenome construction only, deliberately no alignment or tree | 🚧 In development |
| polish_assembly | Uses short reads to improve a long-read assembly | 🚧 In development |
| [pyseer_gwas](workflows/pyseer_gwas) | Standalone pyseer LMM/LRT association over a gene or pangenome-module Rtab, with lineage effects and deliberate covariate-combination scanning | ✅ Available |
| [rna_seq_counts](workflows/rna_seq_counts) | BWA-MEM alignment and featureCounts of paired-end RNA-seq, each sample aligned to its own isolate's assembly (replicates matched to their parent by name). Runs Bakta and Panaroo itself to put isolates on shared ortholog groups, giving one cross-isolate matrix; each step is skipped if you supply its output. Also runs against a single shared reference | 🚧 In development |

## Repository Structure
```
workflows   # Full workflows built from modular tasks
tasks       # Individual reusable WDL task files
docker      # Dockerfiles for custom Docker images
tests       # Test WDLs and pointers to test data
docs        # Per-workflow documentation
```

## Acknowledgements
Workflows are being developed during my time in the [Anahtar Lab](https://anahtarlab.mgh.harvard.edu) and 
[Lemieux Lab](https://www.lemieuxlab.org/) at Massachusetts General Hospital. 
The general structure of this repository is inspired by [Theiagen Genomics](https://github.com/theiagen/public_health_bioinformatics). 
Where specified, Docker images courtesy of [StaPH-B](https://github.com/StaPH-B) and [Quay.io](https://quay.io).