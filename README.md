# microbial-genomics-wdl
WDL workflows for studying the genomes of bacterial and other microbial pathogens: assembly QC, annotation, transcriptomics, GWAS, and more. Designed to be implemented in Terra.

## Workflows
| Workflow | Description |  Status   |
| -------- | ----------- | --------- |
| annotate_assembly | Annotate an assembled genome with optional reference input  | 🔜 Planned |
| assembly_qc | Assess the completeness and quality of an assembled genome | 🚧 In development |
| fetch_sra_reads | Pulls fastq.gz files from Sequence Read Archive into Terra | 🚧 In development |
| [multimodal_gwas](docs/workflows/multimodal_gwas.md) | Use multiple statistical approaches to correlate genes or variants to phenotypes | 🚧 In development |
| polish_assembly | Uses short reads to improve a long-read assembly | 🚧 In development |
| rna_seq_counts | Align trimmed reads using a specified aligner to a reference, and count read pileup | 🔜 Planned |

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

`multimodal_gwas` is a WDL translation of [microGWAS](https://github.com/microbial-pangenomes-lab/microGWAS) 
by the Microbial Pangenomes Lab, MIT-licensed, Copyright (c) 2022, Marco Galardini. 
If you use it, please cite Burgaya et al. 2025, *microGWAS: a computational pipeline to perform 
large-scale bacterial genome-wide association studies*, Microbial Genomics 11:001349. 
See [docs/workflows/multimodal_gwas.md](docs/workflows/multimodal_gwas.md) for full attribution.