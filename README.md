
# Integrated analysis of linker histone dysregulation and genomic alterations across neurodevelopmental disorders (autism) and cancer

This project implements a reproducible bioinformatics pipeline to systematically map **differential expression** and **genomic alterations** of linker histone (H1 family) genes across disease contexts using primary transcriptomic (*GEO*) and cancer genomics (*cBioPortal*) datasets. The pipeline enables automated multi-omics integration to identify linker histones showing convergent evidence of dysregulation or alteration, providing a scalable data-mining framework to support hypothesis generation and experimental prioritization in chromatin biology and disease research.

## Scientific Motivation

Linker histones (H1 family) play a critical role in chromatin organization and gene regulation.
Recent studies suggest their dysregulation may contribute to neurodevelopmental disorders and cancer.

However, systematic analysis across disease datasets is lacking.

This pipeline addresses this gap by:

• Integrating transcriptomic and mutation datasets from GEO and TCGA/cBioPortal  
• Identifying linker histones with consistent dysregulation  
• Prioritizing candidates for experimental validation  

This enables hypothesis generation for chromatin-driven disease mechanisms.

## Key Features

- Automated ingestion of GEO transcriptomic datasets for phenotype-aware linker histone expression analysis
- Direct integration of cancer data from cBioPortal mutation and copy number alteration data.
- Identification of disease-associated linker histone dysregulation across Autism and Diffuse Large B-Cell Lymphoma
- Integrated candidate ranking framework to prioritize linker histones for downstream experimental investigation
- Modular, reproducible Python pipeline designed for extensibility to additional diseases and datasets



## Project Scope
**Diseases:**
- Autism Spectrum Disorder (GEO)
- Diffuse Large B-Cell Lymphoma (cBioPortal)

**Target Genes:**
- H1F0, HIST1H1A, HIST1H1B, HIST1H1C, HIST1H1D, HIST1H1E, H1FX

## Data Sources

This pipeline integrates publicly available datasets:

• GEO – transcriptomic datasets
• cBioPortal – mutation and CNA data

Data is automatically downloaded via API

## Pipeline Overview
Disease → Dataset Fetching → Gene Filtering → Expression Analysis → Mutation Analysis → Evidence Integration → Report Generation

## Pipeline Architecture

data_fetchers/ → downloads datasets  
expression/ → expression analysis  
mutation/ → mutation analysis  
atac/ → ATAC-seq HPC workflow scaffold for GSE236992  
integration/ → integrates multi-omics data  
report/ → generates final candidate reports  
utils/ → shared functions  


## Installation

Try running on google  colab by copy pasting the command given below
```bash
!git clone https://github.com/arifbinazhar/Histone-Linker-Atlas.git
!cd HistoneLinker-Atlas
!pip install -r /content/Histone-Linker-Atlas/requirements.txt

import os
os.chdir('/content/Histone-Linker-Atlas')
!python run_pipeline.py --disease all

```

The code below allows you to run the pipeline in parts

```bash
Usage
!python run_pipeline.py --disease autism
!python run_pipeline.py --disease dlbc
!python run_pipeline.py --disease atac
!python run_pipeline.py --disease all


Outputs
Normalized expression matrices for H1 family genes
Mutation frequency tables and variant summaries
Ranked disease-histone association tables
Bar plot as well as a heatmap
Intended Use
This pipeline is designed as an internal research support tool to assist experimental chromatin biology labs in prioritizing linker histone variants and disease contexts for downstream functional validation.

```

## Graphical Outputs

Graph showing mutation frequency per histone gene in Diffuse Large B-Cell Lymphoma <br><br>
![alt text](https://github.com/arifbinazhar/Histone-Linker-Atlas/blob/main/plots/dlbc_mutation_barplot.png)


Bar plot showing mean expression of the Histone H1 gene. <br><br>
![alt text](https://github.com/arifbinazhar/Histone-Linker-Atlas/blob/main/plots/histone_expression_barplot.png)

## ATAC-seq analysis for GSE236992

This repository includes a reproducible ATAC-seq workflow scaffold for **GSE236992**, the ATAC-seq subseries from a CHD8 autism-associated chromatin study. GSE236992 profiles chromatin accessibility changes caused by hemizygous deletion of **CHD8** in human induced neurons. CHD8 is a high-confidence autism-associated chromatin regulator, making this dataset useful for connecting autism chromatin dysregulation to linker histone biology.

The ATAC module focuses downstream interpretation on the linker histone H1 genes:

`H1F0`, `HIST1H1A`, `HIST1H1B`, `HIST1H1C`, `HIST1H1D`, `HIST1H1E`, and `H1FX`.

### What runs on HPC

Raw-read processing is intended for an HPC server and is documented in `atac/README.md`. The workflow scripts use:

- SRA Toolkit (`prefetch`, `fasterq-dump`) to download reads
- FastQC and MultiQC for FASTQ quality control
- fastp for adapter trimming
- bowtie2 against hg38 for alignment
- samtools and bedtools for BAM filtering, duplicate removal, mitochondrial read removal, and blacklist filtering
- MACS2 for ATAC peak calling and consensus peak generation
- R scripts for differential accessibility and peak annotation

Run the HPC scripts from the repository root after editing `config/atac_config.yaml` and filling true SRA accessions in `metadata/atac_samples_GSE236992.csv`:

```bash
bash atac/download_sra.sh config/atac_config.yaml
bash atac/qc_fastq.sh config/atac_config.yaml
bash atac/trim_reads.sh config/atac_config.yaml
bash atac/align_atac.sh config/atac_config.yaml
bash atac/filter_bam.sh config/atac_config.yaml
bash atac/call_peaks.sh config/atac_config.yaml
Rscript atac/differential_accessibility.R config/atac_config.yaml
Rscript atac/annotate_peaks.R config/atac_config.yaml
python atac/summarize_atac_results.py --config config/atac_config.yaml
```

### What runs in this repo

The normal Python pipeline only runs lightweight ATAC summarization and biological integration. It does not download FASTQ files, align reads, call peaks, or create bigWigs.

```bash
python run_pipeline.py --disease atac
```

This step expects the HPC workflow to have produced:

- `results/atac/GSE236992/qc/atac_qc_summary.csv`
- `results/atac/GSE236992/differential_accessibility/differential_accessibility_GSE236992.csv`
- `results/atac/GSE236992/annotation/peak_gene_links.csv`

It writes linker histone and chromatin regulator peak links to:

`data/processed/integration/atac_histone_links_GSE236992.csv`

### GitHub-lightweight data policy

GitHub should contain scripts, metadata templates, config files, small CSV summaries, and plots. The following are intentionally excluded from GitHub through `.gitignore`:

- FASTQ / FASTQ.GZ / FQ files
- BAM / BAI / SAM files
- bigWig / bedGraph files
- `data/raw/atac/`
- `results/atac/aligned/`
- `results/atac/bigwig/`

The sample metadata file is currently a template. Replace `FILL_FROM_GEO_SRA_RUN_TABLE` values in `metadata/atac_samples_GSE236992.csv` after downloading the GEO/SRA run table.

## Biological Interpretation 

Based on our analysis we were able to identify that **HIST1H1C (H1-2)** was the highest-priority linker histone variant.
This might suggest its central role in autism-associated chromatin dysregulation. This can be seen in the CSV file linked below.

[Link to the CSV result](https://github.com/arifbinazhar/Histone-Linker-Atlas/blob/main/data/processed/integration/ranked_histones.csv)

## Reproducibility

This pipeline ensures reproducibility through:

• Version-controlled code (GitHub)
• Fixed dependency versions
• Modular architecture
• Deterministic execution

All results can be reproduced using the provided config.yaml.


## Citation

If you use this pipeline, please cite:

Azhar AB. Histone-Linker-Atlas: Multi-Omics Pipeline for Linker Histone Dysregulation Analysis.
GitHub, 2026.



