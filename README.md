
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

Graph showing mutation frequency per histone gene in Diffuse Large B-Cell Lymphoma 
![alt text](https://github.com/arifbinazhar/Histone-Linker-Atlas/blob/main/plots/dlbc_mutation_barplot.png)


Bar plot showing mean expression of the Histone H1 gene.
![alt text](https://github.com/arifbinazhar/Histone-Linker-Atlas/blob/main/plots/histone_expression_barplot.png)

## Biological Interpretation 

Based on our analysis we were able to identify that HIST1H1C (H1-2) was the highest-priority linker histone variant.
This might suggest its central role in autism-associated chromatin dysregulation. 
[Download the dataset CSV file](https://github.com/arifbinazhar/Histone-Linker-Atlas/blob/main/data/processed/integration/ranked_histones.csv)


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



