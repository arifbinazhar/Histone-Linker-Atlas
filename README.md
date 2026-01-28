
# HistoneLinker-Atlas

A modular, reproducible data-mining pipeline for identifying disease-associated dysregulation and mutation patterns in linker histone (H1 family) genes across neurodevelopmental disorders and cancer.

## Biological Motivation
Linker histones (H1 variants) play a critical role in higher-order chromatin organization and transcriptional regulation. Increasing evidence suggests variant-specific dysregulation and mutation of H1 proteins in developmental disorders and malignancies, yet a systematic, cross-disease computational framework for prioritizing candidate histones for experimental validation is lacking.

This project addresses this gap by integrating public gene expression and mutation datasets into a unified, disease-centric discovery pipeline.

## Project Scope
**Diseases:**
- Autism Spectrum Disorder (GEO)
- Glioblastoma (TCGA / cBioPortal)

**Target Genes:**
- H1F0, HIST1H1A, HIST1H1B, HIST1H1C, HIST1H1D, HIST1H1E, H1FX

## Pipeline Overview
Disease → Dataset Fetching → Gene Filtering → Expression Analysis → Mutation Analysis → Evidence Integration → Report Generation

## Installation
```bash
git clone https://github.com/yourusername/HistoneLinker-Atlas.git
cd HistoneLinker-Atlas
pip install -r requirements.txt


Usage
python run_pipeline.py --disease autism

Outputs
Normalized expression matrices for H1 family genes
Mutation frequency tables and variant summaries
Ranked disease–histone association tables
Publication-ready figures and reports
Intended Use
This pipeline is designed as an internal research support tool to assist experimental chromatin biology labs in prioritizing linker histone variants and disease contexts for downstream functional validation.


