# ATAC-seq workflow for GSE236992

GSE236992 is the ATAC-seq subseries from a CHD8 autism-associated chromatin study. It profiles chromatin accessibility changes after hemizygous CHD8 deletion in human induced neurons. This module connects differential accessibility to linker histone H1 genes: H1F0, HIST1H1A, HIST1H1B, HIST1H1C, HIST1H1D, HIST1H1E, and H1FX.

## Repository policy

GitHub tracks scripts, configs, metadata templates, small summary CSV files, plots, and downstream integration code. Raw reads and large intermediates belong on HPC storage and are excluded by `.gitignore`.

Do not commit FASTQ, BAM, SAM, BAI, bigWig, bedGraph, or large alignment/coverage directories.

## Required HPC tools

- SRA Toolkit: `prefetch`, `fasterq-dump`
- FastQC and MultiQC
- fastp
- bowtie2
- samtools
- bedtools
- MACS2
- R packages for downstream analysis: DiffBind or DESeq2, ChIPseeker or HOMER-compatible annotation tools

## Configuration

Edit `config/atac_config.yaml` before running on HPC:

- Set `reference.bowtie2_index` to the hg38 bowtie2 index prefix.
- Set `reference.blacklist_bed` to the hg38 ENCODE blacklist BED.
- Set optional `reference.tss_bed` and `reference.chrom_sizes` for deeper QC.
- Fill `metadata/atac_samples_GSE236992.csv` with true SRA run accessions from the GEO/SRA run table.

## HPC workflow

Run scripts from the repository root:

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

The HPC workflow writes into `data/raw/atac/GSE236992/` and `results/atac/GSE236992/`.

## QC metrics

`atac/filter_bam.sh` writes samtools flagstat, idxstats, and markdup logs used by `atac/summarize_atac_results.py` to report total reads, mapped reads, duplicate rate, mitochondrial read fraction, final usable reads, estimated library size, and peak counts.

FRiP and fragment size metrics can be added on HPC after peak calling with commands like:

```bash
samtools view -c results/atac/GSE236992/filtered/SAMPLE.filtered.bam
bedtools intersect -u -abam results/atac/GSE236992/filtered/SAMPLE.filtered.bam -b results/atac/GSE236992/peaks/SAMPLE_peaks.narrowPeak | samtools view -c -
samtools view results/atac/GSE236992/filtered/SAMPLE.filtered.bam | awk '($9 > 0 && $9 < 2000) {print $9}' > results/atac/GSE236992/qc/SAMPLE.fragment_sizes.txt
```

TSS enrichment requires a curated hg38 TSS BED file, configured as `reference.tss_bed` in `config/atac_config.yaml`. Use an HPC ATAC QC tool such as deepTools, ataqv, or custom bedtools windows around TSS sites when that reference is available.

## Lightweight repo integration

After differential accessibility and peak annotation are available, run:

```bash
python run_pipeline.py --disease atac
```

This summarizes available ATAC QC files and writes linker histone/chromatin regulator links to:

`data/processed/integration/atac_histone_links_GSE236992.csv`
