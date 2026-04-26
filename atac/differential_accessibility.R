#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[[1]], "config/atac_config.yaml")
config <- yaml::read_yaml(config_file)

metadata <- read.csv(config$paths$metadata, stringsAsFactors = FALSE)
out_dir <- config$paths$differential_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_sheet <- data.frame(
  SampleID = metadata$sample_id,
  Tissue = metadata$cell_type,
  Factor = "ATAC",
  Condition = metadata$condition,
  Replicate = metadata$replicate,
  bamReads = metadata$bam,
  Peaks = metadata$peak_file,
  PeakCaller = "narrow"
)

sample_sheet_file <- file.path(out_dir, "diffbind_sample_sheet_GSE236992.csv")
write.csv(sample_sheet, sample_sheet_file, row.names = FALSE)

condition_counts <- table(sample_sheet$Condition)
if (length(condition_counts) < 2 || any(condition_counts < 2)) {
  message("Need at least two replicates per condition for DiffBind contrast. Wrote sample sheet only: ", sample_sheet_file)
  empty_output <- data.frame(
    peak_id = character(),
    chr = character(),
    start = integer(),
    end = integer(),
    baseMean = numeric(),
    log2FoldChange = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    direction = character(),
    nearest_gene = character()
  )
  write.csv(
    empty_output,
    file.path(out_dir, "differential_accessibility_GSE236992.csv"),
    row.names = FALSE
  )
  quit(status = 0)
}

if (!requireNamespace("DiffBind", quietly = TRUE)) {
  message("DiffBind is not installed. Wrote sample sheet and empty result template.")
  empty_output <- data.frame(
    peak_id = character(),
    chr = character(),
    start = integer(),
    end = integer(),
    baseMean = numeric(),
    log2FoldChange = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    direction = character(),
    nearest_gene = character()
  )
  write.csv(
    empty_output,
    file.path(out_dir, "differential_accessibility_GSE236992.csv"),
    row.names = FALSE
  )
  quit(status = 0)
}

suppressPackageStartupMessages(library(DiffBind))

dba_obj <- dba(sampleSheet = sample_sheet_file)
dba_obj <- dba.count(dba_obj, peaks = file.path(config$paths$peak_dir, "consensus_peaks.bed"))
dba_obj <- dba.contrast(
  dba_obj,
  categories = DBA_CONDITION,
  minMembers = 2
)
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)
report <- as.data.frame(dba.report(dba_obj, method = DBA_DESEQ2, th = 1))

if (nrow(report) == 0) {
  output <- data.frame(
    peak_id = character(),
    chr = character(),
    start = integer(),
    end = integer(),
    baseMean = numeric(),
    log2FoldChange = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    direction = character(),
    nearest_gene = character()
  )
} else {
  output <- data.frame(
    peak_id = paste0(report$seqnames, ":", report$start, "-", report$end),
    chr = as.character(report$seqnames),
    start = report$start,
    end = report$end,
    baseMean = if ("Conc" %in% names(report)) report$Conc else NA_real_,
    log2FoldChange = report$Fold,
    pvalue = report$p.value,
    padj = report$FDR,
    direction = ifelse(report$Fold > 0, "more_accessible_in_CHD8_deletion", "less_accessible_in_CHD8_deletion"),
    nearest_gene = NA_character_
  )
}

write.csv(
  output,
  file.path(out_dir, "differential_accessibility_GSE236992.csv"),
  row.names = FALSE
)
