#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[[1]], "config/atac_config.yaml")
config <- yaml::read_yaml(config_file)

da_file <- file.path(config$paths$differential_dir, "differential_accessibility_GSE236992.csv")
out_dir <- config$paths$annotation_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

da <- read.csv(da_file, stringsAsFactors = FALSE)

if (!requireNamespace("ChIPseeker", quietly = TRUE) ||
    !requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  message("ChIPseeker/TxDb/org.Hs.eg.db not installed. Writing annotation template.")
  output <- data.frame(
    peak_id = da$peak_id,
    chr = da$chr,
    start = da$start,
    end = da$end,
    nearest_gene = da$nearest_gene,
    distance_to_tss = NA_integer_,
    annotation = NA_character_
  )
  write.csv(output, file.path(out_dir, "peak_gene_links.csv"), row.names = FALSE)
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

gr <- GenomicRanges::GRanges(
  seqnames = da$chr,
  ranges = IRanges::IRanges(start = da$start, end = da$end),
  peak_id = da$peak_id
)

anno <- annotatePeak(
  gr,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db"
)

anno_df <- as.data.frame(anno)
output <- data.frame(
  peak_id = anno_df$peak_id,
  chr = as.character(anno_df$seqnames),
  start = anno_df$start,
  end = anno_df$end,
  nearest_gene = anno_df$SYMBOL,
  distance_to_tss = anno_df$distanceToTSS,
  annotation = anno_df$annotation
)

write.csv(output, file.path(out_dir, "peak_gene_links.csv"), row.names = FALSE)
