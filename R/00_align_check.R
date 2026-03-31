#!/usr/bin/env Rscript
# ============================================================
# 00_align_check.R
# Count matrix と Metadata のサンプル整合性チェック
# ============================================================

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg))
  script_dir <- dirname(script_path)
} else {
  script_dir <- getwd()
}

load_path <- file.path(script_dir, "00_load_data.R")
if (!file.exists(load_path)) {
  stop("[ERROR] File not found: ", load_path)
}
source(load_path)

suppressPackageStartupMessages({
  library(getopt)
})

spec <- matrix(c(
  "counts",   "c", 1, "character",
  "metadata", "m", 1, "character",
  "outdir",   "o", 1, "character",
  "verbose",  "v", 0, "logical",
  "help",     "h", 0, "logical"
), byrow = TRUE, ncol = 4)

print_help <- function() {
  cat("RNA-seq sample alignment check\n")
  cat("EN: Column names of the count matrix must match sample IDs in metadata.\n")
  cat("JP: カウント行列の列名はメタデータのサンプルIDと一致している必要があります。\n\n")
  cat("Usage:\n")
  cat("  Rscript 00_align_check.R -c counts.rds -m metadata.xlsx [-o outdir] [-v]\n\n")
  cat("Options:\n")
  cat("  -c, --counts    Path to count matrix (.rds)\n")
  cat("  -m, --metadata  Path to metadata (.xlsx/.csv/.tsv)\n")
  cat("  -o, --outdir    Output directory (optional)\n")
  cat("  -v, --verbose   Verbose logging\n")
  cat("  -h, --help      Show this help message\n")
}

opt <- getopt(spec)

if (!is.null(opt$help)) {
  print_help()
  quit(status = 0)
}

if (is.null(opt$counts) || is.null(opt$metadata)) {
  print_help()
  stop("[ERROR] --counts と --metadata は必須です。", call. = FALSE)
}

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)

if (!is.null(opt$verbose)) {
  message("[DEBUG] script_dir: ", script_dir)
  message("[DEBUG] counts: ", opt$counts)
  message("[DEBUG] metadata: ", opt$metadata)
}

tryCatch({
  align_data(counts_raw, meta_raw)
  message("[INFO] サンプル整合性チェック完了")
}, error = function(e) {
  message("[ERROR] サンプル整合性チェックに失敗しました: ", e$message)
  stop(e)
})
