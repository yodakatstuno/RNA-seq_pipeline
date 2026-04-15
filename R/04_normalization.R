#!/usr/bin/env Rscript
# ============================================================
# 04_normalization.R
# Expression Normalization / 発現量の正規化
# ============================================================
#
# PURPOSE:
#   Normalizes raw counts using DESeq2 size factors, TMM (edgeR),
#   log2(counts+1), and optional advanced transforms such as
#   VST / rlog / CPM / row-wise Z-score.
#
# INPUTS:  --counts, --metadata (standard)
# PARAMS:  --method (legacy: deseq2|tmm|log2|all; extended: comma-separated)
# OUTPUTS: normalized_*.rds/csv, scaled_zscore.rds/csv, normalization_comparison.pdf
#
# NOTE:
#   The legacy meaning of --method all is preserved intentionally:
#   deseq2 + tmm + log2 only.
#
# ============================================================

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg))
  script_dir <- dirname(script_path)
} else {
  script_dir <- getwd()
}
message("[DEBUG] script_dir: ", script_dir)
load_path <- file.path(script_dir, "00_load_data.R")
if (!file.exists(load_path)) {
  stop("[ERROR] File not found: ", load_path)
}
source(load_path)

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

option_list <- c(base_option_list, list(
  make_option("--method", type = "character", default = "all"),
  make_option("--force_rlog", type = "character", default = "FALSE"),
  make_option("--rlog_max_samples", type = "integer", default = 30),
  make_option("--zscore_source", type = "character", default = "vst")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("04_normalization", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

parse_methods <- function(method_arg) {
  tokens <- trimws(unlist(strsplit(tolower(method_arg), ",")))
  tokens <- tokens[tokens != ""]
  if (length(tokens) == 0) {
    tokens <- "all"
  }

  expanded <- character(0)
  for (token in tokens) {
    if (token == "all") {
      expanded <- c(expanded, "deseq2", "tmm", "log2")
    } else {
      expanded <- c(expanded, token)
    }
  }

  allowed <- c("deseq2", "tmm", "log2", "vst", "rlog", "cpm", "zscore")
  invalid <- setdiff(unique(expanded), allowed)
  if (length(invalid) > 0) {
    stop("[ERROR] 未対応の method: ", paste(invalid, collapse = ", "))
  }
  unique(expanded)
}

save_matrix_outputs <- function(mat, stem) {
  saveRDS(mat, file.path(out_sub, paste0(stem, ".rds")))
  write.csv(mat, file.path(out_sub, paste0(stem, ".csv")))
}

plot_ready_label <- function(method) {
  switch(
    method,
    "deseq2" = "DESeq2 (log2)",
    "tmm" = "TMM CPM (log2)",
    "log2" = "log2(counts + 1)",
    "vst" = "VST",
    "rlog" = "rlog",
    "cpm" = "CPM (log2)",
    NULL
  )
}

results <- list()
methods <- parse_methods(opt$method)
normalization_methods <- intersect(methods, c("deseq2", "tmm", "log2", "vst", "rlog", "cpm"))
downstream_transforms <- intersect(methods, c("zscore"))
force_rlog <- parse_bool(opt$force_rlog)
zscore_source <- tolower(trimws(opt$zscore_source))

message("[INFO] 実行する正規化手法: ", paste(normalization_methods, collapse = ", "))
if (length(downstream_transforms) > 0) {
  message("[INFO] 下流変換として実行: ", paste(downstream_transforms, collapse = ", "))
}

for (method in normalization_methods) {
  message("[INFO] 正規化実行: ", method)
  result <- normalize_counts_matrix(
    counts = counts,
    meta = meta,
    method = method,
    force_rlog = force_rlog,
    rlog_max_samples = opt$rlog_max_samples
  )
  results[[method]] <- result$matrix

  stem <- switch(
    method,
    "deseq2" = "normalized_deseq2",
    "tmm" = "normalized_tmm",
    "log2" = "normalized_log2",
    "vst" = "normalized_vst",
    "rlog" = "normalized_rlog",
    "cpm" = "normalized_cpm"
  )
  save_matrix_outputs(result$matrix, stem)

  if (method == "deseq2" && !is.null(result$extras$size_factors)) {
    sf <- result$extras$size_factors
    write.csv(
      data.frame(sample = names(sf), size_factor = sf),
      file.path(out_sub, "size_factors.csv"),
      row.names = FALSE
    )
  }

  if (method %in% c("tmm", "cpm") && !is.null(result$extras$norm_factors)) {
    write.csv(result$extras$norm_factors,
              file.path(out_sub, "tmm_norm_factors.csv"))
  }
}

if ("zscore" %in% methods) {
  allowed_sources <- c("deseq2", "tmm", "log2", "vst", "rlog", "cpm")
  if (!(zscore_source %in% allowed_sources)) {
    stop("[ERROR] zscore_source は ", paste(allowed_sources, collapse = ", "),
         " のいずれかです。")
  }

  if (is.null(results[[zscore_source]])) {
    message("[INFO] Z-score 用の元データを作成: ", zscore_source)
    source_result <- normalize_counts_matrix(
      counts = counts,
      meta = meta,
      method = zscore_source,
      force_rlog = force_rlog,
      rlog_max_samples = opt$rlog_max_samples
    )
    results[[zscore_source]] <- source_result$matrix
  }

  message("[INFO] 行方向 Z-score 変換を下流処理として実行: source=", zscore_source)
  zscore_mat <- row_zscore_matrix(results[[zscore_source]])
  save_matrix_outputs(zscore_mat, "scaled_zscore")
  writeLines(
    paste0("source_method=", zscore_source),
    con = file.path(out_sub, "scaled_zscore_source.txt")
  )
  results[["zscore"]] <- zscore_mat
}

message("[INFO] 分布比較プロット作成...")
raw_long <- data.frame(
  sample = rep(colnames(counts), each = nrow(counts)),
  value = as.vector(log2(counts + 1)),
  type = "Raw (log2)",
  stringsAsFactors = FALSE
)
raw_long$sample_label <- factor(
  sample_labels[raw_long$sample],
  levels = sample_labels[colnames(counts)]
)

plot_list <- list(raw_long)
for (method in names(results)) {
  method_label <- plot_ready_label(method)
  if (is.null(method_label)) {
    next
  }
  mat <- results[[method]]
  plot_values <- if (method %in% c("deseq2", "tmm", "cpm")) log2(mat + 1) else mat
  method_long <- data.frame(
    sample = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(plot_values),
    type = method_label,
    stringsAsFactors = FALSE
  )
  method_long$sample_label <- factor(
    sample_labels[method_long$sample],
    levels = sample_labels[colnames(counts)]
  )
  plot_list[[length(plot_list) + 1]] <- method_long
}

all_long <- do.call(rbind, plot_list)

p <- ggplot(all_long, aes(x = sample_label, y = value, fill = type)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ type, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1,
    size = sample_label_fontsize(length(sample_labels))
  )) +
  labs(title = "Normalization Comparison", x = "", y = "Expression")

ggsave(file.path(out_sub, "normalization_comparison.pdf"), p,
       width = max(12, ncol(counts) * 0.5), height = 12)
ggsave(file.path(out_sub, "normalization_comparison.png"), p,
       width = max(12, ncol(counts) * 0.5), height = 12, dpi = 200)

message("[完了] 04_normalization.R")
