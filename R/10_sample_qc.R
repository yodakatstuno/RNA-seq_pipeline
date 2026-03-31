#!/usr/bin/env Rscript
# ============================================================
# 10_sample_qc.R
# サンプル品質評価
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
  library(reshape2)
})

option_list <- c(base_option_list, list(
  make_option("--color_by", type = "character", default = "condition")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("10_qc", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata

color_col <- if (opt$color_by %in% colnames(meta)) opt$color_by else colnames(meta)[1]

# ---- Library size ----
lib_sizes <- colSums(counts)
ls_df <- data.frame(
  sample = names(lib_sizes),
  library_size = lib_sizes,
  group = meta[[color_col]]
)
ls_df <- ls_df[order(ls_df$library_size), ]
ls_df$sample <- factor(ls_df$sample, levels = ls_df$sample)

p_lib <- ggplot(ls_df, aes(x = sample, y = library_size / 1e6, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Library Size per Sample", x = "", y = "Library Size (millions)")
ggsave(file.path(out_sub, "library_size.pdf"), p_lib,
       width = 8, height = max(6, nrow(ls_df) * 0.25))
ggsave(file.path(out_sub, "library_size.png"), p_lib,
       width = 8, height = max(6, nrow(ls_df) * 0.25), dpi = 300)

# ---- 検出遺伝子数 ----
n_detected <- colSums(counts > 0)
det_df <- data.frame(
  sample = names(n_detected),
  n_genes = n_detected,
  group = meta[[color_col]]
)
det_df <- det_df[order(det_df$n_genes), ]
det_df$sample <- factor(det_df$sample, levels = det_df$sample)

p_det <- ggplot(det_df, aes(x = sample, y = n_genes, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Detected Genes per Sample", x = "", y = "Number of Detected Genes (count > 0)")
ggsave(file.path(out_sub, "detected_genes.pdf"), p_det,
       width = 8, height = max(6, nrow(det_df) * 0.25))

# ---- 発現分布 (box plot) ----
log_counts <- log2(counts + 1)
log_long <- melt(log_counts)
colnames(log_long) <- c("gene", "sample", "log2_expression")
log_long$group <- meta[as.character(log_long$sample), color_col]

p_box <- ggplot(log_long, aes(x = sample, y = log2_expression, fill = group)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(title = "Expression Distribution per Sample", x = "", y = "log2(count + 1)")
ggsave(file.path(out_sub, "expression_distribution.pdf"), p_box,
       width = max(10, ncol(counts) * 0.4), height = 7)

# ---- 密度プロット ----
p_dens <- ggplot(log_long, aes(x = log2_expression, color = sample)) +
  geom_density() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Expression Density per Sample", x = "log2(count + 1)", y = "Density")
ggsave(file.path(out_sub, "expression_density.pdf"), p_dens, width = 10, height = 6)

# ---- サンプル間相関 ----
cor_mat <- cor(log_counts, method = "spearman")

suppressPackageStartupMessages(library(pheatmap))
anno_col <- data.frame(row.names = colnames(log_counts))
anno_col[[color_col]] <- meta[[color_col]]

pdf(file.path(out_sub, "sample_correlation_heatmap.pdf"), width = 10, height = 10)
pheatmap(cor_mat,
         annotation_col = anno_col,
         annotation_row = anno_col,
         clustering_method = "ward.D2",
         main = "Sample-Sample Spearman Correlation")
dev.off()

# ---- Outlier 検出 ----
mean_cor <- apply(cor_mat, 1, function(x) mean(x[x < 1]))  # 自身を除く
outlier_threshold <- mean(mean_cor) - 2 * sd(mean_cor)
outliers <- names(mean_cor[mean_cor < outlier_threshold])

# QC サマリー
sink(file.path(out_sub, "qc_summary.txt"))
cat("========================================\n")
cat("サンプル品質評価サマリー\n")
cat("========================================\n")
cat("サンプル数:", ncol(counts), "\n")
cat("遺伝子数:", nrow(counts), "\n")
cat("\nLibrary Size:\n")
cat("  Min:", min(lib_sizes), "\n")
cat("  Max:", max(lib_sizes), "\n")
cat("  Mean:", round(mean(lib_sizes)), "\n")
cat("  Median:", round(median(lib_sizes)), "\n")
cat("\n検出遺伝子数:\n")
cat("  Min:", min(n_detected), "\n")
cat("  Max:", max(n_detected), "\n")
cat("  Mean:", round(mean(n_detected)), "\n")
cat("\nサンプル間相関 (Spearman mean):\n")
cat("  Min:", round(min(mean_cor), 4), "\n")
cat("  Max:", round(max(mean_cor), 4), "\n")
cat("\nOutlier候補 (mean cor < mean - 2*SD):\n")
if (length(outliers) > 0) {
  cat("  ", paste(outliers, collapse = ", "), "\n")
} else {
  cat("  なし\n")
}
sink()

write.csv(data.frame(
  sample = names(lib_sizes),
  library_size = lib_sizes,
  detected_genes = n_detected[names(lib_sizes)],
  mean_correlation = mean_cor[names(lib_sizes)],
  potential_outlier = names(lib_sizes) %in% outliers,
  group = meta[[color_col]]
), file.path(out_sub, "qc_metrics.csv"), row.names = FALSE)

message("[完了] 10_sample_qc.R")
