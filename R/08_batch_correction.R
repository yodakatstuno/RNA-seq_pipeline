#!/usr/bin/env Rscript
# ============================================================
# 08_batch_correction.R
# バッチ効果補正: ComBat / limma removeBatchEffect
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
  library(ggrepel)
  library(sva)
  library(limma)
})

option_list <- c(base_option_list, list(
  make_option("--batch_col",     type = "character", default = "batch"),
  make_option("--method",        type = "character", default = "both"),
  make_option("--condition_col", type = "character", default = "condition")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("08_batch_correction", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

if (!(opt$batch_col %in% colnames(meta))) {
  stop("[ERROR] バッチカラム '", opt$batch_col, "' がmetadataに存在しません。")
}

log_counts <- log2(counts + 1)
gene_var <- apply(log_counts, 1, var)
log_counts <- log_counts[gene_var > 0, , drop = FALSE]
if (nrow(log_counts) < 2) {
  stop("[ERROR] 変動のある遺伝子が不足しています。PCA を実行できません。")
}
batch <- meta[[opt$batch_col]]

# ---- PCA before correction ----
pca_before <- prcomp(t(log_counts), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_before$x[, 1], PC2 = pca_before$x[, 2],
  batch = batch,
  condition = if (opt$condition_col %in% colnames(meta)) meta[[opt$condition_col]] else "NA",
  sample = rownames(meta),
  label = sample_labels[rownames(meta)]
)
p_before <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch, shape = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label), size = sample_label_size(length(sample_labels)),
                  max.overlaps = Inf, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCA Before Batch Correction")
ggsave(file.path(out_sub, "pca_before_correction.pdf"), p_before, width = 10, height = 8)
ggsave(file.path(out_sub, "pca_before_correction.png"), p_before, width = 10, height = 8, dpi = 300)

methods <- if (opt$method == "both") c("combat", "limma") else opt$method

# ---- ComBat ----
if ("combat" %in% methods) {
  message("[INFO] ComBat 補正...")
  mod <- NULL
  if (opt$condition_col %in% colnames(meta)) {
    mod <- model.matrix(~ factor(meta[[opt$condition_col]]))
  }
  combat_corrected <- ComBat(dat = log_counts, batch = batch, mod = mod)

  saveRDS(combat_corrected, file.path(out_sub, "corrected_combat.rds"))

  pca_combat <- prcomp(t(combat_corrected), scale. = TRUE)
  pca_df2 <- data.frame(
    PC1 = pca_combat$x[, 1], PC2 = pca_combat$x[, 2],
    batch = batch,
    condition = pca_df$condition,
    sample = rownames(meta),
    label = sample_labels[rownames(meta)]
  )
  p_combat <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = batch, shape = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = label), size = sample_label_size(length(sample_labels)),
                    max.overlaps = Inf, show.legend = FALSE) +
    theme_minimal() +
    labs(title = "PCA After ComBat Correction")
  ggsave(file.path(out_sub, "pca_after_combat.pdf"), p_combat, width = 10, height = 8)
  ggsave(file.path(out_sub, "pca_after_combat.png"), p_combat, width = 10, height = 8, dpi = 300)
}

# ---- limma removeBatchEffect ----
if ("limma" %in% methods) {
  message("[INFO] limma removeBatchEffect 補正...")
  design_bio <- NULL
  if (opt$condition_col %in% colnames(meta)) {
    design_bio <- model.matrix(~ factor(meta[[opt$condition_col]]))
  }
  limma_corrected <- removeBatchEffect(log_counts, batch = batch, design = design_bio)

  saveRDS(limma_corrected, file.path(out_sub, "corrected_limma.rds"))

  pca_limma <- prcomp(t(limma_corrected), scale. = TRUE)
  pca_df3 <- data.frame(
    PC1 = pca_limma$x[, 1], PC2 = pca_limma$x[, 2],
    batch = batch,
    condition = pca_df$condition,
    sample = rownames(meta),
    label = sample_labels[rownames(meta)]
  )
  p_limma <- ggplot(pca_df3, aes(x = PC1, y = PC2, color = batch, shape = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = label), size = sample_label_size(length(sample_labels)),
                    max.overlaps = Inf, show.legend = FALSE) +
    theme_minimal() +
    labs(title = "PCA After limma Batch Correction")
  ggsave(file.path(out_sub, "pca_after_limma.pdf"), p_limma, width = 10, height = 8)
  ggsave(file.path(out_sub, "pca_after_limma.png"), p_limma, width = 10, height = 8, dpi = 300)
}

message("[完了] 08_batch_correction.R")
