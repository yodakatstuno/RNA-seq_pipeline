#!/usr/bin/env Rscript
# ============================================================
# 05_expression_pattern.R
# WGCNA 共発現解析
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

Sys.setenv(OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           KMP_DUPLICATE_LIB_OK = "TRUE",
           KMP_INIT_AT_FORK = "FALSE",
           KMP_SHM_DISABLE = "1")
suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
})

if (exists("disableWGCNAThreads")) {
  disableWGCNAThreads()
}

option_list <- c(base_option_list, list(
  make_option("--power_range",     type = "character", default = "1:20"),
  make_option("--min_module_size", type = "integer",   default = 30),
  make_option("--merge_cut",       type = "double",    default = 0.25),
  make_option("--trait_col",       type = "character", default = "condition")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("05_expression_pattern", opt$outdir)

# ---- データ ----
counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata

log_counts <- log2(counts + 1)
filtered   <- select_top_variable(log_counts, 5000)
datExpr    <- t(filtered)  # samples x genes

# ---- ソフト閾値選択 ----
message("[WGCNA] ソフト閾値を探索中...")
pr <- eval(parse(text = paste0("c(", gsub(":", ",", opt$power_range), ")")))
if (length(pr) == 2 && pr[2] > pr[1]) {
  pr <- seq(pr[1], pr[2])
}

sft <- pickSoftThreshold(datExpr, powerVector = pr, verbose = 3)
power <- sft$powerEstimate
if (is.na(power)) power <- 6
message("[WGCNA] 選択されたソフト閾値: ", power)

# ソフト閾値プロット
pdf(file.path(out_sub, "soft_threshold.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R²)",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = sft$fitIndices[, 1], col = "red")
abline(h = 0.9, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = sft$fitIndices[, 1], col = "red")
dev.off()

# ---- ネットワーク構築 & モジュール検出 ----
message("[WGCNA] ネットワーク構築中...")
net <- blockwiseModules(datExpr,
                         power = power,
                         maxBlockSize = 10000,
                         TOMType = "unsigned",
                         minModuleSize = opt$min_module_size,
                         reassignThreshold = 0,
                         mergeCutHeight = opt$merge_cut,
                         numericLabels = TRUE,
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
message("[WGCNA] モジュール数: ", length(unique(moduleColors)))

# モジュール割り当て保存
module_df <- data.frame(gene = colnames(datExpr),
                         module_number = moduleLabels,
                         module_color = moduleColors)
write.csv(module_df, file.path(out_sub, "module_assignments.csv"), row.names = FALSE)

# デンドログラム
pdf(file.path(out_sub, "gene_dendrogram.pdf"), width = 12, height = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ---- Module-trait 関連 ----
if (opt$trait_col %in% colnames(meta)) {
  message("[WGCNA] Module-trait 関連解析...")

  trait <- meta[[opt$trait_col]]
  if (is.character(trait) || is.factor(trait)) {
    trait_numeric <- as.numeric(as.factor(trait))
  } else {
    trait_numeric <- as.numeric(trait)
  }

  MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  cor_mt <- cor(MEs, trait_numeric, use = "p")
  pval_mt <- corPvalueStudent(cor_mt, nrow(datExpr))

  mt_df <- data.frame(
    module = rownames(cor_mt),
    correlation = as.numeric(cor_mt),
    p_value = as.numeric(pval_mt)
  )
  write.csv(mt_df, file.path(out_sub, "module_trait_correlation.csv"), row.names = FALSE)

  # ヒートマップ
  pdf(file.path(out_sub, "module_trait_heatmap.pdf"), width = 6, height = 10)
  textMatrix <- paste0(signif(cor_mt, 2), "\n(", signif(pval_mt, 1), ")")
  dim(textMatrix) <- dim(cor_mt)
  par(mar = c(6, 8, 3, 3))
  labeledHeatmap(Matrix = cor_mt,
                 xLabels = opt$trait_col,
                 yLabels = rownames(cor_mt),
                 ySymbols = rownames(cor_mt),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1, 1),
                 main = "Module-Trait Relationships")
  dev.off()
}

message("[完了] 05_expression_pattern.R")
