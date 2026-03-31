#!/usr/bin/env Rscript
# ============================================================
# 04_normalization.R
# Expression Normalization / 発現量の正規化
# ============================================================
#
# PURPOSE:
#   Normalizes raw counts using DESeq2 size factors, TMM (edgeR),
#   and/or log2(counts+1). Produces comparison box plots.
#
# INPUTS:  --counts, --metadata (standard)
# PARAMS:  --method (deseq2|tmm|log2|all)
# OUTPUTS: normalized_*.rds/csv, normalization_comparison.pdf
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
  make_option("--method", type = "character", default = "all")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("04_normalization", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

methods <- if (opt$method == "all") c("deseq2", "tmm", "log2") else opt$method

# ---- DESeq2 size factor normalization ----
if ("deseq2" %in% methods) {
  suppressPackageStartupMessages(library(DESeq2))
  message("[INFO] DESeq2 正規化...")

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = meta,
                                 design = ~ 1)
  dds <- estimateSizeFactors(dds)
  norm_deseq2 <- counts(dds, normalized = TRUE)

  saveRDS(norm_deseq2, file.path(out_sub, "normalized_deseq2.rds"))
  write.csv(norm_deseq2, file.path(out_sub, "normalized_deseq2.csv"))

  sf <- sizeFactors(dds)
  write.csv(data.frame(sample = names(sf), size_factor = sf),
            file.path(out_sub, "size_factors.csv"), row.names = FALSE)
}

# ---- TMM normalization ----
if ("tmm" %in% methods) {
  suppressPackageStartupMessages(library(edgeR))
  message("[INFO] TMM 正規化...")

  y <- DGEList(counts = counts)
  y <- calcNormFactors(y, method = "TMM")
  norm_tmm <- cpm(y, normalized.lib.sizes = TRUE)

  saveRDS(norm_tmm, file.path(out_sub, "normalized_tmm.rds"))
  write.csv(norm_tmm, file.path(out_sub, "normalized_tmm.csv"))

  nf <- y$samples
  write.csv(nf, file.path(out_sub, "tmm_norm_factors.csv"))
}

# ---- log2 変換 ----
if ("log2" %in% methods) {
  message("[INFO] log2(counts + 1) 変換...")
  log2_counts <- log2(counts + 1)

  saveRDS(log2_counts, file.path(out_sub, "normalized_log2.rds"))
  write.csv(log2_counts, file.path(out_sub, "normalized_log2.csv"))
}

# ---- 分布比較プロット ----
message("[INFO] 分布比較プロット作成...")

# box plot: raw vs normalized
raw_long <- data.frame(
  sample = rep(colnames(counts), each = nrow(counts)),
  value = as.vector(log2(counts + 1)),
  type = "Raw (log2)"
)
raw_long$sample_label <- factor(sample_labels[raw_long$sample], levels = sample_labels[colnames(counts)])

plot_list <- list(raw_long)

if ("deseq2" %in% methods) {
  deseq2_long <- data.frame(
    sample = rep(colnames(norm_deseq2), each = nrow(norm_deseq2)),
    value = as.vector(log2(norm_deseq2 + 1)),
    type = "DESeq2 (log2)"
  )
  deseq2_long$sample_label <- factor(sample_labels[deseq2_long$sample], levels = sample_labels[colnames(counts)])
  plot_list <- c(plot_list, list(deseq2_long))
}

if ("tmm" %in% methods) {
  tmm_long <- data.frame(
    sample = rep(colnames(norm_tmm), each = nrow(norm_tmm)),
    value = as.vector(log2(norm_tmm + 1)),
    type = "TMM CPM (log2)"
  )
  tmm_long$sample_label <- factor(sample_labels[tmm_long$sample], levels = sample_labels[colnames(counts)])
  plot_list <- c(plot_list, list(tmm_long))
}

all_long <- do.call(rbind, plot_list)

p <- ggplot(all_long, aes(x = sample_label, y = value, fill = type)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   size = sample_label_fontsize(length(sample_labels)))) +
  labs(title = "Normalization Comparison", x = "", y = "log2 Expression")
ggsave(file.path(out_sub, "normalization_comparison.pdf"), p,
       width = max(12, ncol(counts) * 0.5), height = 12)
ggsave(file.path(out_sub, "normalization_comparison.png"), p,
       width = max(12, ncol(counts) * 0.5), height = 12, dpi = 200)

message("[完了] 04_normalization.R")
