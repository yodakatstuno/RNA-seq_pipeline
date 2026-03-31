#!/usr/bin/env Rscript
# ============================================================
# 11_visualization.R
# 総合的な可視化
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
  library(pheatmap)
  library(RColorBrewer)
})

option_list <- c(base_option_list, list(
  make_option("--de_result", type = "character",
              default = file.path(pipeline_output_root, "01_DE", "de_results.csv")),
  make_option("--plots",     type = "character", default = "heatmap,volcano,ma,pca"),
  make_option("--top_n",     type = "integer",   default = 50),
  make_option("--color_by",  type = "character", default = "condition")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("11_visualization", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

log_counts <- log2(counts + 1)
plots <- tolower(trimws(unlist(strsplit(opt$plots, ","))))
color_col <- if (opt$color_by %in% colnames(meta)) opt$color_by else colnames(meta)[1]

# DE 結果読み込み
de <- NULL
if (file.exists(opt$de_result)) {
  de <- read.csv(opt$de_result, stringsAsFactors = FALSE)
}

# ---- Volcano ----
if ("volcano" %in% plots && !is.null(de)) {
  message("[VIS] Volcano plot...")
  de$neg_log10p <- -log10(de$pvalue)
  de$neg_log10p[is.na(de$neg_log10p)] <- 0

  de$label <- ""
  top_up <- head(de[de$significant == TRUE & de$log2FoldChange > 0, ], 10)
  top_down <- head(de[de$significant == TRUE & de$log2FoldChange < 0, ], 10)
  label_genes <- c(top_up$gene, top_down$gene)
  de$label[de$gene %in% label_genes] <- de$gene[de$gene %in% label_genes]

  p <- ggplot(de, aes(x = log2FoldChange, y = neg_log10p)) +
    geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("grey60", "firebrick")) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
    theme_minimal(base_size = 14) +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)")
  ggsave(file.path(out_sub, "volcano_plot.pdf"), p, width = 10, height = 8)
  ggsave(file.path(out_sub, "volcano_plot.png"), p, width = 10, height = 8, dpi = 300)
}

# ---- MA ----
if ("ma" %in% plots && !is.null(de) && "baseMean" %in% colnames(de)) {
  message("[VIS] MA plot...")
  p <- ggplot(de, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("grey60", "firebrick")) +
    geom_hline(yintercept = 0, linetype = "solid") +
    theme_minimal(base_size = 14) +
    labs(title = "MA Plot", x = "log10(Mean Expression)", y = "log2 Fold Change")
  ggsave(file.path(out_sub, "ma_plot.pdf"), p, width = 10, height = 8)
  ggsave(file.path(out_sub, "ma_plot.png"), p, width = 10, height = 8, dpi = 300)
}

# ---- Heatmap ----
if ("heatmap" %in% plots) {
  message("[VIS] Heatmap...")
  if (!is.null(de) && "significant" %in% colnames(de)) {
    sig <- de[de$significant == TRUE, ]
    top_genes <- head(sig[order(sig$padj), "gene"], opt$top_n)
  } else {
    vars <- apply(log_counts, 1, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:min(opt$top_n, length(vars))]
  }
  top_genes <- top_genes[top_genes %in% rownames(log_counts)]

  if (length(top_genes) > 1) {
    anno_col <- data.frame(row.names = colnames(log_counts))
    for (cn in colnames(meta)) {
      if (is.character(meta[[cn]]) || is.factor(meta[[cn]])) {
        anno_col[[cn]] <- meta[[cn]]
      }
    }

    pdf(file.path(out_sub, "heatmap.pdf"), width = 12, height = max(10, length(top_genes) * 0.2))
    pheatmap(log_counts[top_genes, , drop = FALSE],
             scale = "row",
             annotation_col = anno_col,
             clustering_method = "ward.D2",
             labels_col = sample_labels[colnames(log_counts)],
             fontsize_col = sample_label_fontsize(length(sample_labels)),
             show_rownames = (length(top_genes) <= 80),
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             main = paste("Top", length(top_genes), "Genes Heatmap"))
    dev.off()
  }
}

# ---- PCA ----
if ("pca" %in% plots) {
  message("[VIS] PCA plot...")
  filtered <- select_top_variable(log_counts, 2000)
  pca <- prcomp(t(filtered), scale. = TRUE)
  var_exp <- summary(pca)$importance[2, ] * 100

  pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    group = meta[[color_col]],
    sample = rownames(meta),
    label = sample_labels[rownames(meta)]
  )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = label), size = sample_label_size(length(sample_labels)),
                    max.overlaps = Inf) +
    theme_minimal(base_size = 14) +
    labs(title = "PCA Plot",
         x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         color = color_col)
  ggsave(file.path(out_sub, "pca_plot.pdf"), p, width = 10, height = 8)
  ggsave(file.path(out_sub, "pca_plot.png"), p, width = 10, height = 8, dpi = 300)
}

message("[完了] 11_visualization.R")
