#!/usr/bin/env Rscript
# ============================================================
# 12_phenotype_association.R
# 疾患・表現型関連解析
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
  library(pheatmap)
})

option_list <- c(base_option_list, list(
  make_option("--phenotype_col", type = "character", default = "condition"),
  make_option("--method",        type = "character", default = "all"),
  make_option("--top_n",         type = "integer",   default = 100)
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("12_phenotype", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata

log_counts <- log2(filter_low_counts(counts) + 1)

phenotype <- meta[[opt$phenotype_col]]
is_numeric_pheno <- is.numeric(phenotype)

methods <- if (opt$method == "all") c("correlation", "association", "signature") else opt$method

# ---- Correlation (continuous phenotype or converted) ----
if ("correlation" %in% methods) {
  message("[INFO] 遺伝子-表現型相関解析...")

  if (!is_numeric_pheno) {
    pheno_numeric <- as.numeric(as.factor(phenotype))
  } else {
    pheno_numeric <- phenotype
  }

  cors <- apply(log_counts, 1, function(x) {
    ct <- cor.test(x, pheno_numeric, method = "spearman")
    c(rho = ct$estimate, p = ct$p.value)
  })
  cor_df <- data.frame(
    gene = rownames(log_counts),
    rho = cors["rho.rho", ],
    p_value = cors["p", ]
  )
  cor_df$fdr <- p.adjust(cor_df$p_value, method = "BH")
  cor_df <- cor_df[order(cor_df$p_value), ]

  write.csv(cor_df, file.path(out_sub, "phenotype_correlation.csv"), row.names = FALSE)

  # Top genes
  top <- head(cor_df, opt$top_n)
  p <- ggplot(top, aes(x = reorder(gene, rho), y = rho)) +
    geom_bar(stat = "identity", aes(fill = rho > 0)) +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top", opt$top_n, "Phenotype-Correlated Genes"),
         x = "", y = "Spearman rho")
  ggsave(file.path(out_sub, "correlation_barplot.pdf"), p,
         width = 8, height = max(8, opt$top_n * 0.15))
}

# ---- Association (group comparison) ----
if ("association" %in% methods && !is_numeric_pheno) {
  message("[INFO] 表現型グループ間の関連検定...")

  groups <- unique(phenotype)
  if (length(groups) == 2) {
    # Wilcoxon rank sum test
    g1 <- rownames(meta)[phenotype == groups[1]]
    g2 <- rownames(meta)[phenotype == groups[2]]

    pvals <- apply(log_counts, 1, function(x) {
      wilcox.test(x[g1], x[g2])$p.value
    })
    assoc_df <- data.frame(gene = names(pvals), p_value = pvals)
    assoc_df$fdr <- p.adjust(assoc_df$p_value, method = "BH")

    # effect size (mean diff)
    assoc_df$mean_diff <- rowMeans(log_counts[, g2]) - rowMeans(log_counts[, g1])
    assoc_df <- assoc_df[order(assoc_df$p_value), ]

    write.csv(assoc_df, file.path(out_sub, "phenotype_association.csv"), row.names = FALSE)
  } else if (length(groups) > 2) {
    # Kruskal-Wallis
    pvals <- apply(log_counts, 1, function(x) {
      kruskal.test(x ~ factor(phenotype))$p.value
    })
    assoc_df <- data.frame(gene = names(pvals), p_value = pvals)
    assoc_df$fdr <- p.adjust(assoc_df$p_value, method = "BH")
    assoc_df <- assoc_df[order(assoc_df$p_value), ]

    write.csv(assoc_df, file.path(out_sub, "phenotype_association_kw.csv"), row.names = FALSE)
  }
}

# ---- Expression signature ----
if ("signature" %in% methods) {
  message("[INFO] Expression signature 抽出...")

  vars <- apply(log_counts, 1, var)
  top_var_genes <- names(sort(vars, decreasing = TRUE))[1:min(opt$top_n, length(vars))]

  anno_col <- data.frame(row.names = colnames(log_counts))
  anno_col[[opt$phenotype_col]] <- meta[[opt$phenotype_col]]

  pdf(file.path(out_sub, "expression_signature_heatmap.pdf"), width = 12, height = 14)
  pheatmap(log_counts[top_var_genes, , drop = FALSE],
           scale = "row",
           annotation_col = anno_col,
           clustering_method = "ward.D2",
           show_rownames = (length(top_var_genes) <= 80),
           main = paste("Expression Signature (Top", length(top_var_genes), "Variable Genes)"))
  dev.off()

  # sample-level signature score (mean of top genes)
  sig_score <- colMeans(log_counts[top_var_genes, ])
  sig_df <- data.frame(
    sample = names(sig_score),
    signature_score = sig_score,
    phenotype = meta[[opt$phenotype_col]]
  )
  write.csv(sig_df, file.path(out_sub, "signature_scores.csv"), row.names = FALSE)

  p <- ggplot(sig_df, aes(x = phenotype, y = signature_score, fill = phenotype)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2) +
    theme_minimal() +
    labs(title = "Expression Signature Score", x = opt$phenotype_col, y = "Signature Score")
  ggsave(file.path(out_sub, "signature_score_boxplot.pdf"), p, width = 8, height = 6)
}

message("[完了] 12_phenotype_association.R")
