#!/usr/bin/env Rscript
# ============================================================
# 09_regression.R
# 回帰モデル解析
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
  make_option("--formula_str", type = "character", default = "~ condition"),
  make_option("--gene_list",   type = "character", default = ""),
  make_option("--method",      type = "character", default = "limma")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("09_regression", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata

log_counts <- log2(counts + 1)

# 非構文的な列名に対応するため、モデル用に列名を安全化
name_map <- setNames(make.names(colnames(meta), unique = TRUE), colnames(meta))
meta_model <- meta
colnames(meta_model) <- name_map
formula_str_model <- opt$formula_str
for (cn in names(name_map)) {
  if (cn != name_map[[cn]]) {
    if (grepl(cn, formula_str_model, fixed = TRUE)) {
      formula_str_model <- gsub(cn, name_map[[cn]], formula_str_model, fixed = TRUE)
    }
  }
}
if (formula_str_model != opt$formula_str) {
  message("[INFO] 解析用に式を変換しました: ", opt$formula_str, " -> ", formula_str_model)
}

# 対象遺伝子
if (opt$gene_list != "" && opt$gene_list != "NA") {
  genes <- trimws(unlist(strsplit(opt$gene_list, ",")))
  genes <- intersect(genes, rownames(log_counts))
} else {
  # フィルタリングして全遺伝子
  genes <- rownames(filter_low_counts(counts))
}

message("[INFO] 解析対象遺伝子数: ", length(genes))

# ---- limma ベース ----
if (opt$method == "limma") {
  suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
  })

  formula_obj <- as.formula(formula_str_model)
  design <- model.matrix(formula_obj, data = meta_model)
  message("[INFO] Design matrix: ", paste(colnames(design), collapse = ", "))

  y <- DGEList(counts = counts[genes, , drop = FALSE])
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  # 各係数のresults
  for (i in 2:ncol(design)) {
    coef_name <- colnames(design)[i]
    tt <- topTable(fit, coef = i, number = Inf)
    tt$gene <- rownames(tt)
    write.csv(tt, file.path(out_sub, paste0("regression_", coef_name, ".csv")),
              row.names = FALSE)
    message("[INFO] 係数 '", coef_name, "': ",
            sum(tt$adj.P.Val < 0.05), " genes at FDR < 0.05")
  }
}

# ---- lm / glm ベース (遺伝子ごと) ----
if (opt$method %in% c("lm", "glm")) {
  message("[INFO] 遺伝子ごとの ", opt$method, " 回帰...")

  results_list <- list()
  for (g in genes) {
    df_tmp <- meta_model
    df_tmp$expression <- as.numeric(log_counts[g, ])

    formula_full <- as.formula(paste("expression", formula_str_model))

    if (opt$method == "lm") {
      fit <- lm(formula_full, data = df_tmp)
    } else {
      fit <- glm(formula_full, data = df_tmp, family = gaussian())
    }

    s <- summary(fit)
    coefs <- as.data.frame(s$coefficients)
    coefs$gene <- g
    coefs$term <- rownames(coefs)
    results_list[[g]] <- coefs
  }

  results_all <- do.call(rbind, results_list)
  colnames(results_all) <- c("estimate", "std_error", "t_or_z_value", "p_value", "gene", "term")
  results_all$fdr <- p.adjust(results_all$p_value, method = "BH")

  write.csv(results_all, file.path(out_sub, "regression_results.csv"), row.names = FALSE)

  # volcano-like plot for each term
  terms <- unique(results_all$term)
  terms <- terms[terms != "(Intercept)"]

  for (tm in terms) {
    sub_df <- results_all[results_all$term == tm, ]
    sub_df$significant <- sub_df$fdr < 0.05

    p <- ggplot(sub_df, aes(x = estimate, y = -log10(p_value))) +
      geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
      scale_color_manual(values = c("grey60", "red")) +
      theme_minimal() +
      labs(title = paste("Regression:", tm), x = "Estimate", y = "-log10(p-value)")
    ggsave(file.path(out_sub, paste0("regression_volcano_", gsub("[^a-zA-Z0-9]", "_", tm), ".pdf")),
           p, width = 8, height = 6)
  }
}

message("[完了] 09_regression.R")
