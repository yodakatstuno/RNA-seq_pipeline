#!/usr/bin/env Rscript
# ============================================================
# 01_differential_expression.R
# Differential Expression Analysis / 差次発現解析
# ============================================================
#
# PURPOSE:
#   Identifies genes with significant expression differences
#   between experimental conditions using DESeq2, edgeR, or limma.
#
# REQUIRED INPUTS:
#   --counts    Count matrix (.rds) — raw integer counts, genes × samples
#   --metadata  Metadata (.xlsx/.csv) — must contain condition column
#
# KEY PARAMETERS:
#   --tool          deseq2 | edger | limma
#   --condition_col Metadata column for group comparison (e.g., "genotype")
#   --ref_level     Reference group (e.g., "WT", "control")
#   --treat_level   Treatment group (e.g., "Ho", "treatment")
#   --design_mode   simple (~ condition)
#                   additive (~ factor1 + factor2)
#                   interaction (~ factor1 + factor2 + factor1:factor2)
#   --factors       Comma-separated factor columns for additive/interaction
#                   Example: "dpf,genotype"
#   --contrast_mode pairwise | auto | manual
#
# DESIGN FORMULA EXAMPLES:
#   Simple:      ~ genotype          (compare WT vs Ho)
#   Additive:    ~ dpf + genotype    (account for dpf effect)
#   Interaction: ~ dpf + genotype + dpf:genotype
#
# OUTPUTS (in <outdir>/<analysis_name>/):
#   results/de_results.csv          — Full DE results table
#   plots/volcano_plot.pdf/png      — Volcano plot
#   plots/ma_plot.pdf/png           — MA plot
#   plots/heatmap_top_deg.pdf/png   — Top DEG heatmap
#   logs/summary.txt                — Analysis summary
#   logs/sessionInfo.txt            — R session info
#
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
})
# NOTE: "package built under different R version" warnings are usually safe to ignore unless functionality breaks.

# ---- 引数パース ----
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

option_list <- c(base_option_list, list(
  make_option("--tool",          type = "character", default = "deseq2"),
  make_option("--condition_col", type = "character", default = "condition"),
  make_option("--ref_level",     type = "character", default = "control"),
  make_option("--treat_level",   type = "character", default = "treatment"),
  make_option("--fdr_cutoff",    type = "double",    default = 0.05),
  make_option("--lfc_cutoff",    type = "double",    default = 1.0),
  make_option("--design_mode",   type = "character", default = "simple"),
  make_option("--factors",       type = "character", default = ""),
  make_option("--contrast_mode", type = "character", default = "pairwise"),
  make_option("--contrast",      type = "character", default = ""),
  make_option("--test",          type = "character", default = "Wald"),
  make_option("--reduced_design", type = "character", default = ""),
  make_option("--subset_col",    type = "character", default = ""),
  make_option("--subset_val",    type = "character", default = ""),
  make_option("--analysis_name", type = "character", default = "01_DE"),
  make_option("--explain_results", type = "character", default = "FALSE")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("01_DE", opt$outdir, opt$analysis_name)
results_dir <- file.path(out_sub, "results")
plots_dir <- file.path(out_sub, "plots")
interp_dir <- file.path(out_sub, "interpretation")
logs_dir <- file.path(out_sub, "logs")
dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  v <- tolower(as.character(x))
  v %in% c("true", "1", "yes", "y")
}

explain_results <- parse_bool(opt$explain_results)

get_factors <- function() {
  if (opt$factors == "" || is.na(opt$factors)) {
    return(opt$condition_col)
  }
  # Remove surrounding quotes that may be passed from command line
  raw <- gsub('^["\']|["\']$', '', opt$factors)
  fs <- trimws(unlist(strsplit(raw, ",")))
  fs <- fs[fs != ""]
  if (length(fs) == 0) opt$condition_col else fs
}

safe_factor_name <- function(x) {
  if (make.names(x) != x || grepl("[^[:alnum:]_.]", x)) {
    return(paste0("`", x, "`"))
  }
  x
}

build_design <- function(factors, condition_col) {
  factors_safe <- vapply(factors, safe_factor_name, character(1))
  if (opt$design_mode == "simple") {
    return(paste("~", safe_factor_name(condition_col)))
  }
  if (opt$design_mode == "additive") {
    return(paste("~", paste(factors_safe, collapse = " + ")))
  }
  if (opt$design_mode == "interaction") {
    if (length(factors) < 2) {
      message("[WARNING] interaction には2因子以上が必要です。additive に変更します。")
      return(paste("~", paste(factors_safe, collapse = " + ")))
    }
    base <- paste(factors_safe, collapse = " + ")
    inter <- paste(factors_safe[1], factors_safe[2], sep = ":")
    return(paste("~", base, "+", inter))
  }
  message("[WARNING] design_mode が不正なため simple に変更します。")
  return(paste("~", safe_factor_name(condition_col)))
}

parse_subset_values <- function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(character(0))
  }
  vals <- trimws(unlist(strsplit(x, ",")))
  vals[vals != ""]
}

apply_subset <- function(counts, meta, subset_col, subset_val) {
  vals <- parse_subset_values(subset_val)
  if (subset_col == "" || length(vals) == 0) {
    return(list(counts = counts, meta = meta, applied = FALSE, values = character(0)))
  }
  if (!(subset_col %in% colnames(meta))) {
    stop("[ERROR] subset_col が metadata に存在しません: ", subset_col)
  }

  meta_vals <- trimws(as.character(meta[[subset_col]]))
  keep <- !is.na(meta_vals) & meta_vals %in% vals
  if (!any(keep)) {
    stop("[ERROR] subset 条件に一致するサンプルがありません: ",
         subset_col, " in {", paste(vals, collapse = ", "), "}")
  }

  meta_sub <- meta[keep, , drop = FALSE]
  counts_sub <- counts[, rownames(meta_sub), drop = FALSE]
  if (ncol(counts_sub) != nrow(meta_sub)) {
    stop("[ERROR] subset 後の counts / metadata 整合に失敗しました。")
  }
  message("[INFO] Subsetting applied: ", subset_col, " in {",
          paste(vals, collapse = ", "), "} -> ", nrow(meta_sub), " samples retained")
  list(counts = counts_sub, meta = meta_sub, applied = TRUE, values = vals)
}

normalize_test_name <- function(x) {
  x_lower <- tolower(trimws(as.character(x)))
  if (x_lower == "wald") return("Wald")
  if (x_lower == "lrt") return("LRT")
  stop("[ERROR] --test は Wald または LRT を指定してください。")
}

# ---- データ読み込み ----
counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
subset_info <- apply_subset(counts, meta, opt$subset_col, opt$subset_val)
counts <- subset_info$counts
meta <- subset_info$meta
opt$test <- normalize_test_name(opt$test)

# 因子の準備
factors <- get_factors()
missing_factors <- setdiff(factors, colnames(meta))
if (length(missing_factors) > 0) {
  stop("[ERROR] Metadata に存在しない因子: ", paste(missing_factors, collapse = ", "))
}

for (f in factors) {
  if (is.numeric(meta[[f]])) {
    message("[WARNING] 数値が条件として使われています。カテゴリ変数として扱います。")
  }
  meta[[f]] <- as.factor(meta[[f]])
}

# 条件カラムを参照レベルで固定（他のレベルは保持）
if (opt$condition_col %in% colnames(meta)) {
  meta[[opt$condition_col]] <- as.factor(meta[[opt$condition_col]])
  if (!(opt$ref_level %in% levels(meta[[opt$condition_col]]))) {
    stop("[ERROR] 参照レベルが条件カラムに存在しません: ", opt$ref_level)
  }
  if (!(opt$treat_level %in% levels(meta[[opt$condition_col]]))) {
    stop("[ERROR] 処理レベルが条件カラムに存在しません: ", opt$treat_level)
  }
  meta[[opt$condition_col]] <- relevel(meta[[opt$condition_col]], ref = opt$ref_level)
}
message("[DEBUG] condition_col: ", opt$condition_col)
if (opt$condition_col %in% colnames(meta)) {
  message("[DEBUG] condition levels: ", paste(levels(meta[[opt$condition_col]]), collapse = ", "))
}

validate_condition <- function(col) {
  if (!(col %in% colnames(meta))) {
    stop("[ERROR] 指定された列 '", col, "' はmetadataに存在しません")
  }
  if (any(is.na(meta[[col]]))) {
    stop("[ERROR] 条件カラムに NA が含まれています: ", col)
  }
  if (any(meta[[col]] == "")) {
    stop("[ERROR] 条件カラムに空文字が含まれています: ", col)
  }
  if (length(unique(meta[[col]])) < 2) {
    stop("[ERROR] 条件のレベルが2未満です: ", col)
  }
  if (any(table(meta[[col]]) < 2)) {
    stop("[ERROR] 各群に少なくとも2サンプル必要です: ", col)
  }
  message("[INFO] 条件: ", col)
  message("[INFO] レベル: ", paste(levels(meta[[col]]), collapse = ", "))
  message("[INFO] サンプル数: ", paste(table(meta[[col]]), collapse = ", "))
}

validate_condition(opt$condition_col)

# モデル用に安全な列名へ変換（DESeq2 が非構文的な列名を扱えないため）
name_map <- setNames(make.names(colnames(meta), unique = TRUE), colnames(meta))
meta_model <- meta
colnames(meta_model) <- name_map
condition_col_model <- name_map[[opt$condition_col]]
factors_model <- unname(name_map[factors])
if (any(is.na(factors_model)) || is.na(condition_col_model)) {
  stop("[ERROR] モデル用の列名変換に失敗しました。")
}
if (any(names(name_map) != name_map)) {
  message("[INFO] モデル用に列名を変換しました: ",
          paste0(names(name_map), "→", name_map, collapse = ", "))
}

design_formula_str <- build_design(factors_model, condition_col_model)
message("[INFO] Design formula: ", design_formula_str)
message("[INFO] Factors: ", paste(factors, collapse = ", "))

# 反復数とバランスのチェック
comb <- interaction(meta[, factors, drop = FALSE], drop = TRUE)
rep_counts <- table(comb)
if (any(rep_counts < 2)) {
  message("[WARNING] 反復数が2未満のグループがあります: ",
          paste(names(rep_counts)[rep_counts < 2], collapse = ", "))
}
if (max(rep_counts) / max(1, min(rep_counts)) > 2) {
  message("[WARNING] デザインが不均衡です（最大/最小の比が大きい）。")
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

annotate_significance <- function(df) {
  df$significant <- (abs(df$log2FoldChange) >= opt$lfc_cutoff &
                      df$padj < opt$fdr_cutoff)
  df$significant[is.na(df$significant)] <- FALSE
  df
}

annotate_interaction <- function(df, term) {
  df$effect_direction <- ifelse(is.na(df$log2FoldChange), NA,
                                ifelse(df$log2FoldChange > 0,
                                       "stronger_in_treatment", "weaker_in_treatment"))
  df$biological_meaning <- paste0(
    "Positive log2FC means stronger response in mutant relative to control under treatment. (", term, ")"
  )
  df$interpretation_note <- "Interaction terms require careful interpretation; main effects are not equivalent."
  df
}

interaction_explanation <- function(term) {
  # Simple bilingual explanation based on term name
  en <- paste0("This coefficient represents how the effect of ",
               term, " differs between the reference and comparison levels.")
  jp <- paste0("この係数は、", term, " の効果が参照レベルと比較レベルでどのように異なるかを示します。")
  list(en = en, jp = jp)
}

write_interpretation_report <- function(path) {
  cat("========================================\n", file = path)
  cat("Interpretation Report\n", file = path, append = TRUE)
  cat("========================================\n", file = path, append = TRUE)
  cat("EN: Positive log2FC means higher expression in the comparison group.\n", file = path, append = TRUE)
  cat("EN: Negative log2FC means lower expression in the comparison group.\n", file = path, append = TRUE)
  cat("EN: Main effect describes the average effect of a factor.\n", file = path, append = TRUE)
  cat("EN: Interaction effect describes how one factor’s effect changes across another factor.\n", file = path, append = TRUE)
  cat("EN: Volcano plots show log2FC (x) vs -log10(p-value) (y).\n\n", file = path, append = TRUE)
  cat("JP: log2FC が正なら比較群で発現が高いことを示します。\n", file = path, append = TRUE)
  cat("JP: log2FC が負なら比較群で発現が低いことを示します。\n", file = path, append = TRUE)
  cat("JP: 主効果は因子の平均的な影響を示します。\n", file = path, append = TRUE)
  cat("JP: 交互作用効果は因子の影響が別の因子によってどのように変わるかを示します。\n", file = path, append = TRUE)
  cat("JP: Volcano plot は log2FC (x) と -log10(p値) (y) を表示します。\n", file = path, append = TRUE)
}

# ---- DESeq2 ----
run_deseq2 <- function() {
  suppressPackageStartupMessages(library(DESeq2))
  message("[DESeq2] 解析開始...")

  formula_obj <- as.formula(design_formula_str)
  # DESeq2 requires integer counts; round if necessary
  counts_int <- round(counts)
  storage.mode(counts_int) <- "integer"
  dds <- DESeqDataSetFromMatrix(countData = counts_int,
                                colData = meta_model,
                                design = formula_obj)
  if (opt$test == "LRT") {
    if (opt$reduced_design == "" || is.na(opt$reduced_design)) {
      stop("[ERROR] --test=LRT では --reduced_design が必要です。")
    }
    reduced_formula <- as.formula(opt$reduced_design)
    message("[INFO] DESeq2 test: LRT")
    message("[INFO] Reduced formula: ", opt$reduced_design)
    dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
  } else {
    if (opt$reduced_design != "") {
      message("[WARNING] --reduced_design は --test=Wald では使用されません。")
    }
    message("[INFO] DESeq2 test: Wald")
    dds <- DESeq(dds, test = "Wald")
  }
  rn <- resultsNames(dds)
  message("[INFO] resultsNames: ", paste(rn, collapse = ", "))

  if (opt$design_mode == "interaction" && opt$contrast_mode == "pairwise") {
    message("[WARNING] You are using an interaction design. Main effects and pairwise comparisons are not equivalent.")
    message("[WARNING] 交互作用モデルでは、単純な群間比較は主効果と同じ意味にはなりません。")
  }

  if (opt$test == "LRT") {
    if (opt$contrast_mode != "pairwise" || opt$contrast != "") {
      message("[WARNING] LRT では contrast_mode / contrast は使用されません。full vs reduced モデルの結果を返します。")
    }
    res <- results(dds, alpha = opt$fdr_cutoff)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- annotate_significance(res_df)
    message("[DESeq2] DEG数 (LRT): ", sum(res_df$significant, na.rm = TRUE))
    return(list(results = res_df, dds = dds))
  }

  if (opt$contrast_mode == "auto") {
    summary_path <- NULL
    if (opt$design_mode == "interaction") {
      summary_path <- file.path(interp_dir, "interaction_summary.txt")
      cat("Design formula: ", design_formula_str, "\n", file = summary_path)
      cat("Factors: ", paste(factors, collapse = ", "), "\n\n", file = summary_path, append = TRUE)
      cat("Interaction terms:\n", file = summary_path, append = TRUE)
    }
    for (nm in rn) {
      res <- results(dds, name = nm, alpha = opt$fdr_cutoff)
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df <- annotate_significance(res_df)
      if (grepl(":", nm)) {
        res_df <- annotate_interaction(res_df, nm)
        if (!is.null(summary_path)) {
          expl <- interaction_explanation(nm)
          cat("- ", nm, "\n", file = summary_path, append = TRUE)
          cat("  EN: ", expl$en, "\n", file = summary_path, append = TRUE)
          cat("  JP: ", expl$jp, "\n", file = summary_path, append = TRUE)
          cat("  Significant genes: ", sum(res_df$significant, na.rm = TRUE), "\n\n",
              file = summary_path, append = TRUE)
        }
      }
      write.csv(res_df, file.path(results_dir, paste0("de_results_", safe_name(nm), ".csv")),
                row.names = FALSE)

      # Volcano plot per term
      p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
        scale_color_manual(values = c("grey60", "red")) +
        geom_hline(yintercept = -log10(opt$fdr_cutoff), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = c(-opt$lfc_cutoff, opt$lfc_cutoff),
                   linetype = "dashed", color = "blue") +
        theme_minimal() +
        labs(title = ifelse(grepl(":", nm),
                            paste("Interaction Effect Volcano Plot:", nm),
                            paste("Main Effect Volcano Plot:", nm)),
             x = "log2 Fold Change", y = "-log10(p-value)")
      ggsave(file.path(plots_dir, paste0("volcano_", safe_name(nm), ".pdf")), p_vol, width = 8, height = 6)
      ggsave(file.path(plots_dir, paste0("volcano_", safe_name(nm), ".png")), p_vol, width = 8, height = 6, dpi = 300)
    }

    # Primary result (first term)
    res_primary <- results(dds, name = rn[1], alpha = opt$fdr_cutoff)
    res_df <- as.data.frame(res_primary)
    res_df$gene <- rownames(res_df)
    res_df <- annotate_significance(res_df)
    message("[DESeq2] DEG数 (primary): ", sum(res_df$significant, na.rm = TRUE))
    return(list(results = res_df, dds = dds))
  }

  if (opt$contrast_mode == "manual") {
    if (opt$contrast == "") stop("contrast_mode=manual には --contrast が必要です。")
    if (!(opt$contrast %in% rn)) {
      stop("指定した contrast が resultsNames に存在しません: ", opt$contrast)
    }
    res <- results(dds, name = opt$contrast, alpha = opt$fdr_cutoff)
  } else {
    res <- results(dds,
                   contrast = c(condition_col_model, opt$treat_level, opt$ref_level),
                   alpha = opt$fdr_cutoff)
  }

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- annotate_significance(res_df)

  # Interaction terms (additional outputs)
  if (opt$design_mode == "interaction") {
    interaction_terms <- rn[grepl(":", rn)]
    if (length(interaction_terms) > 0) {
      summary_path <- file.path(interp_dir, "interaction_summary.txt")
      cat("Design formula: ", design_formula_str, "\n", file = summary_path)
      cat("Factors: ", paste(factors, collapse = ", "), "\n\n", file = summary_path, append = TRUE)
      cat("Interaction terms:\n", file = summary_path, append = TRUE)
      for (nm in interaction_terms) {
        res_int <- results(dds, name = nm, alpha = opt$fdr_cutoff)
        res_int_df <- as.data.frame(res_int)
        res_int_df$gene <- rownames(res_int_df)
        res_int_df <- annotate_significance(res_int_df)
        res_int_df <- annotate_interaction(res_int_df, nm)
        write.csv(res_int_df,
                  file.path(results_dir, paste0("interaction_results_", safe_name(nm), ".csv")),
                  row.names = FALSE)

        p_int <- ggplot(res_int_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
          geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
          scale_color_manual(values = c("grey60", "red")) +
          geom_hline(yintercept = -log10(opt$fdr_cutoff), linetype = "dashed", color = "blue") +
          geom_vline(xintercept = c(-opt$lfc_cutoff, opt$lfc_cutoff),
                     linetype = "dashed", color = "blue") +
          theme_minimal() +
          labs(title = paste("Interaction Effect Volcano Plot:", nm),
               x = "log2 Fold Change", y = "-log10(p-value)")
        ggsave(file.path(plots_dir, paste0("interaction_volcano_", safe_name(nm), ".pdf")),
               p_int, width = 8, height = 6)
        ggsave(file.path(plots_dir, paste0("interaction_volcano_", safe_name(nm), ".png")),
               p_int, width = 8, height = 6, dpi = 300)

        expl <- interaction_explanation(nm)
        cat("- ", nm, "\n", file = summary_path, append = TRUE)
        cat("  EN: ", expl$en, "\n", file = summary_path, append = TRUE)
        cat("  JP: ", expl$jp, "\n\n", file = summary_path, append = TRUE)
        cat("  Significant genes: ", sum(res_int_df$significant, na.rm = TRUE), "\n\n",
            file = summary_path, append = TRUE)
      }
    }
  }

  message("[DESeq2] DEG数: ", sum(res_df$significant, na.rm = TRUE))
  return(list(results = res_df, dds = dds))
}

# ---- edgeR ----
run_edger <- function() {
  suppressPackageStartupMessages(library(edgeR))
  message("[edgeR] 解析開始...")
  if (opt$design_mode != "simple" || opt$contrast_mode != "pairwise") {
    message("[WARNING] edgeR では現在 simple + pairwise のみをサポートします。")
  }

  group <- meta[[opt$condition_col]]
  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~ group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  res$significant <- (abs(res$logFC) >= opt$lfc_cutoff &
                       res$FDR < opt$fdr_cutoff)

  # 列名を統一
  colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
  colnames(res)[colnames(res) == "PValue"] <- "pvalue"
  colnames(res)[colnames(res) == "FDR"] <- "padj"

  message("[edgeR] DEG数: ", sum(res$significant, na.rm = TRUE))
  return(list(results = res))
}

# ---- limma ----
run_limma <- function() {
  suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
  })
  message("[limma-voom] 解析開始...")
  if (opt$design_mode != "simple" || opt$contrast_mode != "pairwise") {
    message("[WARNING] limma では現在 simple + pairwise のみをサポートします。")
  }

  group <- meta[[opt$condition_col]]
  design <- model.matrix(~ group)

  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf)
  res$gene <- rownames(res)
  res$significant <- (abs(res$logFC) >= opt$lfc_cutoff &
                       res$adj.P.Val < opt$fdr_cutoff)

  colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
  colnames(res)[colnames(res) == "P.Value"] <- "pvalue"
  colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"

  message("[limma] DEG数: ", sum(res$significant, na.rm = TRUE))
  return(list(results = res))
}

# ---- 実行 ----
result <- tryCatch({
  switch(opt$tool,
    "deseq2" = run_deseq2(),
    "edger"  = run_edger(),
    "limma"  = run_limma(),
    stop("Unknown tool: ", opt$tool)
  )
}, error = function(e) {
  message("[ERROR] DE解析に失敗しました: ", e$message)
  stop(e)
})

res_df <- result$results

# ---- 有意遺伝子チェック ----
if (sum(res_df$significant, na.rm = TRUE) == 0) {
  message("[WARNING] 有意なDEGが検出されませんでした")
}

# ---- 結果保存 ----
write.csv(res_df, file.path(results_dir, "de_results.csv"), row.names = FALSE)
write.csv(res_df, file.path(out_sub, "de_results.csv"), row.names = FALSE)
message("[INFO] DE結果保存: ", file.path(results_dir, "de_results.csv"))

# ---- 可視化 ----

# Volcano plot
p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("grey60", "red")) +
  geom_hline(yintercept = -log10(opt$fdr_cutoff), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-opt$lfc_cutoff, opt$lfc_cutoff),
             linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = paste("Volcano Plot (", opt$tool, ")"),
       x = "log2 Fold Change", y = "-log10(p-value)")
ggsave(file.path(out_sub, "volcano_plot.pdf"), p_vol, width = 8, height = 6)
ggsave(file.path(out_sub, "volcano_plot.png"), p_vol, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "volcano_plot.pdf"), p_vol, width = 8, height = 6)
ggsave(file.path(plots_dir, "volcano_plot.png"), p_vol, width = 8, height = 6, dpi = 300)

# MA plot
if ("baseMean" %in% colnames(res_df)) {
  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("grey60", "red")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    theme_minimal() +
    labs(title = paste("MA Plot (", opt$tool, ")"),
         x = "log10(mean expression)", y = "log2 Fold Change")
  ggsave(file.path(out_sub, "ma_plot.pdf"), p_ma, width = 8, height = 6)
  ggsave(file.path(out_sub, "ma_plot.png"), p_ma, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, "ma_plot.pdf"), p_ma, width = 8, height = 6)
  ggsave(file.path(plots_dir, "ma_plot.png"), p_ma, width = 8, height = 6, dpi = 300)
}

# Heatmap (上位DEG)

sig_genes <- res_df[res_df$significant == TRUE, ]
if (nrow(sig_genes) > 0) {
  top_genes <- head(sig_genes[order(sig_genes$padj), "gene"], 50)
  top_genes <- top_genes[!is.na(top_genes)]

  if (length(top_genes) < 2) {
    message("[INFO] ヒートマップをスキップ（遺伝子数不足）")
  } else {
    log_counts <- log2(counts[top_genes, , drop = FALSE] + 1)
    anno_col <- data.frame(row.names = colnames(counts))
    anno_col[[opt$condition_col]] <- meta[[opt$condition_col]]

    pdf(file.path(plots_dir, "heatmap_top_deg.pdf"), width = 10, height = 12)
    pheatmap(log_counts,
             scale = "row",
             annotation_col = anno_col,
             show_rownames = (length(top_genes) <= 50),
             clustering_method = "ward.D2",
             main = paste("Top DEGs Heatmap (", opt$tool, ")"))
    dev.off()

    png(file.path(plots_dir, "heatmap_top_deg.png"), width = 1200, height = 1400, res = 150)
    pheatmap(log_counts,
             scale = "row",
             annotation_col = anno_col,
             show_rownames = (length(top_genes) <= 50),
             clustering_method = "ward.D2",
             main = paste("Top DEGs Heatmap (", opt$tool, ")"))
    dev.off()
  }
}

# サマリー出力
sink(file.path(logs_dir, "summary.txt"))
cat("========================================\n")
cat("差次発現解析サマリー\n")
cat("========================================\n")
cat("ツール:", opt$tool, "\n")
cat("DESeq2 test:", opt$test, "\n")
if (opt$reduced_design != "") cat("Reduced design:", opt$reduced_design, "\n")
cat("条件:", opt$condition_col, "\n")
cat("比較:", opt$treat_level, "vs", opt$ref_level, "\n")
cat("Design:", design_formula_str, "\n")
cat("Factors:", paste(factors, collapse = ", "), "\n")
cat("Contrast mode:", opt$contrast_mode, "\n")
if (opt$contrast != "") cat("Contrast:", opt$contrast, "\n")
if (subset_info$applied) {
  cat("Subset:", opt$subset_col, " in ", paste(subset_info$values, collapse = ","), "\n", sep = "")
}
cat("FDR閾値:", opt$fdr_cutoff, "\n")
cat("log2FC閾値:", opt$lfc_cutoff, "\n")
cat("総遺伝子数:", nrow(res_df), "\n")
cat("有意な DEG 数:", sum(res_df$significant, na.rm = TRUE), "\n")
cat("  Up-regulated:", sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Down-regulated:", sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE), "\n")
sink()

if (explain_results) {
  write_interpretation_report(file.path(interp_dir, "interpretation_report.txt"))
}

# 再現性情報
sink(file.path(logs_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("[完了] 01_differential_expression.R")
