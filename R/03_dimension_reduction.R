#!/usr/bin/env Rscript
# ============================================================
# 03_dimension_reduction.R
# PCA / t-SNE / UMAP — Diagnostic variance analysis
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
})

option_list <- c(base_option_list, list(
  make_option("--methods",         type = "character", default = "pca,umap"),
  make_option("--color_by",        type = "character", default = "condition"),
  make_option("--top_n",           type = "integer",   default = 2000),
  make_option("--perplexity",      type = "integer",   default = 30),
  make_option("--n_neighbors",     type = "integer",   default = 15),
  # New diagnostic options
  make_option("--color_by_all",    type = "character", default = "FALSE"),
  make_option("--split_by",        type = "character", default = ""),
  make_option("--transform",       type = "character", default = "log2"),
  make_option("--feature_method",  type = "character", default = "variance"),
  make_option("--perplexity_list", type = "character", default = ""),
  make_option("--n_neighbors_list",type = "character", default = ""),
  make_option("--label_samples",   type = "character", default = "TRUE")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("03_dimreduc", opt$outdir)

methods <- tolower(trimws(unlist(strsplit(opt$methods, ","))))

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

color_by_all   <- parse_bool(opt$color_by_all)
label_samples  <- parse_bool(opt$label_samples)
split_by       <- if (is.null(opt$split_by) || opt$split_by == "") NULL else opt$split_by
transform_mode <- tolower(opt$transform)
feature_method <- tolower(opt$feature_method)

# ---- データ読み込み ----
counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

# ---- Logging summary ----
message("============================================================")
message("[INFO] 次元削減解析サマリー")
message("[INFO]   サンプル数: ", ncol(counts))
message("[INFO]   総遺伝子数: ", nrow(counts))
message("[INFO]   top_n: ", opt$top_n)
message("[INFO]   color_by: ", opt$color_by)
message("[INFO]   color_by_all: ", color_by_all)
message("[INFO]   split_by: ", ifelse(is.null(split_by), "(none)", split_by))
message("[INFO]   transform: ", transform_mode)
message("[INFO]   feature_method: ", feature_method)
message("[INFO]   label_samples: ", label_samples)
message("[INFO]   methods: ", paste(methods, collapse = ", "))
message("============================================================")

# ---- Transform ----
transform_data <- function(raw_counts, mode) {
  if (mode == "vst" || mode == "rlog") {
    suppressPackageStartupMessages(library(DESeq2))
    counts_int <- round(raw_counts)
    storage.mode(counts_int) <- "integer"
    # minimal DESeq2 object for transform only
    col_data <- data.frame(row.names = colnames(counts_int), dummy = rep(1, ncol(counts_int)))
    dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = col_data, design = ~ 1)
    if (mode == "vst") {
      message("[INFO] VST 変換中...")
      vsd <- vst(dds, blind = TRUE)
      return(assay(vsd))
    } else {
      message("[INFO] rlog 変換中...")
      rld <- rlog(dds, blind = TRUE)
      return(assay(rld))
    }
  }
  # default: log2
  message("[INFO] log2(counts+1) 変換中...")
  return(log2(raw_counts + 1))
}

transformed <- transform_data(counts, transform_mode)

# ---- Feature selection ----
select_features <- function(mat, top_n, method = "variance") {
  if (method == "mad") {
    scores <- apply(mat, 1, mad)
    message("[INFO] MAD ベースで上位 ", min(top_n, length(scores)), " 遺伝子を選択")
  } else {
    scores <- apply(mat, 1, var)
    message("[INFO] 分散ベースで上位 ", min(top_n, length(scores)), " 遺伝子を選択")
  }
  top_genes <- names(sort(scores, decreasing = TRUE))[1:min(top_n, length(scores))]
  return(mat[top_genes, , drop = FALSE])
}

filtered <- select_features(transformed, opt$top_n, feature_method)

color_col <- if (opt$color_by %in% colnames(meta)) opt$color_by else colnames(meta)[1]

# ============================================================
# Helper: plot builder
# ============================================================

make_scatter <- function(df, x_col, y_col, color_col_name, color_vals,
                         title, xlab, ylab, show_labels = TRUE) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                       color = .data[[color_col_name]])) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = xlab, y = ylab, color = color_col_name)
  if (show_labels && "label" %in% colnames(df)) {
    p <- p + geom_text_repel(aes(label = label),
                             size = sample_label_size(nrow(df)),
                             max.overlaps = Inf)
  }
  return(p)
}

save_plot <- function(p, path_no_ext, w = 10, h = 8) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = w, height = h)
  ggsave(paste0(path_no_ext, ".png"), p, width = w, height = h, dpi = 300)
}

# ============================================================
# Core analysis function (can be called per-split or once)
# ============================================================

run_dimreduc <- function(data_mat, meta_sub, out_dir, suffix = "") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Select features for this subset
  sub_filtered <- select_features(data_mat, opt$top_n, feature_method)

  pca_obj <- NULL  # keep for variance analysis

  # ---- PCA ----
  if ("pca" %in% methods) {
    message("[INFO] PCA 実行...", suffix)
    pca_obj <- prcomp(t(sub_filtered), scale. = TRUE)
    var_exp <- summary(pca_obj)$importance[2, ] * 100

    pca_df <- data.frame(
      PC1 = pca_obj$x[, 1], PC2 = pca_obj$x[, 2],
      PC3 = if (ncol(pca_obj$x) >= 3) pca_obj$x[, 3] else 0,
      sample = rownames(meta_sub),
      label = sample_labels[rownames(meta_sub)],
      stringsAsFactors = FALSE
    )
    # attach all metadata columns for multi-factor plotting
    for (cn in colnames(meta_sub)) {
      pca_df[[cn]] <- meta_sub[[cn]]
    }

    # Primary PCA plot
    cc <- if (color_col %in% colnames(pca_df)) color_col else colnames(meta_sub)[1]
    p <- make_scatter(pca_df, "PC1", "PC2", cc, NULL,
                      paste0("PCA", suffix),
                      paste0("PC1 (", round(var_exp[1], 1), "%)"),
                      paste0("PC2 (", round(var_exp[2], 1), "%)"),
                      label_samples)
    save_plot(p, file.path(out_dir, "pca_plot"))

    # PC1 vs PC3
    if (ncol(pca_obj$x) >= 3) {
      p13 <- make_scatter(pca_df, "PC1", "PC3", cc, NULL,
                          paste0("PCA (PC1 vs PC3)", suffix),
                          paste0("PC1 (", round(var_exp[1], 1), "%)"),
                          paste0("PC3 (", round(var_exp[3], 1), "%)"),
                          label_samples)
      save_plot(p13, file.path(out_dir, "pca_pc1_pc3"))
    }

    # Scree plot
    n_pcs <- min(20, length(var_exp))
    scree_df <- data.frame(PC = 1:n_pcs, Variance = var_exp[1:n_pcs])
    p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_line(color = "coral", linewidth = 0.8) +
      geom_point(color = "coral", size = 2) +
      theme_minimal() +
      labs(title = "Scree Plot", x = "Principal Component", y = "% Variance Explained")
    save_plot(p_scree, file.path(out_dir, "scree_plot"), w = 8, h = 5)

    # Loadings & coordinates
    write.csv(pca_obj$rotation[, 1:min(10, ncol(pca_obj$rotation))],
              file.path(out_dir, "pca_loadings.csv"))
    write.csv(pca_df, file.path(out_dir, "pca_coordinates.csv"), row.names = FALSE)

    # ---- Multi-factor PCA ----
    if (color_by_all) {
      message("[INFO] Multi-factor PCA (全メタデータカラム)...")
      multi_dir <- file.path(out_dir, "pca_by_factor")
      dir.create(multi_dir, recursive = TRUE, showWarnings = FALSE)
      for (cn in colnames(meta_sub)) {
        if (cn %in% colnames(pca_df) && length(unique(meta_sub[[cn]])) > 1 &&
            length(unique(meta_sub[[cn]])) <= 50) {
          pca_df[[cn]] <- as.factor(pca_df[[cn]])
          pm <- make_scatter(pca_df, "PC1", "PC2", cn, NULL,
                             paste0("PCA colored by ", cn, suffix),
                             paste0("PC1 (", round(var_exp[1], 1), "%)"),
                             paste0("PC2 (", round(var_exp[2], 1), "%)"),
                             label_samples)
          save_plot(pm, file.path(multi_dir, paste0("pca_by_", make.names(cn))))
        }
      }
    }

    # ---- Variance explained by metadata (CRITICAL) ----
    message("[INFO] PCA分散説明解析（メタデータ因子別）...")
    n_test_pcs <- min(5, ncol(pca_obj$x))
    var_results <- data.frame()

    for (cn in colnames(meta_sub)) {
      vals <- meta_sub[[cn]]
      # Skip columns with too many unique values (likely IDs) or all same
      n_uniq <- length(unique(vals))
      if (n_uniq < 2 || n_uniq > nrow(meta_sub) * 0.8) next

      for (pc_i in 1:n_test_pcs) {
        pc_vals <- pca_obj$x[, pc_i]
        tryCatch({
          # Convert to factor for categorical, keep numeric for continuous
          if (is.numeric(vals) && n_uniq > 10) {
            fit <- lm(pc_vals ~ vals)
          } else {
            fit <- lm(pc_vals ~ as.factor(vals))
          }
          ss <- summary(fit)
          aov_res <- anova(fit)
          r2 <- ss$r.squared
          adj_r2 <- ss$adj.r.squared
          f_stat <- ss$fstatistic[1]
          p_val <- aov_res[1, "Pr(>F)"]

          var_results <- rbind(var_results, data.frame(
            metadata_column = cn,
            PC = paste0("PC", pc_i),
            PC_variance_pct = round(var_exp[pc_i], 2),
            R_squared = round(r2, 4),
            Adj_R_squared = round(adj_r2, 4),
            F_statistic = round(f_stat, 4),
            p_value = signif(p_val, 4),
            significance = ifelse(p_val < 0.001, "***",
                           ifelse(p_val < 0.01, "**",
                           ifelse(p_val < 0.05, "*", "ns"))),
            stringsAsFactors = FALSE
          ))
        }, error = function(e) {
          message("[WARNING] ", cn, " vs PC", pc_i, " スキップ: ", e$message)
        })
      }
    }

    if (nrow(var_results) > 0) {
      var_results <- var_results[order(var_results$PC, -var_results$R_squared), ]
      write.csv(var_results, file.path(out_dir, "pca_variance_explained.csv"),
                row.names = FALSE)
      message("[INFO] 分散説明結果を保存: pca_variance_explained.csv")

      # Print top associations
      message("[INFO] === PCA分散説明サマリー ===")
      for (pc_i in 1:n_test_pcs) {
        pc_name <- paste0("PC", pc_i)
        pc_rows <- var_results[var_results$PC == pc_name, ]
        if (nrow(pc_rows) > 0) {
          top <- pc_rows[1, ]
          message(sprintf("[INFO]   %s (%.1f%%): top factor = %s (R²=%.3f, p=%s %s)",
                          pc_name, top$PC_variance_pct,
                          top$metadata_column, top$R_squared,
                          format(top$p_value, scientific = TRUE),
                          top$significance))
        }
      }

      # Variance explained heatmap
      if (nrow(var_results) > 2) {
        tryCatch({
          r2_mat <- reshape(var_results[, c("metadata_column", "PC", "R_squared")],
                            idvar = "metadata_column", timevar = "PC",
                            direction = "wide")
          rownames(r2_mat) <- r2_mat$metadata_column
          r2_mat$metadata_column <- NULL
          colnames(r2_mat) <- sub("R_squared\\.", "", colnames(r2_mat))

          pdf(file.path(out_dir, "pca_variance_explained_heatmap.pdf"), width = 8, height = 6)
          pheatmap(as.matrix(r2_mat),
                   display_numbers = TRUE, number_format = "%.3f",
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   main = "R² of metadata factors vs PCs",
                   color = colorRampPalette(c("white", "steelblue", "darkblue"))(50))
          dev.off()
          message("[INFO] R²ヒートマップ保存: pca_variance_explained_heatmap.pdf")
        }, error = function(e) {
          message("[WARNING] R²ヒートマップ作成をスキップ: ", e$message)
        })
      }
    }
  }

  # ---- Sample distance heatmap ----
  message("[INFO] サンプル間距離ヒートマップ作成...")
  tryCatch({
    dist_mat <- dist(t(sub_filtered), method = "euclidean")
    dist_matrix <- as.matrix(dist_mat)

    # Annotation
    anno_df <- data.frame(row.names = rownames(meta_sub))
    for (cn in colnames(meta_sub)) {
      n_uniq <- length(unique(meta_sub[[cn]]))
      if (n_uniq >= 2 && n_uniq <= 20) {
        anno_df[[cn]] <- as.factor(meta_sub[[cn]])
      }
    }
    if (ncol(anno_df) == 0) anno_df <- NA

    pdf(file.path(out_dir, "sample_distance_heatmap.pdf"), width = 12, height = 10)
    pheatmap(dist_matrix,
             clustering_distance_rows = dist_mat,
             clustering_distance_cols = dist_mat,
             clustering_method = "ward.D2",
             annotation_row = anno_df,
             annotation_col = anno_df,
             labels_row = sample_labels[rownames(meta_sub)],
             labels_col = sample_labels[rownames(meta_sub)],
             fontsize_row = sample_label_fontsize(nrow(meta_sub)),
             fontsize_col = sample_label_fontsize(nrow(meta_sub)),
             main = paste0("Sample Distance Heatmap (", transform_mode, ")"),
             show_rownames = (ncol(data_mat) <= 60),
             show_colnames = (ncol(data_mat) <= 60))
    dev.off()

    png(file.path(out_dir, "sample_distance_heatmap.png"), width = 1400, height = 1200, res = 150)
    pheatmap(dist_matrix,
             clustering_distance_rows = dist_mat,
             clustering_distance_cols = dist_mat,
             clustering_method = "ward.D2",
             annotation_row = anno_df,
             annotation_col = anno_df,
             labels_row = sample_labels[rownames(meta_sub)],
             labels_col = sample_labels[rownames(meta_sub)],
             fontsize_row = sample_label_fontsize(nrow(meta_sub)),
             fontsize_col = sample_label_fontsize(nrow(meta_sub)),
             main = paste0("Sample Distance Heatmap (", transform_mode, ")"),
             show_rownames = (ncol(data_mat) <= 60),
             show_colnames = (ncol(data_mat) <= 60))
    dev.off()
    message("[INFO] サンプル距離ヒートマップ保存完了")
  }, error = function(e) {
    message("[WARNING] サンプル距離ヒートマップ作成をスキップ: ", e$message)
  })

  # ---- t-SNE ----
  if ("tsne" %in% methods) {
    suppressPackageStartupMessages(library(Rtsne))

    perp_values <- opt$perplexity
    if (!is.null(opt$perplexity_list) && opt$perplexity_list != "") {
      perp_values <- as.integer(trimws(unlist(strsplit(opt$perplexity_list, ","))))
    }

    for (perp_raw in perp_values) {
      perp <- min(perp_raw, floor((ncol(sub_filtered) - 1) / 3))
      if (perp < 1) {
        message("[WARNING] perplexity=", perp_raw, " スキップ (サンプル数不足)")
        next
      }
      message("[INFO] t-SNE 実行 (perplexity=", perp, ")...", suffix)
      set.seed(42)
      tsne <- Rtsne(t(sub_filtered), dims = 2, perplexity = perp, verbose = FALSE)

      tsne_df <- data.frame(
        tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2],
        sample = rownames(meta_sub),
        label = sample_labels[rownames(meta_sub)],
        stringsAsFactors = FALSE
      )
      for (cn in colnames(meta_sub)) tsne_df[[cn]] <- meta_sub[[cn]]

      cc <- if (color_col %in% colnames(tsne_df)) color_col else colnames(meta_sub)[1]
      p_tsne <- make_scatter(tsne_df, "tSNE1", "tSNE2", cc, NULL,
                             paste0("t-SNE (perplexity=", perp, ")", suffix),
                             "tSNE1", "tSNE2", label_samples)

      tag <- if (length(perp_values) > 1) paste0("_perp_", perp) else ""
      save_plot(p_tsne, file.path(out_dir, paste0("tsne_plot", tag)))
      write.csv(tsne_df, file.path(out_dir, paste0("tsne_coordinates", tag, ".csv")),
                row.names = FALSE)
    }
  }

  # ---- UMAP ----
  if ("umap" %in% methods) {
    suppressPackageStartupMessages(library(uwot))

    nn_values <- opt$n_neighbors
    if (!is.null(opt$n_neighbors_list) && opt$n_neighbors_list != "") {
      nn_values <- as.integer(trimws(unlist(strsplit(opt$n_neighbors_list, ","))))
    }

    for (nn in nn_values) {
      if (nn >= ncol(sub_filtered)) {
        message("[WARNING] n_neighbors=", nn, " スキップ (サンプル数不足)")
        next
      }
      message("[INFO] UMAP 実行 (n_neighbors=", nn, ")...", suffix)
      set.seed(42)
      umap_res <- uwot::umap(t(sub_filtered), n_neighbors = nn, min_dist = 0.3)

      umap_df <- data.frame(
        UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2],
        sample = rownames(meta_sub),
        label = sample_labels[rownames(meta_sub)],
        stringsAsFactors = FALSE
      )
      for (cn in colnames(meta_sub)) umap_df[[cn]] <- meta_sub[[cn]]

      cc <- if (color_col %in% colnames(umap_df)) color_col else colnames(meta_sub)[1]
      p_umap <- make_scatter(umap_df, "UMAP1", "UMAP2", cc, NULL,
                             paste0("UMAP (n_neighbors=", nn, ")", suffix),
                             "UMAP1", "UMAP2", label_samples)

      tag <- if (length(nn_values) > 1) paste0("_nn_", nn) else ""
      save_plot(p_umap, file.path(out_dir, paste0("umap_plot", tag)))
      write.csv(umap_df, file.path(out_dir, paste0("umap_coordinates", tag, ".csv")),
                row.names = FALSE)
    }
  }
}

# ============================================================
# Run analysis: full dataset
# ============================================================

message("[INFO] 全データで次元削減を実行...")
run_dimreduc(transformed, meta, out_sub, suffix = "")

# ============================================================
# Run analysis: split by factor
# ============================================================

if (!is.null(split_by)) {
  if (!(split_by %in% colnames(meta))) {
    message("[WARNING] split_by='", split_by, "' はメタデータに存在しません。スキップします。")
  } else {
    levels_split <- unique(meta[[split_by]])
    message("[INFO] split_by='", split_by, "' で分割解析 (", length(levels_split), " groups)...")

    for (lvl in levels_split) {
      lvl_safe <- gsub("[^A-Za-z0-9_\\-]", "_", as.character(lvl))
      sub_dir <- file.path(out_sub, paste0(split_by, "_", lvl_safe))
      idx <- meta[[split_by]] == lvl
      n_sub <- sum(idx)

      if (n_sub < 3) {
        message("[WARNING] ", split_by, "=", lvl, " サンプル数不足 (n=", n_sub, ")。スキップ。")
        next
      }

      message("[INFO] --- ", split_by, " = ", lvl, " (n=", n_sub, ") ---")
      sub_data <- transformed[, idx, drop = FALSE]
      sub_meta <- meta[idx, , drop = FALSE]
      run_dimreduc(sub_data, sub_meta, sub_dir, suffix = paste0(" [", split_by, "=", lvl, "]"))
    }
  }
}

message("[完了] 03_dimension_reduction.R")
