#!/usr/bin/env Rscript
# ============================================================
# 02_sample_clustering.R
# Sample & Gene Clustering / サンプル・遺伝子クラスタリング解析
# ============================================================
#
# PURPOSE:
#   Mode-selective clustering analysis.
#   Modes:
#     sample       — Clusters samples (hierarchical, k-means, PCA).
#                    Preserves original Step 2 behavior.
#     gene_module  — Clusters genes into modules via k-means,
#                    visualized with pheatmap.
#     trajectory   — Calculates per-gene Z-scores and plots
#                    faceted trend lines grouped by a metadata column.
#
# INPUTS:  --counts, --metadata (standard)
# PARAMS:  --analysis_mode (sample|gene_module|trajectory),
#          --method, --k, --dist_method, --top_n,
#          --target_genes, --trend_x_col, --group_col
# OUTPUTS:
#   sample:      dendrogram, distance heatmap, k-means PCA, elbow plot
#   gene_module: gene_module_heatmap.png/pdf, gene_module_assignments.csv
#   trajectory:  cluster_trends.png/pdf
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
  library(pheatmap)
  library(ggdendro)
  library(ggrepel)
  library(stats)
})

option_list <- c(base_option_list, list(
  make_option("--analysis_mode", type = "character", default = "sample",
              help = "Analysis mode: sample | gene_module | trajectory"),
  make_option("--method",      type = "character", default = "both"),
  make_option("--k",           type = "integer",   default = 3),
  make_option("--dist_method", type = "character", default = "euclidean"),
  make_option("--top_n",       type = "integer",   default = 1000),
  make_option("--target_genes", type = "character", default = ""),
  make_option("--trend_x_col", type = "character", default = ""),
  make_option("--group_col",   type = "character", default = "")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("02_clustering", opt$outdir)

# ---- データ読み込み ----
counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

# ---- Target genes helper ----
get_target_genes <- function(path) {
  if (path == "" || is.na(path)) {
    return(list(genes = character(0), symbol_map = NULL))
  }
  if (!file.exists(path)) {
    stop("[ERROR] target_genes file not found: ", path)
  }

  target_df <- read.csv(path, stringsAsFactors = FALSE)
  if (!("gene" %in% colnames(target_df))) {
    stop("[ERROR] target_genes CSV には 'gene' 列が必要です。")
  }

  keep <- rep(TRUE, nrow(target_df))
  if ("significant" %in% colnames(target_df)) {
    sig_vals <- tolower(trimws(as.character(target_df$significant)))
    keep <- sig_vals %in% c("true", "1", "yes", "y")
  } else if ("padj" %in% colnames(target_df)) {
    keep <- !is.na(target_df$padj) & target_df$padj < 0.05
  } else if ("pvalue" %in% colnames(target_df)) {
    keep <- !is.na(target_df$pvalue) & target_df$pvalue < 0.05
  } else {
    message("[WARNING] significant/padj/pvalue 列がないため gene 列の全遺伝子を使用します。")
  }

  filtered_df <- target_df[keep, , drop = FALSE]
  genes <- unique(filtered_df$gene)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) {
    stop("[ERROR] target_genes から使用可能な遺伝子を抽出できませんでした。")
  }

  # Build symbol map if available
  symbol_map <- NULL
  if ("symbol" %in% colnames(filtered_df)) {
    map_rows <- filtered_df[!is.na(filtered_df$symbol) & filtered_df$symbol != "", ]
    if (nrow(map_rows) > 0) {
      symbol_map <- setNames(map_rows$symbol, map_rows$gene)
      # Deduplicate: keep first occurrence
      symbol_map <- symbol_map[!duplicated(names(symbol_map))]
    }
  }

  list(genes = genes, symbol_map = symbol_map)
}

# 正規化 & 上位変動遺伝子選択
log_counts <- log2(counts + 1)
target_result <- get_target_genes(opt$target_genes)
target_genes <- target_result$genes
symbol_map <- target_result$symbol_map

if (length(target_genes) > 0) {
  matched_genes <- intersect(target_genes, rownames(log_counts))
  if (length(matched_genes) == 0) {
    stop("[ERROR] target_genes と count matrix の遺伝子名が一致しません。")
  }
  filtered <- log_counts[matched_genes, , drop = FALSE]
  message("[INFO] target_genes 由来の ", length(matched_genes),
          " genes をクラスタリングに使用します。")

  # Apply symbol mapping to row names if available
  if (!is.null(symbol_map)) {
    current_names <- rownames(filtered)
    new_names <- ifelse(current_names %in% names(symbol_map),
                        symbol_map[current_names],
                        current_names)
    # Avoid duplicate row names
    if (!any(duplicated(new_names))) {
      rownames(filtered) <- new_names
      message("[INFO] 遺伝子名を symbol に変換しました (", sum(current_names != new_names), " genes)")
    } else {
      message("[WARNING] symbol 変換で重複が生じるため、元の遺伝子名を維持します。")
    }
  }
} else {
  filtered <- select_top_variable(log_counts, opt$top_n)
}

# ── Validate analysis_mode ──
analysis_mode <- tolower(trimws(opt$analysis_mode))
if (!(analysis_mode %in% c("sample", "gene_module", "trajectory"))) {
  stop("[ERROR] --analysis_mode は sample / gene_module / trajectory のいずれかです。")
}
message("[INFO] Analysis mode: ", analysis_mode)


# ==============================================================
# MODE: sample  (preserves existing behavior)
# ==============================================================
if (analysis_mode == "sample") {

  # ---- 距離行列 ----
  if (opt$dist_method == "pearson") {
    d <- as.dist(1 - cor(filtered, method = "pearson"))
  } else if (opt$dist_method == "spearman") {
    d <- as.dist(1 - cor(filtered, method = "spearman"))
  } else {
    d <- dist(t(filtered), method = "euclidean")
  }

  # ---- 階層クラスタリング ----
  if (opt$method %in% c("hierarchical", "both")) {
    message("[INFO] 階層クラスタリング実行...")
    hc <- hclust(d, method = "ward.D2")
    dendro_cex <- sample_label_cex(length(sample_labels))

    # デンドログラム
    pdf(file.path(out_sub, "dendrogram.pdf"), width = 12, height = 8)
    plot(hc, main = "Sample Hierarchical Clustering", xlab = "", sub = "",
         labels = sample_labels[hc$labels], cex = dendro_cex)
    dev.off()

    png(file.path(out_sub, "dendrogram.png"), width = 1400, height = 800, res = 150)
    plot(hc, main = "Sample Hierarchical Clustering", xlab = "", sub = "",
         labels = sample_labels[hc$labels], cex = dendro_cex)
    dev.off()

    # 距離行列ヒートマップ
    anno_col <- data.frame(row.names = colnames(filtered))
    for (cn in colnames(meta)) {
      if (is.character(meta[[cn]]) || is.factor(meta[[cn]])) {
        anno_col[[cn]] <- meta[[cn]]
      }
    }

    pdf(file.path(out_sub, "distance_heatmap.pdf"), width = 10, height = 10)
    pheatmap(as.matrix(d),
             clustering_method = "ward.D2",
             annotation_col = anno_col,
             annotation_row = anno_col,
             labels_col = sample_labels[colnames(filtered)],
             labels_row = sample_labels[colnames(filtered)],
             fontsize_col = sample_label_fontsize(length(sample_labels)),
             fontsize_row = sample_label_fontsize(length(sample_labels)),
             main = "Sample Distance Matrix")
    dev.off()
  }

  # ---- k-means クラスタリング ----
  if (opt$method %in% c("kmeans", "both")) {
    unique_profiles <- nrow(unique(as.data.frame(t(filtered))))
    k_centers <- max(1L, min(as.integer(opt$k), nrow(meta), unique_profiles))
    if (!identical(k_centers, as.integer(opt$k))) {
      message("[WARNING] 指定した k=", opt$k, " を k=", k_centers,
              " に調整しました（サンプル数/一意プロファイル数の制約）。")
    }
    message("[INFO] k-means クラスタリング実行 (k=", k_centers, ")...")
    set.seed(42)
    km <- kmeans(t(filtered), centers = k_centers, nstart = 25)

    meta$kmeans_cluster <- as.factor(km$cluster)

    # PCA で可視化
    pca <- prcomp(t(filtered), scale. = TRUE)
    pca_df <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      cluster = meta$kmeans_cluster,
      sample = rownames(meta),
      label = sample_labels[rownames(meta)]
    )
    var_exp <- summary(pca)$importance[2, 1:2] * 100

    p_km <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = label),
                      size = sample_label_size(length(sample_labels)),
                      max.overlaps = Inf) +
      theme_minimal() +
      labs(title = paste("k-means Clustering (k =", k_centers, ")"),
           x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
           y = paste0("PC2 (", round(var_exp[2], 1), "%)"))
    ggsave(file.path(out_sub, "kmeans_pca.pdf"), p_km, width = 10, height = 8)
    ggsave(file.path(out_sub, "kmeans_pca.png"), p_km, width = 10, height = 8, dpi = 300)

    # Elbow plot
    max_k <- min(10, nrow(meta), unique_profiles)
    wss <- sapply(seq_len(max_k), function(k) {
      tryCatch(
        kmeans(t(filtered), centers = k, nstart = 10)$tot.withinss,
        error = function(e) NA_real_
      )
    })
    elbow_df <- data.frame(k = seq_len(max_k), wss = wss)
    elbow_df <- elbow_df[!is.na(elbow_df$wss), , drop = FALSE]
    p_elbow <- ggplot(elbow_df, aes(x = k, y = wss)) +
      geom_line() + geom_point() +
      theme_minimal() +
      labs(title = "Elbow Plot", x = "Number of Clusters (k)",
           y = "Total Within Sum of Squares")
    ggsave(file.path(out_sub, "elbow_plot.pdf"), p_elbow, width = 8, height = 6)

    # 結果保存
    write.csv(data.frame(sample = rownames(meta), cluster = km$cluster),
              file.path(out_sub, "kmeans_assignments.csv"), row.names = FALSE)

    if (!is.null(opt$trend_x_col) && trimws(opt$trend_x_col) != "") {
      if (!(opt$trend_x_col %in% colnames(meta))) {
        stop("[ERROR] trend_x_col が metadata に存在しません: ", opt$trend_x_col)
      }

      trend_x <- meta[[opt$trend_x_col]]
      if (all(is.na(trend_x))) {
        stop("[ERROR] trend_x_col='", opt$trend_x_col, "' がすべて NA です。")
      }

      trend_df <- data.frame(
        sample = rownames(meta),
        label = sample_labels[rownames(meta)],
        cluster = meta$kmeans_cluster,
        x_value = trend_x,
        mean_expression = colMeans(filtered[, rownames(meta), drop = FALSE], na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      if (is.numeric(trend_df$x_value)) {
        trend_df <- trend_df[order(trend_df$x_value, trend_df$cluster), , drop = FALSE]
      } else {
        lvl_order <- unique(as.character(trend_df$x_value))
        trend_df$x_value <- factor(trend_df$x_value, levels = lvl_order)
        trend_df <- trend_df[order(trend_df$x_value, trend_df$cluster), , drop = FALSE]
      }

      p_trend <- ggplot(trend_df, aes(x = x_value, y = mean_expression,
                                      color = cluster, group = cluster)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1) +
        stat_summary(fun = mean, geom = "point", size = 3) +
        geom_point(alpha = 0.5, position = position_jitter(width = 0.08, height = 0)) +
        geom_text_repel(aes(label = label),
                        size = sample_label_size(length(sample_labels)) - 0.3,
                        max.overlaps = 25,
                        show.legend = FALSE) +
        theme_minimal() +
        labs(title = paste("Expression Trend by Cluster:", opt$trend_x_col),
             x = opt$trend_x_col,
             y = "Mean log2 expression across clustering genes")
      ggsave(file.path(out_sub, "cluster_expression_trend.pdf"), p_trend, width = 11, height = 8)
      ggsave(file.path(out_sub, "cluster_expression_trend.png"), p_trend, width = 11, height = 8, dpi = 300)
    } else {
      message("[INFO] trend_x_col が未指定のため trend plot はスキップします。")
    }
  }


# ==============================================================
# MODE: gene_module  (NEW — gene-wise k-means + pheatmap)
# ==============================================================
} else if (analysis_mode == "gene_module") {

  message("[INFO] Gene module clustering (gene-wise k-means)")

  # Gene-wise k-means: rows = genes, cols = samples
  gene_mat <- as.matrix(filtered)  # genes x samples
  n_genes <- nrow(gene_mat)
  k_centers <- max(1L, min(as.integer(opt$k), n_genes))
  if (!identical(k_centers, as.integer(opt$k))) {
    message("[WARNING] k=", opt$k, " を k=", k_centers, " に調整しました（遺伝子数の制約）。")
  }

  message("[INFO] Gene k-means: ", n_genes, " genes, k=", k_centers)
  set.seed(42)
  km_genes <- kmeans(gene_mat, centers = k_centers, nstart = 25)

  # Assign modules
  gene_modules <- data.frame(
    gene = rownames(gene_mat),
    module = km_genes$cluster,
    stringsAsFactors = FALSE
  )
  gene_modules <- gene_modules[order(gene_modules$module, gene_modules$gene), ]
  write.csv(gene_modules, file.path(out_sub, "gene_module_assignments.csv"), row.names = FALSE)
  message("[INFO] Gene module assignments saved: gene_module_assignments.csv")

  # Reorder matrix by module for heatmap
  ordered_genes <- gene_modules$gene
  heatmap_mat <- gene_mat[ordered_genes, , drop = FALSE]

  # Scale rows (Z-score) for visualization
  heatmap_scaled <- t(scale(t(heatmap_mat)))
  # Clamp extreme values for visualization
  heatmap_scaled[heatmap_scaled > 3] <- 3
  heatmap_scaled[heatmap_scaled < -3] <- -3
  heatmap_scaled[is.na(heatmap_scaled)] <- 0

  # Annotation data frames
  anno_row <- data.frame(
    Module = as.factor(gene_modules$module),
    row.names = gene_modules$gene
  )
  anno_col <- data.frame(row.names = colnames(heatmap_mat))
  for (cn in colnames(meta)) {
    if (is.character(meta[[cn]]) || is.factor(meta[[cn]])) {
      anno_col[[cn]] <- meta[[cn]]
    }
  }

  # Decide fontsize for gene labels
  gene_fontsize <- if (n_genes <= 50) 6 else if (n_genes <= 100) 4 else 2
  show_gene_names <- n_genes <= 200

  # PDF version
  pdf(file.path(out_sub, "gene_module_heatmap.pdf"),
      width = 10, height = max(8, min(20, n_genes * 0.06)))
  pheatmap(heatmap_scaled,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           clustering_method = "ward.D2",
           annotation_row = anno_row,
           annotation_col = anno_col,
           labels_col = sample_labels[colnames(heatmap_mat)],
           fontsize_col = sample_label_fontsize(length(sample_labels)),
           fontsize_row = gene_fontsize,
           show_rownames = show_gene_names,
           main = paste("Gene Modules (k =", k_centers, ")"),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()

  # PNG version
  png(file.path(out_sub, "gene_module_heatmap.png"),
      width = 1200, height = max(800, min(2400, n_genes * 8)), res = 150)
  pheatmap(heatmap_scaled,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           clustering_method = "ward.D2",
           annotation_row = anno_row,
           annotation_col = anno_col,
           labels_col = sample_labels[colnames(heatmap_mat)],
           fontsize_col = sample_label_fontsize(length(sample_labels)),
           fontsize_row = gene_fontsize,
           show_rownames = show_gene_names,
           main = paste("Gene Modules (k =", k_centers, ")"),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()

  message("[INFO] Gene module heatmap saved.")


# ==============================================================
# MODE: trajectory  (NEW — Z-score faceted trend plot)
# ==============================================================
} else if (analysis_mode == "trajectory") {

  # Validate required arguments
  if (is.null(opt$trend_x_col) || trimws(opt$trend_x_col) == "") {
    stop("[ERROR] trajectory モードでは --trend_x_col が必須です。")
  }
  if (!(opt$trend_x_col %in% colnames(meta))) {
    stop("[ERROR] trend_x_col が metadata に存在しません: ", opt$trend_x_col)
  }
  if (is.null(opt$group_col) || trimws(opt$group_col) == "") {
    stop("[ERROR] trajectory モードでは --group_col が必須です。")
  }
  if (!(opt$group_col %in% colnames(meta))) {
    stop("[ERROR] group_col が metadata に存在しません: ", opt$group_col)
  }

  message("[INFO] Trajectory trend plot (Z-score faceted)")

  # Gene-wise k-means to assign clusters (reuse k)
  gene_mat <- as.matrix(filtered)  # genes x samples
  n_genes <- nrow(gene_mat)
  k_centers <- max(1L, min(as.integer(opt$k), n_genes))
  message("[INFO] Gene k-means for trajectory: ", n_genes, " genes, k=", k_centers)
  set.seed(42)
  km_genes <- kmeans(gene_mat, centers = k_centers, nstart = 25)

  gene_clusters <- data.frame(
    gene = rownames(gene_mat),
    cluster = km_genes$cluster,
    stringsAsFactors = FALSE
  )
  write.csv(gene_clusters, file.path(out_sub, "trajectory_gene_clusters.csv"), row.names = FALSE)

  # Per-gene Z-scores
  gene_means <- rowMeans(gene_mat)
  gene_sds   <- apply(gene_mat, 1, sd)
  gene_sds[gene_sds == 0] <- 1  # avoid division by zero
  z_mat <- (gene_mat - gene_means) / gene_sds

  # Reshape to long format using base R
  n_samples <- ncol(z_mat)
  gene_names <- rownames(z_mat)
  sample_names <- colnames(z_mat)

  long_df <- data.frame(
    gene = rep(gene_names, times = n_samples),
    sample = rep(sample_names, each = n_genes),
    zscore = as.vector(z_mat),
    stringsAsFactors = FALSE
  )

  # Add cluster assignment
  long_df$cluster <- gene_clusters$cluster[match(long_df$gene, gene_clusters$gene)]
  long_df$cluster <- paste0("Module ", long_df$cluster)

  # Add metadata columns
  long_df$x_value <- meta[long_df$sample, opt$trend_x_col]
  long_df$group   <- meta[long_df$sample, opt$group_col]

  # Handle x_value ordering
  if (is.numeric(long_df$x_value)) {
    long_df <- long_df[order(long_df$x_value), ]
  } else {
    lvl_order <- unique(as.character(meta[[opt$trend_x_col]]))
    long_df$x_value <- factor(long_df$x_value, levels = lvl_order)
  }
  long_df$group <- as.factor(long_df$group)

  # Faceted trend plot
  p_trend <- ggplot(long_df, aes(x = x_value, y = zscore, color = group)) +
    # Individual gene trajectories (low alpha)
    geom_line(aes(group = interaction(gene, group)), alpha = 0.08, linewidth = 0.3) +
    # Mean group trajectory (bold)
    stat_summary(aes(group = group), fun = mean, geom = "line",
                 linewidth = 1.2) +
    stat_summary(aes(group = group), fun = mean, geom = "point",
                 size = 2.5) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = "Gene Expression Trajectories by Module",
      x = opt$trend_x_col,
      y = "Z-score (per-gene normalized)",
      color = opt$group_col
    )

  # Adjust width based on number of facets
  plot_width <- max(10, k_centers * 4)
  plot_height <- max(6, ceiling(k_centers / 3) * 4)

  ggsave(file.path(out_sub, "cluster_trends.pdf"), p_trend,
         width = plot_width, height = plot_height)
  ggsave(file.path(out_sub, "cluster_trends.png"), p_trend,
         width = plot_width, height = plot_height, dpi = 300)

  message("[INFO] Trajectory trend plot saved: cluster_trends.png / cluster_trends.pdf")
}

message("[完了] 02_sample_clustering.R (mode: ", analysis_mode, ")")
