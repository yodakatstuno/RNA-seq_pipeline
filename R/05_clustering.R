#!/usr/bin/env Rscript
# ============================================================
# 05_clustering.R
# Time-series gene clustering / 時系列遺伝子クラスタリング
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
  make_option("--normalization_method", type = "character", default = "vst"),
  make_option("--cluster_method", type = "character", default = "both"),
  make_option("--k", type = "integer", default = 6),
  make_option("--dist_method", type = "character", default = "euclidean"),
  make_option("--top_n", type = "integer", default = 1000),
  make_option("--top_method", type = "character", default = "variance"),
  make_option("--target_genes", type = "character", default = ""),
  make_option("--time_col", type = "character", default = ""),
  make_option("--group_col", type = "character", default = ""),
  make_option("--force_rlog", type = "character", default = "FALSE"),
  make_option("--rlog_max_samples", type = "integer", default = 30),
  make_option("--zscore_clamp", type = "double", default = 3)
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("05_clustering", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata
sample_labels <- build_sample_labels(meta, opt$label_mode, opt$symbol_col)

normalization_method <- tolower(trimws(opt$normalization_method))
cluster_method <- tolower(trimws(opt$cluster_method))
force_rlog <- parse_bool(opt$force_rlog)

allowed_norm <- c("deseq2", "tmm", "log2", "vst", "rlog", "cpm")
if (!(normalization_method %in% allowed_norm)) {
  stop("[ERROR] normalization_method は ", paste(allowed_norm, collapse = ", "),
       " のいずれかです。")
}

allowed_cluster <- c("hierarchical", "kmeans", "both")
if (!(cluster_method %in% allowed_cluster)) {
  stop("[ERROR] cluster_method は hierarchical / kmeans / both のいずれかです。")
}

top_method <- tolower(trimws(opt$top_method))
allowed_top_methods <- c("variance", "mean", "deg")
if (!(top_method %in% allowed_top_methods)) {
  stop("[ERROR] top_method は variance / mean / deg のいずれかです。")
}

get_target_genes <- function(path) {
  if (path == "" || is.na(path)) {
    return(character(0))
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
  }

  genes <- unique(target_df$gene[keep])
  genes[!is.na(genes) & genes != ""]
}

normalized_stem <- switch(
  normalization_method,
  "deseq2" = "normalized_deseq2",
  "tmm" = "normalized_tmm",
  "log2" = "normalized_log2",
  "vst" = "normalized_vst",
  "rlog" = "normalized_rlog",
  "cpm" = "normalized_cpm"
)

existing_normalized_path <- file.path(opt$outdir, "04_normalization", paste0(normalized_stem, ".rds"))
if (file.exists(existing_normalized_path)) {
  message("[INFO] 既存の正規化データを使用: ", existing_normalized_path)
  normalized_mat <- readRDS(existing_normalized_path)
} else {
  message("[INFO] 正規化データを新規作成: ", normalization_method)
  normalized_mat <- normalize_counts_matrix(
    counts = counts,
    meta = meta,
    method = normalization_method,
    force_rlog = force_rlog,
    rlog_max_samples = opt$rlog_max_samples
  )$matrix
}

target_genes <- get_target_genes(opt$target_genes)
select_genes_for_clustering <- function(mat, top_n, top_method, target_genes = character(0)) {
  message("[INFO] Gene selection method: ", top_method, " (top_n=", top_n, ")")

  if (top_method == "deg") {
    if (length(target_genes) == 0) {
      message("[WARN] top_method=deg が指定されましたが target_genes が未指定のため、variance selection にフォールバックします。")
      top_method <- "variance"
    } else {
      matched_genes <- intersect(target_genes, rownames(mat))
      if (length(matched_genes) == 0) {
        stop("[ERROR] target_genes と正規化行列の遺伝子名が一致しません。")
      }
      message("[INFO] DEG-based selection: ", length(matched_genes), " genes")
      return(mat[matched_genes, , drop = FALSE])
    }
  }

  if (top_method == "mean") {
    gene_scores <- rowMeans(mat, na.rm = TRUE)
    top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:min(top_n, length(gene_scores))]
    message("[INFO] Mean-expression based selection: ", length(top_genes), " genes")
    return(mat[top_genes, , drop = FALSE])
  }

  gene_scores <- apply(mat, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:min(top_n, length(gene_scores))]
  message("[INFO] Variance-based selection: ", length(top_genes), " genes")
  mat[top_genes, , drop = FALSE]
}

filtered_mat <- select_genes_for_clustering(
  mat = normalized_mat,
  top_n = opt$top_n,
  top_method = top_method,
  target_genes = target_genes
)

message("[INFO] Z-score is treated as a downstream clustering transform and is recomputed internally.")
z_mat <- row_zscore_matrix(filtered_mat)

detect_time_col <- function(meta, requested) {
  if (!is.null(requested) && trimws(requested) != "") {
    if (!(requested %in% colnames(meta))) {
      stop("[ERROR] time_col が metadata に存在しません: ", requested)
    }
    return(requested)
  }
  autodetect_metadata_column(
    meta,
    preferred = c("dpf", "time", "timepoint", "stage", "day", "days"),
    patterns = c("(^|_)dpf$", "time", "stage", "day")
  )
}

detect_group_col <- function(meta, requested) {
  if (!is.null(requested) && trimws(requested) != "") {
    if (!(requested %in% colnames(meta))) {
      stop("[ERROR] group_col が metadata に存在しません: ", requested)
    }
    return(requested)
  }
  autodetect_metadata_column(
    meta,
    preferred = c("genotype", "condition", "group", "treatment"),
    patterns = c("genotype", "condition", "group", "treat")
  )
}

time_col <- detect_time_col(meta, opt$time_col)
group_col <- detect_group_col(meta, opt$group_col)
overlay_col <- autodetect_metadata_column(
  meta,
  preferred = c("condition"),
  patterns = c("condition", "treat", "group")
)
if (!is.null(overlay_col) && identical(overlay_col, group_col)) {
  overlay_col <- NULL
}
message("[INFO] time_col: ", ifelse(is.null(time_col), "(none)", time_col))
message("[INFO] group_col: ", ifelse(is.null(group_col), "(none)", group_col))
message("[INFO] overlay_col: ", ifelse(is.null(overlay_col), "(none)", overlay_col))
message("[INFO] Clustering parameters: normalization=", normalization_method,
        ", method=", cluster_method,
        ", k=", opt$k,
        ", dist=", opt$dist_method,
        ", top_method=", top_method)

order_samples <- function(meta, time_col = NULL, group_col = NULL) {
  order_df <- data.frame(sample = rownames(meta), stringsAsFactors = FALSE)

  if (!is.null(group_col)) {
    order_df$group_val <- as.character(meta[[group_col]])
  } else {
    order_df$group_val <- order_df$sample
  }

  if (!is.null(time_col)) {
    time_vals <- meta[[time_col]]
    if (is.numeric(time_vals)) {
      order_df$time_val <- time_vals
    } else {
      order_df$time_val <- factor(as.character(time_vals), levels = unique(as.character(time_vals)))
    }
  } else {
    order_df$time_val <- seq_len(nrow(order_df))
  }

  order_df <- order_df[order(order_df$group_val, order_df$time_val, order_df$sample), , drop = FALSE]
  order_df$sample
}

sample_order <- order_samples(meta, time_col = time_col, group_col = group_col)
z_mat <- z_mat[, sample_order, drop = FALSE]
meta <- meta[sample_order, , drop = FALSE]

make_distance <- function(mat, dist_method) {
  if (dist_method == "pearson") {
    return(as.dist(1 - cor(t(mat), method = "pearson", use = "pairwise.complete.obs")))
  }
  if (dist_method == "spearman") {
    return(as.dist(1 - cor(t(mat), method = "spearman", use = "pairwise.complete.obs")))
  }
  dist(mat, method = "euclidean")
}

safe_k <- function(k, mat) {
  unique_profiles <- nrow(unique(as.data.frame(mat)))
  max(1L, min(as.integer(k), nrow(mat), unique_profiles))
}

build_annotation_col <- function(meta) {
  anno_col <- data.frame(row.names = rownames(meta))
  for (cn in colnames(meta)) {
    if (is.character(meta[[cn]]) || is.factor(meta[[cn]])) {
      anno_col[[cn]] <- meta[[cn]]
    }
  }
  anno_col
}

save_cluster_outputs <- function(assignments, method_name) {
  ordered_genes <- assignments$gene
  clustered_mat <- z_mat[ordered_genes, , drop = FALSE]
  saveRDS(clustered_mat, file.path(out_sub, paste0("clustered_matrix_", method_name, ".rds")))
  write.csv(clustered_mat, file.path(out_sub, paste0("clustered_matrix_", method_name, ".csv")))
  write.csv(assignments, file.path(out_sub, paste0("cluster_assignments_", method_name, ".csv")),
            row.names = FALSE)

  heatmap_mat <- clustered_mat
  heatmap_mat[heatmap_mat > opt$zscore_clamp] <- opt$zscore_clamp
  heatmap_mat[heatmap_mat < -opt$zscore_clamp] <- -opt$zscore_clamp
  heatmap_mat[is.na(heatmap_mat)] <- 0

  anno_row <- data.frame(
    Cluster = factor(assignments$cluster),
    row.names = assignments$gene
  )
  anno_col <- build_annotation_col(meta)
  sample_font <- sample_label_fontsize(length(sample_labels))
  gene_font <- if (nrow(heatmap_mat) <= 50) 6 else if (nrow(heatmap_mat) <= 100) 4 else 2
  show_rownames <- nrow(heatmap_mat) <= 200

  pdf(file.path(out_sub, paste0("cluster_heatmap_", method_name, ".pdf")),
      width = max(10, ncol(heatmap_mat) * 0.4),
      height = max(8, min(20, nrow(heatmap_mat) * 0.06)))
  pheatmap(
    heatmap_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = anno_row,
    annotation_col = anno_col,
    labels_col = sample_labels[colnames(heatmap_mat)],
    fontsize_col = sample_font,
    fontsize_row = gene_font,
    show_rownames = show_rownames,
    main = paste("Time-series Clustering:", method_name),
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  dev.off()

  png(file.path(out_sub, paste0("cluster_heatmap_", method_name, ".png")),
      width = max(1200, ncol(heatmap_mat) * 80),
      height = max(900, min(2800, nrow(heatmap_mat) * 8)),
      res = 150)
  pheatmap(
    heatmap_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = anno_row,
    annotation_col = anno_col,
    labels_col = sample_labels[colnames(heatmap_mat)],
    fontsize_col = sample_font,
    fontsize_row = gene_font,
    show_rownames = show_rownames,
    main = paste("Time-series Clustering:", method_name),
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  dev.off()

  sample_cluster_mean <- do.call(
    rbind,
    lapply(sort(unique(assignments$cluster)), function(cluster_id) {
      genes_in_cluster <- assignments$gene[assignments$cluster == cluster_id]
      data.frame(
        sample = colnames(z_mat),
        cluster = paste0("Cluster ", cluster_id),
        mean_expression = colMeans(z_mat[genes_in_cluster, , drop = FALSE], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  if (!is.null(time_col)) {
    sample_cluster_mean$time_value <- meta[sample_cluster_mean$sample, time_col]
    if (is.numeric(meta[[time_col]])) {
      sample_cluster_mean$time_numeric <- as.numeric(sample_cluster_mean$time_value)
    } else {
      sample_cluster_mean$time_value <- factor(
        as.character(sample_cluster_mean$time_value),
        levels = unique(as.character(meta[[time_col]]))
      )
    }
  } else {
    sample_cluster_mean$time_value <- factor(
      sample_cluster_mean$sample,
      levels = colnames(z_mat),
      labels = sample_labels[colnames(z_mat)]
    )
  }

  if (!is.null(group_col)) {
    sample_cluster_mean$group_value <- as.factor(meta[sample_cluster_mean$sample, group_col])
  } else {
    sample_cluster_mean$group_value <- factor("all")
  }
  if (!is.null(overlay_col)) {
    sample_cluster_mean$overlay_value <- as.factor(meta[sample_cluster_mean$sample, overlay_col])
  } else {
    sample_cluster_mean$overlay_value <- factor("all")
  }

  mean_labs <- labs(
    title = paste("Cluster Mean Expression:", method_name),
    subtitle = paste0("color=", ifelse(is.null(group_col), "group", group_col),
                      ifelse(is.null(overlay_col), "", paste0(", shape/linetype=", overlay_col))),
    y = "Mean Z-score",
    color = ifelse(is.null(group_col), "group", group_col),
    shape = ifelse(is.null(overlay_col), "", overlay_col),
    linetype = ifelse(is.null(overlay_col), "", overlay_col)
  )

  if (!is.null(time_col) && is.numeric(meta[[time_col]])) {
    p_mean <- ggplot(sample_cluster_mean,
                     aes(x = time_numeric, y = mean_expression,
                         color = group_value,
                         shape = overlay_value,
                         linetype = overlay_value,
                         group = interaction(group_value, overlay_value, cluster))) +
      stat_summary(fun = mean, geom = "line", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 2.2) +
      facet_wrap(~ cluster, scales = "free_y") +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      mean_labs +
      labs(x = time_col)
  } else {
    p_mean <- ggplot(sample_cluster_mean,
                     aes(x = time_value, y = mean_expression,
                         color = group_value,
                         shape = overlay_value,
                         linetype = overlay_value,
                         group = interaction(group_value, overlay_value, cluster))) +
      stat_summary(fun = mean, geom = "line", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 2.2) +
      facet_wrap(~ cluster, scales = "free_y") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      mean_labs +
      labs(x = ifelse(is.null(time_col), "sample", time_col))
  }

  if (is.null(overlay_col)) {
    p_mean <- p_mean + guides(shape = "none", linetype = "none")
  }

  ggsave(file.path(out_sub, paste0("cluster_mean_expression_", method_name, ".pdf")),
         p_mean,
         width = max(10, length(unique(assignments$cluster)) * 3),
         height = max(6, ceiling(length(unique(assignments$cluster)) / 3) * 4))
  ggsave(file.path(out_sub, paste0("cluster_mean_expression_", method_name, ".png")),
         p_mean,
         width = max(10, length(unique(assignments$cluster)) * 3),
         height = max(6, ceiling(length(unique(assignments$cluster)) / 3) * 4),
         dpi = 300)
}

if (cluster_method %in% c("hierarchical", "both")) {
  message("[INFO] 階層クラスタリング実行...")
  d <- make_distance(z_mat, opt$dist_method)
  hc <- hclust(d, method = "ward.D2")
  k_hclust <- safe_k(opt$k, z_mat)
  hc_clusters <- cutree(hc, k = k_hclust)
  hc_assignments <- data.frame(
    gene = hc$labels,
    cluster = as.integer(hc_clusters[hc$labels]),
    stringsAsFactors = FALSE
  )
  hc_assignments <- hc_assignments[order(hc_assignments$cluster, hc_assignments$gene), , drop = FALSE]
  save_cluster_outputs(hc_assignments, "hierarchical")
}

if (cluster_method %in% c("kmeans", "both")) {
  message("[INFO] k-means クラスタリング実行...")
  k_kmeans <- safe_k(opt$k, z_mat)
  set.seed(123)
  km <- kmeans(z_mat, centers = k_kmeans, nstart = 25)
  km_assignments <- data.frame(
    gene = rownames(z_mat),
    cluster = km$cluster,
    stringsAsFactors = FALSE
  )
  km_assignments <- km_assignments[order(km_assignments$cluster, km_assignments$gene), , drop = FALSE]
  save_cluster_outputs(km_assignments, "kmeans")
}

writeLines(
  c(
    paste0("normalization_method=", normalization_method),
    paste0("cluster_method=", cluster_method),
    paste0("top_method=", top_method),
    paste0("time_col=", ifelse(is.null(time_col), "", time_col)),
    paste0("group_col=", ifelse(is.null(group_col), "", group_col)),
    paste0("overlay_col=", ifelse(is.null(overlay_col), "", overlay_col)),
    paste0("genes_used=", nrow(z_mat))
  ),
  con = file.path(out_sub, "clustering_run_info.txt")
)

message("[完了] 05_clustering.R")
