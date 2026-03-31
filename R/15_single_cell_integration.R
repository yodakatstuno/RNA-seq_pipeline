#!/usr/bin/env Rscript
# ============================================================
# 15_single_cell_integration.R
# 疑似bulk / single-cell 統合解析
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
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

option_list <- c(base_option_list, list(
  make_option("--sc_data",     type = "character", default = ""),
  make_option("--method",      type = "character", default = "all"),
  make_option("--resolution",  type = "double",    default = 0.5),
  make_option("--n_neighbors", type = "integer",   default = 30)
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("15_singlecell", opt$outdir)

methods <- if (opt$method == "all") c("umap", "clustering", "marker") else opt$method

# ---- データ読み込み ----
# Single-cell データがある場合はそちらを使用
# ない場合はbulk count matrixをSeuratオブジェクトとして解析

if (opt$sc_data != "" && file.exists(opt$sc_data)) {
  message("[INFO] Single-cell データ読み込み: ", opt$sc_data)
  sc <- readRDS(opt$sc_data)
  if (!inherits(sc, "Seurat")) {
    # count matrixとして扱う
    sc <- CreateSeuratObject(counts = sc, min.cells = 3, min.features = 200)
  }
} else {
  message("[INFO] Bulk count matrixをSeuratオブジェクトとして処理...")
  counts_raw <- load_counts(opt$counts)
  meta_raw   <- load_metadata(opt$metadata)
  dat        <- align_data(counts_raw, meta_raw)

  sc <- CreateSeuratObject(counts = dat$counts, meta.data = dat$metadata)
}

# ---- 前処理 ----
message("[INFO] 正規化 & 変動遺伝子選択...")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
sc <- ScaleData(sc)
sc <- RunPCA(sc, npcs = 50, verbose = FALSE)

# ---- UMAP ----
if ("umap" %in% methods) {
  message("[INFO] UMAP 実行...")
  sc <- RunUMAP(sc, dims = 1:30, n.neighbors = opt$n_neighbors)

  p <- DimPlot(sc, reduction = "umap", pt.size = 1.5) +
    ggtitle("UMAP")
  ggsave(file.path(out_sub, "umap.pdf"), p, width = 10, height = 8)
  ggsave(file.path(out_sub, "umap.png"), p, width = 10, height = 8, dpi = 300)

  # metadata カラムでの色分け
  for (cn in colnames(sc@meta.data)) {
    if (is.character(sc@meta.data[[cn]]) || is.factor(sc@meta.data[[cn]])) {
      n_unique <- length(unique(sc@meta.data[[cn]]))
      if (n_unique > 1 && n_unique < 50) {
        p_col <- DimPlot(sc, reduction = "umap", group.by = cn, pt.size = 1.5) +
          ggtitle(paste("UMAP colored by", cn))
        ggsave(file.path(out_sub, paste0("umap_", cn, ".pdf")), p_col, width = 10, height = 8)
      }
    }
  }
}

# ---- Clustering ----
if ("clustering" %in% methods) {
  message("[INFO] クラスタリング (resolution=", opt$resolution, ")...")
  sc <- FindNeighbors(sc, dims = 1:30)
  sc <- FindClusters(sc, resolution = opt$resolution)

  if ("umap" %in% Reductions(sc)) {
    p_cl <- DimPlot(sc, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.5) +
      ggtitle(paste("Clusters (resolution =", opt$resolution, ")"))
    ggsave(file.path(out_sub, "clusters_umap.pdf"), p_cl, width = 10, height = 8)
    ggsave(file.path(out_sub, "clusters_umap.png"), p_cl, width = 10, height = 8, dpi = 300)
  }

  # クラスタ割り当て保存
  cluster_df <- data.frame(
    cell = colnames(sc),
    cluster = Idents(sc)
  )
  write.csv(cluster_df, file.path(out_sub, "cluster_assignments.csv"), row.names = FALSE)
}

# ---- Marker gene detection ----
if ("marker" %in% methods) {
  message("[INFO] Marker gene 検出...")

  if (!"seurat_clusters" %in% colnames(sc@meta.data)) {
    sc <- FindNeighbors(sc, dims = 1:30)
    sc <- FindClusters(sc, resolution = opt$resolution)
  }

  markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25,
                             logfc.threshold = 0.25)
  write.csv(markers, file.path(out_sub, "marker_genes.csv"), row.names = FALSE)

  # Top markers per cluster
  top_markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 10, order_by = avg_log2FC)
  write.csv(top_markers, file.path(out_sub, "top_markers_per_cluster.csv"), row.names = FALSE)

  # Dot plot
  top5 <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 5, order_by = avg_log2FC)

  if (nrow(top5) > 0) {
    p_dot <- DotPlot(sc, features = unique(top5$gene)) +
      RotatedAxis() +
      ggtitle("Top Marker Genes per Cluster")
    ggsave(file.path(out_sub, "marker_dotplot.pdf"), p_dot,
           width = max(12, length(unique(top5$gene)) * 0.4), height = 8)
  }

  # Heatmap (top markers)
  top3 <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 3, order_by = avg_log2FC)

  if (nrow(top3) > 0) {
    p_heat <- DoHeatmap(sc, features = unique(top3$gene)) +
      ggtitle("Marker Gene Heatmap")
    ggsave(file.path(out_sub, "marker_heatmap.pdf"), p_heat, width = 14, height = 12)
  }
}

# Seurat object 保存
saveRDS(sc, file.path(out_sub, "seurat_object.rds"))

message("[完了] 15_single_cell_integration.R")
