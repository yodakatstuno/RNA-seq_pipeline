#!/usr/bin/env Rscript
# ============================================================
# 14_network_analysis.R
# 遺伝子間ネットワーク解析
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
  library(ggplot2)
  library(igraph)
  library(pheatmap)
})

option_list <- c(base_option_list, list(
  make_option("--method",        type = "character", default = "both"),
  make_option("--top_n",         type = "integer",   default = 500),
  make_option("--cor_threshold", type = "double",    default = 0.8),
  make_option("--hub_top",       type = "integer",   default = 20),
  make_option("--detect_modules", type = "character", default = "FALSE"),
  make_option("--trait_cols",    type = "character", default = "")
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("14_network", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata

log_counts <- log2(counts + 1)
filtered <- select_top_variable(log_counts, opt$top_n)

methods <- if (opt$method == "both") c("coexpression", "wgcna") else opt$method
parse_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}
parse_csv_arg <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",")))
  vals[vals != ""]
}
opt$detect_modules <- parse_bool(opt$detect_modules)
trait_cols <- parse_csv_arg(opt$trait_cols)
if (opt$detect_modules && !("wgcna" %in% methods)) {
  methods <- unique(c(methods, "wgcna"))
}

# ---- Co-expression network ----
if ("coexpression" %in% methods) {
  message("[NET] Co-expression network 構築...")

  cor_mat <- cor(t(filtered), method = "pearson")

  # 閾値でエッジ作成
  adj_mat <- abs(cor_mat)
  adj_mat[adj_mat < opt$cor_threshold] <- 0
  diag(adj_mat) <- 0

  # igraph
  g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
  g <- simplify(g)

  message("[NET] ノード数: ", vcount(g), " エッジ数: ", ecount(g))

  # Hub gene (degree centrality)
  deg <- degree(g)
  hub_df <- data.frame(gene = names(deg), degree = deg)
  hub_df <- hub_df[order(-hub_df$degree), ]
  write.csv(hub_df, file.path(out_sub, "hub_genes_degree.csv"), row.names = FALSE)

  # Betweenness centrality
  btw <- betweenness(g)
  hub_df$betweenness <- btw[hub_df$gene]
  write.csv(hub_df, file.path(out_sub, "network_centrality.csv"), row.names = FALSE)

  # Top hub genes
  top_hubs <- head(hub_df, opt$hub_top)
  p_hub <- ggplot(top_hubs, aes(x = reorder(gene, degree), y = degree)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top", opt$hub_top, "Hub Genes (Degree)"), x = "", y = "Degree")
  ggsave(file.path(out_sub, "hub_genes_barplot.pdf"), p_hub, width = 8, height = 8)
  ggsave(file.path(out_sub, "hub_genes_barplot.png"), p_hub, width = 8, height = 8, dpi = 300)

  # ネットワーク可視化 (上位ノードのサブグラフ)
  if (vcount(g) > 0 && ecount(g) > 0) {
    top_nodes <- head(hub_df$gene, min(100, vcount(g)))
    sub_g <- induced_subgraph(g, top_nodes)

    pdf(file.path(out_sub, "coexpression_network.pdf"), width = 14, height = 14)
    plot(sub_g,
         vertex.size = sqrt(degree(sub_g)) * 2,
         vertex.label.cex = 0.5,
         vertex.color = "lightblue",
         edge.width = 0.5,
         main = "Co-expression Network (Top Genes)")
    dev.off()
  }

  # エッジリスト保存
  el <- as_data_frame(g, what = "edges")
  write.csv(el, file.path(out_sub, "edge_list.csv"), row.names = FALSE)
}

# ---- WGCNA-based network ----
if ("wgcna" %in% methods) {
  suppressPackageStartupMessages(library(WGCNA))
  if (exists("disableWGCNAThreads")) {
    disableWGCNAThreads()
  }
  message("[NET] WGCNA ネットワーク解析...")

  datExpr <- t(filtered)
  sft <- pickSoftThreshold(datExpr, powerVector = c(1:20), verbose = 0)
  power <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

  # TOM
  adj <- adjacency(datExpr, power = power)
  TOM <- TOMsimilarity(adj)
  rownames(TOM) <- colnames(TOM) <- colnames(datExpr)

  # Connectivity
  connectivity <- rowSums(adj) - 1
  conn_df <- data.frame(gene = names(connectivity), connectivity = connectivity)
  conn_df <- conn_df[order(-conn_df$connectivity), ]
  write.csv(conn_df, file.path(out_sub, "wgcna_connectivity.csv"), row.names = FALSE)

  top_conn <- head(conn_df, opt$hub_top)
  p_conn <- ggplot(top_conn, aes(x = reorder(gene, connectivity), y = connectivity)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top", opt$hub_top, "Hub Genes (WGCNA Connectivity)"),
         x = "", y = "Connectivity")
  ggsave(file.path(out_sub, "wgcna_hub_genes.pdf"), p_conn, width = 8, height = 8)

  if (opt$detect_modules) {
    message("[NET] WGCNA module detection 実行...")
    net <- blockwiseModules(
      datExpr,
      power = power,
      networkType = "signed",
      TOMType = "signed",
      minModuleSize = min(30, ncol(datExpr)),
      reassignThreshold = 0,
      mergeCutHeight = 0.25,
      numericLabels = FALSE,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    )

    module_df <- data.frame(
      gene = colnames(datExpr),
      module = net$colors,
      stringsAsFactors = FALSE
    )
    write.csv(module_df, file.path(out_sub, "wgcna_module_assignments.csv"), row.names = FALSE)

    module_counts <- sort(table(net$colors), decreasing = TRUE)
    write.csv(data.frame(module = names(module_counts), n_genes = as.integer(module_counts)),
              file.path(out_sub, "wgcna_module_sizes.csv"), row.names = FALSE)

    if (length(trait_cols) > 0) {
      missing_traits <- setdiff(trait_cols, colnames(meta))
      if (length(missing_traits) > 0) {
        stop("[ERROR] trait_cols が metadata に存在しません: ",
             paste(missing_traits, collapse = ", "))
      }

      trait_df <- meta[, trait_cols, drop = FALSE]
      trait_numeric_list <- list()
      for (cn in trait_cols) {
        vals <- trait_df[[cn]]
        if (is.numeric(vals)) {
          trait_numeric_list[[cn]] <- vals
        } else {
          mm <- model.matrix(~ as.factor(vals) - 1)
          colnames(mm) <- paste(cn, sub("^as.factor\\(vals\\)", "", colnames(mm)), sep = "_")
          trait_numeric_list[[cn]] <- mm
        }
      }
      trait_matrix <- as.data.frame(do.call(cbind, trait_numeric_list))
      rownames(trait_matrix) <- rownames(meta)

      MEs <- orderMEs(net$MEs)
      if (ncol(MEs) == 0) {
        message("[WARNING] モジュール eigengene が計算できなかったため、module-trait 相関をスキップします。")
      } else {
      cor_mat <- cor(MEs, trait_matrix, use = "pairwise.complete.obs")
      p_mat <- corPvalueStudent(cor_mat, nSamples = nrow(datExpr))

      write.csv(cor_mat, file.path(out_sub, "wgcna_module_trait_correlations.csv"))
      write.csv(p_mat, file.path(out_sub, "wgcna_module_trait_pvalues.csv"))

      display_mat <- matrix(
        paste0(sprintf("%.2f", cor_mat), "\n(p=", sprintf("%.3g", p_mat), ")"),
        nrow = nrow(cor_mat),
        dimnames = dimnames(cor_mat)
      )

      pdf(file.path(out_sub, "wgcna_module_trait_heatmap.pdf"), width = 10, height = 8)
      pheatmap(
        cor_mat,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        display_numbers = display_mat,
        cluster_rows = nrow(cor_mat) > 1,
        cluster_cols = ncol(cor_mat) > 1,
        main = "Module-Trait Correlation"
      )
      dev.off()
      }
    }
  }
}

message("[完了] 14_network_analysis.R")
