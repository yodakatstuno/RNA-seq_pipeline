#!/usr/bin/env Rscript
# ============================================================
# 07_gsea.R
# Gene Set Enrichment Analysis (GSEA)
# ============================================================
#
# PURPOSE:
#   Ranks ALL genes by log2FoldChange and tests for enrichment
#   of gene sets (GO/KEGG/Reactome) using clusterProfiler.
#   Unlike ORA (Step 6), GSEA does not require a significance cutoff.
#
# INPUTS:
#   --de_result   CSV from Step 1 (must have 'gene' and 'log2FoldChange')
#   --organism    OrgDb package: org.Hs.eg.db | org.Mm.eg.db | org.Dr.eg.db
#   --kegg_org    KEGG code: hsa | mmu | dre
#   --geneset     Which gene sets: GO | KEGG | Reactome | all
#   --min_gs      Minimum gene set size (default: 15)
#   --max_gs      Maximum gene set size (default: 500)
#
# OUTPUTS: GSEA_*_results.csv, dotplots, running score plots (PDF)
#
# NOTE: Reactome is human-only; automatically skipped for other organisms.
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

Sys.setenv(OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           KMP_DUPLICATE_LIB_OK = "TRUE",
           KMP_INIT_AT_FORK = "FALSE",
           KMP_SHM_DISABLE = "1")

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
})

option_list <- c(base_option_list, list(
  make_option("--de_result", type = "character",
              default = file.path(pipeline_output_root, "01_DE", "de_results.csv")),
  make_option("--organism",  type = "character", default = "org.Hs.eg.db"),
  make_option("--kegg_org",  type = "character", default = "hsa"),
  make_option("--geneset",   type = "character", default = "all"),
  make_option("--ranking_metric", type = "character", default = "log2FoldChange"),
  make_option("--min_gs",    type = "integer",   default = 15),
  make_option("--max_gs",    type = "integer",   default = 500)
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("07_gsea", opt$outdir)

suppressPackageStartupMessages({
  if (!require(opt$organism, character.only = TRUE)) {
    stop("[ERROR] OrgDb パッケージを読み込めません: ", opt$organism)
  }
})
orgdb <- get(opt$organism)

strip_version <- function(x) sub("\\\\..*$", "", x)
detect_gene_id_type <- function(genes) {
  g <- strip_version(genes)
  ens_frac <- mean(grepl("^ENS", g), na.rm = TRUE)
  if (is.nan(ens_frac)) ens_frac <- 0
  if (ens_frac >= 0.6) {
    return("ENSEMBL")
  } else if (ens_frac > 0) {
    return("ENSEMBL")
  } else {
    return("SYMBOL")
  }
}

map_de_results <- function(df, orgdb) {
  df$gene_clean <- strip_version(df$gene)
  id_type <- detect_gene_id_type(df$gene_clean)
  mapping <- suppressWarnings(bitr(unique(df$gene_clean),
                                   fromType = id_type,
                                   toType = "ENTREZID",
                                   OrgDb = orgdb))
  mapping <- unique(na.omit(mapping))
  rate <- ifelse(length(unique(df$gene_clean)) > 0,
                 nrow(mapping) / length(unique(df$gene_clean)), 0)
  if (rate < 0.6) {
    message("[WARNING] Gene ID mapping rate is low: ",
            round(rate * 100, 1), "% (fromType=", id_type, ")")
  } else {
    message("[INFO] Gene ID mapping rate: ",
            round(rate * 100, 1), "% (fromType=", id_type, ")")
  }
  df_mapped <- merge(df, mapping, by.x = "gene_clean", by.y = id_type)
  df_mapped <- df_mapped[!is.na(df_mapped$ENTREZID), ]
  list(df = df_mapped, id_type = id_type, rate = rate)
}

normalize_ranking_metric <- function(x) {
  metric <- trimws(as.character(x))
  allowed <- c("log2FoldChange", "stat", "custom")
  if (!(metric %in% allowed)) {
    stop("[ERROR] --ranking_metric は log2FoldChange / stat / custom のいずれかです。")
  }
  metric
}

build_ranking_metric <- function(df, metric) {
  if (metric == "log2FoldChange") {
    if (!("log2FoldChange" %in% colnames(df))) {
      stop("[ERROR] DE結果に log2FoldChange 列がありません。")
    }
    vals <- as.numeric(df$log2FoldChange)
  } else if (metric == "stat") {
    if (!("stat" %in% colnames(df))) {
      stop("[ERROR] DE結果に stat 列がありません。")
    }
    vals <- as.numeric(df$stat)
  } else {
    required <- c("pvalue", "log2FoldChange")
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0) {
      stop("[ERROR] custom ranking には以下の列が必要です: ",
           paste(missing_cols, collapse = ", "))
    }
    pvals <- as.numeric(df$pvalue)
    lfc <- as.numeric(df$log2FoldChange)
    pvals[is.na(pvals) | pvals <= 0] <- .Machine$double.xmin
    vals <- -log10(pvals) * sign(lfc)
  }
  vals[!is.finite(vals)] <- NA_real_
  vals
}

# ---- DE結果 → ランキング ----
de <- read.csv(opt$de_result, stringsAsFactors = FALSE)
opt$ranking_metric <- normalize_ranking_metric(opt$ranking_metric)

mapped <- map_de_results(de, orgdb)
de_mapped <- mapped$df
message("[INFO] 検出された gene ID type: ", mapped$id_type)
if (nrow(de_mapped) == 0) {
  message("[WARNING] DE結果のID変換に失敗しました。遺伝子IDを確認してください。")
  quit(status = 0)
}

# ランキング
de_mapped$ranking_metric_value <- build_ranking_metric(de_mapped, opt$ranking_metric)
agg <- stats::aggregate(ranking_metric_value ~ ENTREZID, data = de_mapped,
                        FUN = function(x) {
                          x <- as.numeric(x)
                          x <- x[is.finite(x)]
                          if (length(x) == 0) return(NA_real_)
                          mean(x, na.rm = TRUE)
                        })
gene_list <- agg$ranking_metric_value
names(gene_list) <- agg$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[is.finite(gene_list)]
if (length(gene_list) == 0) {
  message("[WARNING] GSEA 用のランク付け可能な遺伝子がありません。IDや ranking_metric を確認してください。")
  quit(status = 0)
}

message("[INFO] Ranking metric: ", opt$ranking_metric)
message("[INFO] GSEA ランキング遺伝子数: ", length(gene_list))

genesets <- if (opt$geneset == "all") c("GO", "KEGG", "Reactome") else opt$geneset

# ---- GSEA: GO ----
if ("GO" %in% genesets) {
  message("[GSEA] GO解析...")
  gsea_go <- gseGO(geneList = gene_list,
                    OrgDb = orgdb,
                    ont = "BP",
                    minGSSize = opt$min_gs,
                    maxGSSize = opt$max_gs,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
  if (nrow(as.data.frame(gsea_go)) > 0) {
    write.csv(as.data.frame(gsea_go), file.path(out_sub, "GSEA_GO_results.csv"), row.names = FALSE)

    p <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
      facet_grid(~.sign) + ggtitle("GSEA GO BP")
    ggsave(file.path(out_sub, "GSEA_GO_dotplot.pdf"), p, width = 14, height = 10)

    # 上位のランニングスコアプロット
    top_ids <- head(as.data.frame(gsea_go)$ID, 4)
    for (i in seq_along(top_ids)) {
      p_run <- gseaplot2(gsea_go, geneSetID = top_ids[i], title = top_ids[i])
      ggsave(file.path(out_sub, paste0("GSEA_GO_running_", i, ".pdf")), p_run,
             width = 8, height = 6)
    }
  } else {
    message("[INFO] GSEA GO: 有意な結果なし")
  }
}

# ---- GSEA: KEGG ----
if ("KEGG" %in% genesets) {
  message("[GSEA] KEGG解析...")
  gsea_kegg <- gseKEGG(geneList = gene_list,
                        organism = opt$kegg_org,
                        minGSSize = opt$min_gs,
                        maxGSSize = opt$max_gs,
                        pvalueCutoff = 0.05,
                        verbose = FALSE)
  if (nrow(as.data.frame(gsea_kegg)) > 0) {
    write.csv(as.data.frame(gsea_kegg), file.path(out_sub, "GSEA_KEGG_results.csv"),
              row.names = FALSE)
    p <- dotplot(gsea_kegg, showCategory = 20) + ggtitle("GSEA KEGG")
    ggsave(file.path(out_sub, "GSEA_KEGG_dotplot.pdf"), p, width = 12, height = 10)
  }
}

# ---- GSEA: Reactome ----
if ("Reactome" %in% genesets) {
  if (tolower(opt$organism) != "org.hs.eg.db") {
    message("[INFO] Reactome はヒト専用のためスキップ (organism=", opt$organism, ")")
  } else {
    tryCatch({
      suppressPackageStartupMessages(library(ReactomePA))
      message("[GSEA] Reactome解析...")
      gsea_react <- gsePathway(geneList = gene_list,
                                minGSSize = opt$min_gs,
                                maxGSSize = opt$max_gs,
                                pvalueCutoff = 0.05,
                                verbose = FALSE)
      if (nrow(as.data.frame(gsea_react)) > 0) {
        write.csv(as.data.frame(gsea_react), file.path(out_sub, "GSEA_Reactome_results.csv"),
                  row.names = FALSE)
        p <- dotplot(gsea_react, showCategory = 20) + ggtitle("GSEA Reactome")
        ggsave(file.path(out_sub, "GSEA_Reactome_dotplot.pdf"), p, width = 12, height = 10)
      } else {
        message("[INFO] GSEA Reactome: 有意な結果なし")
      }
    }, error = function(e) message("[WARNING] ReactomePA GSEA failed: ", e$message))
  }
}

message("[完了] 07_gsea.R")
