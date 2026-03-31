#!/usr/bin/env Rscript
# ============================================================
# 06_functional_enrichment.R
# GO / KEGG / Reactome Functional Enrichment
# ============================================================
#
# PURPOSE:
#   Performs over-representation analysis (ORA) on DE genes:
#   GO (BP/MF/CC), KEGG, and Reactome (human only).
#
# INPUTS:
#   --de_result   CSV from Step 1 (must have 'gene' and 'significant' columns)
#   --organism    OrgDb package: org.Hs.eg.db | org.Mm.eg.db | org.Dr.eg.db
#   --kegg_org    KEGG code: hsa | mmu | dre
#   --pval_cutoff P-value threshold for enrichment (default: 0.05)
#
# GENE ID HANDLING:
#   Auto-detects ENSEMBL vs SYMBOL from gene names.
#   Maps to ENTREZID internally via bitr().
#   Logs mapping rate — low rates indicate ID mismatch.
#
# OUTPUTS: GO_*_results.csv, KEGG_results.csv, Reactome_results.csv,
#          dotplots, barplots, emapplots (PDF/PNG)
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
})

option_list <- c(base_option_list, list(
  make_option("--de_result",   type = "character",
              default = file.path(pipeline_output_root, "01_DE", "de_results.csv")),
  make_option("--organism",    type = "character", default = "org.Hs.eg.db"),
  make_option("--kegg_org",    type = "character", default = "hsa"),
  make_option("--pval_cutoff", type = "double",    default = 0.05)
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("06_enrichment", opt$outdir)

# OrgDb ロード（安全に）
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

map_to_entrez <- function(genes, orgdb) {
  clean <- strip_version(genes)
  id_type <- detect_gene_id_type(clean)
  mapping <- suppressWarnings(bitr(clean,
                                   fromType = id_type,
                                   toType = "ENTREZID",
                                   OrgDb = orgdb))
  mapping <- unique(na.omit(mapping))
  uniq_in <- length(unique(clean))
  rate <- ifelse(uniq_in > 0, nrow(mapping) / uniq_in, 0)
  if (rate < 0.6) {
    message("[WARNING] Gene ID mapping rate is low: ",
            round(rate * 100, 1), "% (fromType=", id_type, ")")
  } else {
    message("[INFO] Gene ID mapping rate: ",
            round(rate * 100, 1), "% (fromType=", id_type, ")")
  }
  list(entrez = unique(mapping$ENTREZID), id_type = id_type, rate = rate)
}

# ---- DE結果読み込み ----
message("[INFO] DE結果読み込み: ", opt$de_result)
de <- read.csv(opt$de_result, stringsAsFactors = FALSE)

sig_genes <- de$gene[de$significant == TRUE]
all_genes <- de$gene
message("[INFO] 有意な遺伝子数: ", length(sig_genes))

if (length(sig_genes) == 0) {
  message("[WARNING] 有意な遺伝子が0個です。閾値を確認してください。")
  quit(status = 0)
}

# ---- Entrez ID 変換 ----
sig_map <- map_to_entrez(sig_genes, orgdb)
bg_map  <- map_to_entrez(all_genes, orgdb)
message("[INFO] 検出された gene ID type: ", sig_map$id_type)

if (length(sig_map$entrez) == 0) {
  message("[WARNING] 有意遺伝子の Entrez 変換に失敗しました。遺伝子IDを確認してください。")
  quit(status = 0)
}

# ---- GO解析 ----
message("[INFO] GO enrichment 解析...")

for (ont in c("BP", "MF", "CC")) {
  ego <- enrichGO(gene = sig_map$entrez,
                   universe = bg_map$entrez,
                   OrgDb = orgdb,
                   ont = ont,
                   pAdjustMethod = "BH",
                   pvalueCutoff = opt$pval_cutoff,
                   readable = TRUE)

  if (nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(out_sub, paste0("GO_", ont, "_results.csv")),
              row.names = FALSE)

    p_dot <- dotplot(ego, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment"))
    ggsave(file.path(out_sub, paste0("GO_", ont, "_dotplot.pdf")), p_dot, width = 10, height = 8)
    ggsave(file.path(out_sub, paste0("GO_", ont, "_dotplot.png")), p_dot, width = 10, height = 8, dpi = 300)

    p_bar <- barplot(ego, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment"))
    ggsave(file.path(out_sub, paste0("GO_", ont, "_barplot.pdf")), p_bar, width = 10, height = 8)

    if (nrow(as.data.frame(ego)) >= 2) {
      tryCatch({
        p_emap <- emapplot(pairwise_termsim(ego), showCategory = 30)
        ggsave(file.path(out_sub, paste0("GO_", ont, "_emap.pdf")), p_emap, width = 12, height = 10)
      }, error = function(e) message("[WARNING] emap plot failed for GO ", ont))
    }
  } else {
    message("[INFO] GO ", ont, ": 有意なtermなし")
  }
}

# ---- KEGG解析 ----
message("[INFO] KEGG enrichment 解析...")
ekegg <- enrichKEGG(gene = sig_map$entrez,
                     universe = bg_map$entrez,
                     organism = opt$kegg_org,
                     pvalueCutoff = opt$pval_cutoff)

if (nrow(as.data.frame(ekegg)) > 0) {
  write.csv(as.data.frame(ekegg), file.path(out_sub, "KEGG_results.csv"), row.names = FALSE)
  p_kegg <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
  ggsave(file.path(out_sub, "KEGG_dotplot.pdf"), p_kegg, width = 10, height = 8)
  ggsave(file.path(out_sub, "KEGG_dotplot.png"), p_kegg, width = 10, height = 8, dpi = 300)
} else {
  message("[INFO] KEGG: 有意なpathwayなし")
}

# ---- Reactome解析 ----
if (tolower(opt$organism) == "org.hs.eg.db") {
  tryCatch({
    suppressPackageStartupMessages(library(ReactomePA))
    message("[INFO] Reactome enrichment 解析...")
    ereactome <- enrichPathway(gene = sig_map$entrez,
                                universe = bg_map$entrez,
                                pvalueCutoff = opt$pval_cutoff,
                                readable = TRUE)
    if (nrow(as.data.frame(ereactome)) > 0) {
      write.csv(as.data.frame(ereactome), file.path(out_sub, "Reactome_results.csv"),
                row.names = FALSE)
      p_react <- dotplot(ereactome, showCategory = 20) + ggtitle("Reactome Pathway Enrichment")
      ggsave(file.path(out_sub, "Reactome_dotplot.pdf"), p_react, width = 10, height = 8)
    } else {
      message("[INFO] Reactome: 有意なpathwayなし")
    }
  }, error = function(e) {
    message("[WARNING] ReactomePA が利用できません: ", e$message)
  })
} else {
  message("[INFO] Reactome はヒト専用のためスキップ (organism=", opt$organism, ")")
}

message("[完了] 06_functional_enrichment.R")
