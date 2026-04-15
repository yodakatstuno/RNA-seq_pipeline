#!/usr/bin/env Rscript
# ============================================================
# 00_load_data.R — Common data loading utilities
# 共通データ読み込みユーティリティ
# ============================================================
#
# PURPOSE:
#   Shared utility functions sourced by all analysis scripts.
#   Provides count matrix loading, metadata loading, and
#   sample alignment between counts and metadata.
#
# FUNCTIONS:
#   load_counts(path)    - Load .rds count matrix (genes × samples)
#   load_metadata(path)  - Load .xlsx/.csv/.tsv metadata
#   align_data(counts, meta) - Match samples between counts and metadata
#   filter_low_counts()  - Remove low-expression genes
#   select_top_variable() - Select top variable genes
#
# METADATA FORMAT:
#   Must contain a column whose values match count matrix column names.
#   Auto-detected: "ID", "sample", "sample_id", "sampleID", or best match.
#   Other columns (condition, genotype, dpf, batch, etc.) are preserved.
#
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
})

PIPELINE_STEP_DIRS <- c(
  "01_DE", "02_clustering", "03_dimreduc", "04_normalization",
  "05_clustering", "05_expression_pattern", "06_enrichment", "07_gsea",
  "08_batch_correction", "09_regression", "10_qc",
  "11_visualization", "12_phenotype", "13_ml",
  "14_network", "15_singlecell"
)

pipeline_root_dir <- local({
  if (exists("script_dir", inherits = TRUE)) {
    normalizePath(file.path(get("script_dir", inherits = TRUE), ".."), mustWork = FALSE)
  } else {
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  }
})

pipeline_output_root <- local({
  configured_root <- Sys.getenv("PIPELINE_OUTPUT_DIR", unset = "")
  if (nzchar(configured_root)) {
    normalizePath(configured_root, mustWork = FALSE)
  } else {
    normalizePath(path.expand("~/Output"), mustWork = FALSE)
  }
})

#' コマンドライン引数をパースする基本オプションリスト
base_option_list <- list(
  make_option("--counts",   type = "character", help = "Count matrix RDS path"),
  make_option("--metadata", type = "character", help = "Metadata xlsx/csv path"),
  make_option("--outdir",   type = "character", default = "", help = "Output directory (retained for compatibility; results are written under ~/Output or PIPELINE_OUTPUT_DIR)"),
  make_option("--label_mode", type = "character", default = "sampleID_only",
              help = "Sample label display: sampleID_only | symbol_only | both"),
  make_option("--symbol_col", type = "character", default = "",
              help = "Metadata column used when label_mode is symbol_only or both")
)

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

normalize_label_mode <- function(x) {
  mode <- trimws(as.character(x))
  allowed <- c("sampleID_only", "symbol_only", "both")
  if (!(mode %in% allowed)) {
    stop("[ERROR] --label_mode は sampleID_only / symbol_only / both のいずれかです。")
  }
  mode
}

ensure_pipeline_output_root <- function(requested_outdir = "") {
  requested_outdir <- trimws(as.character(requested_outdir))
  effective_root <- pipeline_output_root
  if (!identical(requested_outdir, "")) {
    effective_root <- normalizePath(requested_outdir, mustWork = FALSE)
  }
  dir.create(effective_root, recursive = TRUE, showWarnings = FALSE)
  if (!identical(requested_outdir, "") &&
      !identical(effective_root, pipeline_output_root)) {
    message("[INFO] --outdir='", requested_outdir,
            "' を使用します (PIPELINE_OUTPUT_DIR または ~/Output を上書き)。")
  }
  assign("pipeline_output_root", effective_root, envir = parent.env(environment()))
  effective_root
}

resolve_step_outdir <- function(step_dir, requested_outdir = "", analysis_name = NULL) {
  if (!(step_dir %in% PIPELINE_STEP_DIRS)) {
    stop("[ERROR] 未知の step_dir: ", step_dir)
  }
  base_dir <- file.path(ensure_pipeline_output_root(requested_outdir), step_dir)
  out_dir <- base_dir
  if (!is.null(analysis_name) && !is.na(analysis_name) && analysis_name != "" &&
      analysis_name != step_dir) {
    out_dir <- file.path(base_dir, analysis_name)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

validate_label_options <- function(meta, label_mode = "sampleID_only", symbol_col = "") {
  label_mode <- normalize_label_mode(label_mode)
  if (label_mode %in% c("symbol_only", "both")) {
    if (is.null(symbol_col) || is.na(symbol_col) || trimws(symbol_col) == "") {
      stop("[ERROR] --label_mode=", label_mode,
           " の場合は --symbol_col が必須です。")
    }
    if (!(symbol_col %in% colnames(meta))) {
      stop("[ERROR] symbol_col が metadata に存在しません: ", symbol_col)
    }
    symbol_vals <- trimws(as.character(meta[[symbol_col]]))
    if (any(is.na(symbol_vals) | symbol_vals == "")) {
      stop("[ERROR] symbol_col='", symbol_col,
           "' に空欄または NA が含まれています。")
    }
  }
  invisible(TRUE)
}

build_sample_labels <- function(meta, label_mode = "sampleID_only", symbol_col = "") {
  label_mode <- normalize_label_mode(label_mode)
  validate_label_options(meta, label_mode, symbol_col)
  sample_ids <- rownames(meta)
  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("[ERROR] metadata の行名にサンプルIDが設定されていません。")
  }

  if (label_mode == "sampleID_only") {
    labels <- sample_ids
  } else {
    symbols <- trimws(as.character(meta[[symbol_col]]))
    if (label_mode == "symbol_only") {
      labels <- symbols
    } else {
      labels <- paste0(sample_ids, " | ", symbols)
    }
  }
  names(labels) <- sample_ids
  labels
}

sample_label_cex <- function(n_labels) {
  if (n_labels <= 12) return(0.9)
  if (n_labels <= 24) return(0.8)
  if (n_labels <= 40) return(0.7)
  if (n_labels <= 60) return(0.6)
  0.5
}

sample_label_size <- function(n_labels) {
  if (n_labels <= 12) return(3.0)
  if (n_labels <= 24) return(2.7)
  if (n_labels <= 40) return(2.4)
  if (n_labels <= 60) return(2.1)
  1.8
}

sample_label_fontsize <- function(n_labels) {
  if (n_labels <= 12) return(10)
  if (n_labels <= 24) return(9)
  if (n_labels <= 40) return(8)
  if (n_labels <= 60) return(7)
  6
}

#' count matrix を読み込む
load_counts <- function(path) {
  message("[INFO] Count matrix 読み込み: ", path)
  counts <- readRDS(path)
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  message("[INFO] Count matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")
  return(counts)
}

#' metadata を読み込む
load_metadata <- function(path) {
  message("[INFO] Metadata 読み込み: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    meta <- as.data.frame(readxl::read_excel(path))
  } else if (ext == "csv") {
    meta <- read.csv(path, stringsAsFactors = FALSE)
  } else if (ext == "tsv" || ext == "txt") {
    meta <- read.delim(path, stringsAsFactors = FALSE)
  } else {
    stop("対応していないメタデータ形式: ", ext)
  }
  message("[INFO] Metadata: ", nrow(meta), " samples x ", ncol(meta), " columns")
  message("[INFO] Metadata columns: ", paste(colnames(meta), collapse = ", "))
  return(meta)
}

#' count matrix と metadata のサンプル名を整合させる
align_data <- function(counts, meta, sample_col = NULL) {
  # metadata のサンプルIDカラムを推定
  if (is.null(sample_col)) {
    # 優先候補（完全一致 + 実際にマッチする列を優先）
    priority <- c("ID", "sample", "sample_id", "sampleID")
    for (cn in priority) {
      if (cn %in% colnames(meta)) {
        overlap <- sum(as.character(meta[[cn]]) %in% colnames(counts))
        message("[DEBUG] Checking column: ", cn, " → matches: ", overlap)
        if (overlap > 0) {
          sample_col <- cn
          break
        }
      }
    }

    # フォールバック: 文字列候補から最大一致の列を選択
    if (is.null(sample_col)) {
      candidates <- grep("sample|Sample|SAMPLE|SampleID|sample_id",
                         colnames(meta), value = TRUE)
      if (length(candidates) > 0) {
        overlaps <- sapply(candidates, function(cn) {
          sum(as.character(meta[[cn]]) %in% colnames(counts))
        })
        for (i in seq_along(candidates)) {
          message("[DEBUG] Checking column: ", candidates[i], " → matches: ", overlaps[i])
        }
        if (max(overlaps) > 0) {
          sample_col <- candidates[which.max(overlaps)]
        }
      }
    }

    # 最終フォールバック: 全列から最大一致を選択（0一致は不可）
    if (is.null(sample_col)) {
      overlaps <- sapply(colnames(meta), function(cn) {
        sum(as.character(meta[[cn]]) %in% colnames(counts))
      })
      for (i in seq_along(overlaps)) {
        message("[DEBUG] Checking column: ", colnames(meta)[i], " → matches: ", overlaps[i])
      }
      if (max(overlaps) > 0) {
        sample_col <- colnames(meta)[which.max(overlaps)]
      }
    }
  }

  if (is.null(sample_col)) {
    stop("[ERROR] サンプルID列を特定できませんでした（一致する列がありません）。")
  }

  message("[INFO] サンプルIDカラム: ", sample_col)

  if (any(is.na(colnames(counts)))) {
    stop("[ERROR] Count matrix のサンプル名に NA が含まれています。")
  }
  if (any(colnames(counts) == "")) {
    stop("[ERROR] Count matrix のサンプル名に空文字が含まれています。")
  }
  if (any(duplicated(colnames(counts)))) {
    stop("[ERROR] Count matrix に重複したサンプルIDが存在します。")
  }

  if (any(is.na(meta[[sample_col]]))) {
    stop("[ERROR] Metadata のサンプルIDに NA が含まれています: ", sample_col)
  }
  if (any(meta[[sample_col]] == "")) {
    stop("[ERROR] Metadata のサンプルIDに空文字が含まれています: ", sample_col)
  }
  rownames(meta) <- meta[[sample_col]]
  if (any(duplicated(rownames(meta)))) {
    stop("[ERROR] Metadata に重複したサンプルIDが存在します。")
  }
  common_samples <- intersect(colnames(counts), rownames(meta))

  if (length(common_samples) == 0) {
    stop("[ERROR] Count matrix と Metadata の間に共通サンプルがありません。\n",
         "  Count matrix サンプル (先頭5): ",
         paste(head(colnames(counts), 5), collapse = ", "), "\n",
         "  Metadata サンプル (先頭5): ",
         paste(head(rownames(meta), 5), collapse = ", "))
  }

  missing_in_counts <- setdiff(rownames(meta), colnames(counts))
  missing_in_meta <- setdiff(colnames(counts), rownames(meta))

  message("[INFO] 共通サンプル数: ", length(common_samples))
  message("[INFO] Count matrix にのみ存在するサンプル数: ", length(missing_in_meta))
  message("[INFO] Metadata にのみ存在するサンプル数: ", length(missing_in_counts))

  if (length(missing_in_meta) > 0) {
    message("[INFO] Count matrix のみ (先頭5): ",
            paste(head(missing_in_meta, 5), collapse = ", "))
  }
  if (length(missing_in_counts) > 0) {
    message("[INFO] Metadata のみ (先頭5): ",
            paste(head(missing_in_counts, 5), collapse = ", "))
  }

  counts <- counts[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]

  if (!all(colnames(counts) == rownames(meta))) {
    stop("[ERROR] counts と metadata のサンプル順序が一致していません。")
  }

  return(list(counts = counts, metadata = meta))
}

#' 基本的なフィルタリング
filter_low_counts <- function(counts, min_count = 10, min_samples = 2) {
  keep <- rowSums(counts >= min_count) >= min_samples
  message("[INFO] フィルタリング: ", sum(keep), " / ", nrow(counts),
          " genes retained (min_count=", min_count,
          ", min_samples=", min_samples, ")")
  return(counts[keep, , drop = FALSE])
}

#' 上位変動遺伝子を選択
select_top_variable <- function(counts, top_n = 2000) {
  vars <- apply(counts, 1, var)
  top_genes <- names(sort(vars, decreasing = TRUE))[1:min(top_n, length(vars))]
  message("[INFO] 上位 ", length(top_genes), " 変動遺伝子を選択")
  return(counts[top_genes, , drop = FALSE])
}

infer_transform_design <- function(meta) {
  if (all(c("genotype", "dpf") %in% colnames(meta))) {
    message("[INFO] DESeq2 transform design: ~ genotype + dpf")
    return(stats::as.formula("~ genotype + dpf"))
  }
  message("[INFO] DESeq2 transform design: ~ 1")
  stats::as.formula("~ 1")
}

prepare_transform_metadata <- function(meta) {
  meta_prepared <- meta
  if (!is.null(rownames(meta_prepared)) && any(rownames(meta_prepared) == "")) {
    stop("[ERROR] metadata の行名にサンプルIDが必要です。")
  }
  for (cn in colnames(meta_prepared)) {
    if (cn == "dpf" && !is.factor(meta_prepared[[cn]])) {
      meta_prepared[[cn]] <- as.factor(meta_prepared[[cn]])
      next
    }
    if (is.character(meta_prepared[[cn]]) || is.logical(meta_prepared[[cn]])) {
      meta_prepared[[cn]] <- as.factor(meta_prepared[[cn]])
    }
  }
  meta_prepared
}

create_transform_dds <- function(counts, meta) {
  suppressPackageStartupMessages(library(DESeq2))
  counts_int <- round(as.matrix(counts))
  storage.mode(counts_int) <- "integer"
  meta_prepared <- prepare_transform_metadata(meta)
  design_formula <- infer_transform_design(meta_prepared)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_int,
    colData = meta_prepared,
    design = design_formula
  )
  DESeq2::estimateSizeFactors(dds)
}

row_zscore_matrix <- function(mat, center = TRUE, scale = TRUE) {
  z_mat <- t(scale(t(as.matrix(mat)), center = center, scale = scale))
  z_mat[is.na(z_mat)] <- 0
  z_mat
}

normalize_counts_matrix <- function(counts, meta, method,
                                    force_rlog = FALSE,
                                    rlog_max_samples = 30) {
  method <- tolower(trimws(as.character(method)))

  if (method == "deseq2") {
    dds <- create_transform_dds(counts, meta)
    return(list(
      matrix = DESeq2::counts(dds, normalized = TRUE),
      extras = list(size_factors = DESeq2::sizeFactors(dds))
    ))
  }

  if (method %in% c("tmm", "cpm")) {
    suppressPackageStartupMessages(library(edgeR))
    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    return(list(
      matrix = edgeR::cpm(y, normalized.lib.sizes = TRUE),
      extras = list(norm_factors = y$samples)
    ))
  }

  if (method == "log2") {
    return(list(matrix = log2(as.matrix(counts) + 1), extras = list()))
  }

  if (method == "vst") {
    dds <- create_transform_dds(counts, meta)
    vsd <- tryCatch(
      DESeq2::vst(dds, blind = FALSE),
      error = function(e) {
        if (grepl("less than 'nsub' rows", e$message, fixed = TRUE)) {
          message("[WARN] VST failed (nsub issue on a small matrix). Falling back to DESeq2::varianceStabilizingTransformation(blind=FALSE).")
          return(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
        }
        stop(e)
      }
    )
    return(list(matrix = SummarizedExperiment::assay(vsd), extras = list()))
  }

  if (method == "rlog") {
    if (ncol(counts) > rlog_max_samples && !isTRUE(force_rlog)) {
      stop("[ERROR] rlog はサンプル数が多いと非常に重くなります。--force_rlog TRUE を指定するか、サンプル数を減らしてください。")
    }
    dds <- create_transform_dds(counts, meta)
    rld <- DESeq2::rlog(dds, blind = FALSE)
    return(list(matrix = SummarizedExperiment::assay(rld), extras = list()))
  }

  if (method == "zscore") {
    return(list(matrix = row_zscore_matrix(counts), extras = list()))
  }

  stop("[ERROR] 未対応の normalization method: ", method)
}

autodetect_metadata_column <- function(meta, preferred = character(0), patterns = character(0)) {
  preferred <- preferred[preferred %in% colnames(meta)]
  if (length(preferred) > 0) {
    return(preferred[[1]])
  }

  if (length(patterns) > 0) {
    for (pat in patterns) {
      hits <- grep(pat, colnames(meta), ignore.case = TRUE, value = TRUE)
      if (length(hits) > 0) {
        return(hits[[1]])
      }
    }
  }

  NULL
}

message("[INFO] 00_load_data.R ユーティリティ読み込み完了")
