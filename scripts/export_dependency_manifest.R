#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

suppressWarnings(options(stringsAsFactors = FALSE))

repo_root <- normalizePath(getwd(), mustWork = TRUE)
manifest_dir <- file.path(repo_root, "dependencies")
manifest_path <- file.path(manifest_dir, "r-package-manifest.csv")
runtime_path <- file.path(manifest_dir, "runtime-versions.txt")

dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

ordered_packages <- c(
  "BiocManager",
  "Rcpp", "RcppArmadillo", "BH", "cpp11", "data.table",
  "png", "XML", "curl", "openssl", "httr", "jsonlite", "stringi",
  "systemfonts", "textshaping", "farver", "isoband",
  "future", "future.apply", "globals", "listenv", "parallelly",
  "matrixStats", "SQUAREM", "spatstat.utils", "prodlim",
  "ggfun", "ggplotify", "aplot", "yulab.utils", "tidytree",
  "getopt", "optparse", "yaml", "readxl", "ggplot2", "pheatmap",
  "ggdendro", "ggrepel", "Rtsne", "uwot", "reshape2", "RColorBrewer",
  "igraph", "caret", "randomForest", "kernlab", "dplyr", "rlang",
  "Seurat",
  "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb",
  "GenomicRanges", "Biobase", "SummarizedExperiment", "AnnotationDbi",
  "preprocessCore", "impute", "limma", "edgeR", "DESeq2", "sva",
  "DOSE", "GOSemSim", "fgsea", "treeio", "ggtree", "enrichplot",
  "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "org.Mm.eg.db",
  "org.Dr.eg.db", "WGCNA"
)

ip <- installed.packages(fields = c("Repository", "RemoteType", "RemoteRepo"))
missing <- ordered_packages[!ordered_packages %in% rownames(ip)]
if (length(missing) > 0) {
  stop("The following packages are not installed locally and cannot be captured exactly: ",
       paste(missing, collapse = ", "))
}

repo_value <- function(pkg) {
  repo <- ip[pkg, "Repository"]
  if (!is.na(repo) && nzchar(repo)) {
    return(repo)
  }
  if (grepl("^org\\..*\\.eg\\.db$", pkg) || pkg %in% c(
    "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges",
    "Biobase", "SummarizedExperiment", "AnnotationDbi", "preprocessCore",
    "impute", "limma", "edgeR", "DESeq2", "sva", "DOSE", "GOSemSim",
    "fgsea", "treeio", "ggtree", "enrichplot", "clusterProfiler", "ReactomePA"
  )) {
    return(paste("Bioconductor", as.character(BiocManager::version())))
  }
  remote_type <- ip[pkg, "RemoteType"]
  remote_repo <- ip[pkg, "RemoteRepo"]
  if (!is.na(remote_type) && nzchar(remote_type)) {
    return(paste(remote_type, remote_repo))
  }
  "unknown"
}

manifest <- data.frame(
  order = seq_along(ordered_packages),
  package = ordered_packages,
  version = unname(ip[ordered_packages, "Version"]),
  repository = vapply(ordered_packages, repo_value, character(1)),
  stringsAsFactors = FALSE
)

write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)

writeLines(
  c(
    paste("R=", R.version$major, ".", R.version$minor, sep = ""),
    paste("Bioconductor=", as.character(BiocManager::version()), sep = ""),
    paste("Generated=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "")
  ),
  runtime_path
)

message("[INFO] Wrote ", manifest_path)
message("[INFO] Wrote ", runtime_path)
