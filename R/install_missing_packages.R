#!/usr/bin/env Rscript
# ============================================================
# install_missing_packages.R
# 必要パッケージのチェックとインストール
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 1)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

message <- function(..., appendLF = TRUE, domain = NULL) {
  cat(..., if (appendLF) "\n" else "", sep = "")
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
default_user_lib <- file.path(repo_root, ".r_libs")
user_lib <- default_user_lib
exact_installer <- file.path(repo_root, "R", "install_exact_dependencies.R")
manifest_path <- file.path(repo_root, "dependencies", "r-package-manifest.csv")
rscript_bin <- Sys.which("Rscript")
if (!nzchar(rscript_bin)) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
}

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_LIBS_USER = user_lib)
.libPaths(unique(c(user_lib, .libPaths())))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  KMP_DUPLICATE_LIB_OK = "TRUE",
  KMP_INIT_AT_FORK = "FALSE",
  KMP_SHM_DISABLE = "1"
)

fallback_required <- c(
  "optparse", "readxl", "ggplot2", "pheatmap", "ggdendro", "ggrepel",
  "Rtsne", "uwot", "reshape2", "RColorBrewer", "igraph", "caret",
  "randomForest", "kernlab", "Seurat", "dplyr", "rlang", "getopt",
  "DESeq2", "edgeR", "limma", "clusterProfiler", "enrichplot",
  "ReactomePA", "sva", "impute", "preprocessCore",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db", "WGCNA"
)

required_packages <- if (file.exists(manifest_path)) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  unique(manifest$package)
} else {
  fallback_required
}

openmp_sensitive <- c("caret", "Seurat", "clusterProfiler", "enrichplot", "ReactomePA", "WGCNA")
quarantine_root <- file.path(user_lib, "_quarantine")

message("[INFO] Repository root: ", repo_root)
message("[INFO] Active user library: ", user_lib)
message("[INFO] Library paths: ", paste(.libPaths(), collapse = " | "))

is_native_load_error <- function(text) {
  grepl(
    "unable to load shared object|invalid mach-o file|slice is not valid mach-o file|wrong architecture|image not found|undefined symbol",
    text,
    ignore.case = TRUE
  )
}

probe_package <- function(pkg, lib_paths = .libPaths()) {
  lib_expr <- paste(sprintf("'%s'", lib_paths), collapse = ", ")
  expr <- paste0(
    "libpaths <- c(", lib_expr, ");",
    ".libPaths(unique(libpaths));",
    "suppressPackageStartupMessages(library('", pkg, "', character.only = TRUE));",
    "cat('OK\\n')"
  )
  script_file <- tempfile("probe-", fileext = ".R")
  stdout_file <- tempfile("probe-", fileext = ".out")
  stderr_file <- tempfile("probe-", fileext = ".err")
  writeLines(expr, script_file)
  on.exit(unlink(c(script_file, stdout_file, stderr_file), force = TRUE), add = TRUE)
  status <- suppressWarnings(system2(
    rscript_bin,
    c("--vanilla", script_file),
    stdout = stdout_file,
    stderr = stderr_file
  )) %||% 0L
  out <- c(readLines(stdout_file, warn = FALSE), readLines(stderr_file, warn = FALSE))
  text <- paste(out, collapse = "\n")
  if (status == 0L) {
    return(list(status = "ok", output = text))
  }
  if (grepl("OMP: Error #179", text, fixed = TRUE) ||
      grepl("Abort trap", text, fixed = TRUE) ||
      identical(status, 134L)) {
    return(list(status = "openmp_error", output = text))
  }
  if (is_native_load_error(text)) {
    return(list(status = "native_load_error", output = text))
  }
  list(status = "load_error", output = text)
}

quarantine_package <- function(pkg, reason) {
  pkg_dir <- file.path(user_lib, pkg)
  if (!dir.exists(pkg_dir)) {
    return(invisible(FALSE))
  }
  dir.create(quarantine_root, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  target <- file.path(quarantine_root, paste0(pkg, "-", stamp))
  ok <- file.rename(pkg_dir, target)
  if (ok) {
    message("[WARN] Quarantined incompatible package from user library: ", pkg)
    message("[WARN]   Reason: ", reason)
    message("[WARN]   Moved to: ", target)
    return(invisible(TRUE))
  }
  message("[WARN] Failed to quarantine package: ", pkg)
  invisible(FALSE)
}

user_library_has_platform_mismatch <- function() {
  native_files <- list.files(
    user_lib,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(so|dylib|dll)$"
  )
  native_files <- native_files[!grepl("/_quarantine/", native_files, fixed = FALSE)]
  if (length(native_files) == 0) {
    return(FALSE)
  }

  file_out <- suppressWarnings(system2("file", native_files, stdout = TRUE, stderr = TRUE))
  any(grepl(": .*ELF ", file_out))
}

sanitize_user_library <- function() {
  user_dirs <- list.dirs(user_lib, recursive = FALSE, full.names = FALSE)
  user_dirs <- user_dirs[!grepl("^00LOCK", user_dirs)]
  user_dirs <- setdiff(user_dirs, "_quarantine")
  if (length(user_dirs) == 0) {
    return(invisible(NULL))
  }

  if (user_library_has_platform_mismatch()) {
    message("[WARN] Detected non-native compiled artifacts in ", user_lib)
    message("[WARN] Quarantining project-local library to avoid shadowing system R packages.")
    for (pkg in user_dirs) {
      quarantine_package(pkg, reason = "Project-local library contains cross-platform compiled artifacts")
    }
    .libPaths(unique(c(user_lib, .Library)))
    return(invisible(NULL))
  }

  for (pkg in user_dirs) {
    probe <- probe_package(pkg, lib_paths = unique(c(user_lib, .Library)))
    if (probe$status != "ok" &&
        grepl(normalizePath(user_lib, winslash = "/", mustWork = TRUE), probe$output, fixed = TRUE)) {
      quarantine_package(pkg, reason = strsplit(probe$output, "\n", fixed = TRUE)[[1]][1])
    }
  }

  .libPaths(unique(c(user_lib, .libPaths())))
}

install_exact_dependencies <- function() {
  if (!file.exists(exact_installer)) {
    message("[WARN] Exact dependency installer not found: ", exact_installer)
    return(invisible(FALSE))
  }

  message("[INFO] Installing missing packages via install_exact_dependencies.R")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)
  out <- suppressWarnings(system2(
    rscript_bin,
    c("--vanilla", exact_installer),
    stdout = TRUE,
    stderr = TRUE,
    env = c(sprintf("R_LIBS_USER=%s", user_lib))
  ))
  status <- attr(out, "status") %||% 0L
  if (length(out) > 0) {
    for (line in out) {
      message(line)
    }
  }
  result <- identical(status, 0L)
  if (!result) {
    summary_line <- if (length(out) > 0) tail(out, 1) else "unknown error"
    message("[WARN] Exact dependency installation did not complete: ", summary_line)
  }

  .libPaths(unique(c(user_lib, .libPaths())))
  invisible(result)
}

classify_required_packages <- function() {
  status_map <- setNames(rep("unknown", length(required_packages)), required_packages)
  details <- setNames(as.list(rep("", length(required_packages))), required_packages)

  for (pkg in required_packages) {
    probe <- probe_package(pkg)
    status_map[[pkg]] <- probe$status
    details[[pkg]] <- probe$output
  }

  list(
    ok = names(status_map)[status_map == "ok"],
    openmp = names(status_map)[status_map == "openmp_error"],
    missing = names(status_map)[status_map %in% c("native_load_error", "load_error")],
    status = status_map,
    details = details
  )
}

sanitize_user_library()
initial <- classify_required_packages()

if (length(initial$missing) > 0) {
  install_exact_dependencies()
}

final <- classify_required_packages()
openmp_failures <- intersect(final$openmp, openmp_sensitive)
hard_failures <- setdiff(final$missing, openmp_sensitive)

if (length(hard_failures) > 0) {
  message("[ERROR] The following packages could not be loaded: ",
          paste(hard_failures, collapse = ", "))
  for (pkg in hard_failures) {
    pkg_output <- final$details[[pkg]]
    if (nzchar(pkg_output)) {
      message("[ERROR] ", pkg, ": ", strsplit(pkg_output, "\n", fixed = TRUE)[[1]][1])
    }
  }
  quit(status = 1)
}

if (length(openmp_failures) > 0) {
  message("[WARN] OpenMP runtime is unavailable for: ",
          paste(openmp_failures, collapse = ", "))
  message("[WARN] The main pipeline will handle this by skipping OpenMP-sensitive steps when needed.")
}

message("[INFO] Package check complete.")
