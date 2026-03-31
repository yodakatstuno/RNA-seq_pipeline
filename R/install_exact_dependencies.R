#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 1)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
manifest_path <- file.path(repo_root, "dependencies", "r-package-manifest.csv")
runtime_path <- file.path(repo_root, "dependencies", "runtime-versions.txt")
user_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(repo_root, ".r_libs"))

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

if (!file.exists(manifest_path)) {
  stop("Manifest not found: ", manifest_path)
}
if (!file.exists(runtime_path)) {
  stop("Runtime version file not found: ", runtime_path)
}

runtime_lines <- readLines(runtime_path, warn = FALSE)
runtime_map <- stats::setNames(
  vapply(strsplit(runtime_lines, "=", fixed = TRUE), `[`, character(1), 2),
  vapply(strsplit(runtime_lines, "=", fixed = TRUE), `[`, character(1), 1)
)

expected_r <- runtime_map[["R"]]
expected_bioc <- runtime_map[["Bioconductor"]]
if (is.na(expected_r) || is.na(expected_bioc)) {
  stop("Failed to parse ", runtime_path)
}

current_r <- paste(R.version$major, R.version$minor, sep = ".")
if (!identical(current_r, expected_r)) {
  stop("R version mismatch: expected ", expected_r, " but found ", current_r)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
manifest <- manifest[order(manifest$order), ]
compatibility_overrides <- c(
  ggplot2 = "3.5.1",
  ggtree = "4.1.1.006"
)
github_overrides <- c(
  ggtree = "1a3d0285536625825ccd597d42ca66a6ce8c17c5"
)
get_effective_version <- function(pkg, version) {
  if (pkg %in% names(compatibility_overrides)) {
    return(unname(compatibility_overrides[[pkg]]))
  }
  version
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = user_lib)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}

BiocManager::install(version = expected_bioc, ask = FALSE, update = FALSE)

pkg_version_or_na <- function(pkg) {
  ip <- installed.packages()
  if (pkg %in% rownames(ip)) {
    return(unname(ip[pkg, "Version"]))
  }
  NA_character_
}

download_cran_source <- function(pkg, version) {
  dest <- file.path(tempdir(), sprintf("%s_%s.tar.gz", pkg, version))
  urls <- c(
    sprintf("%s/src/contrib/%s_%s.tar.gz", getOption("repos")[["CRAN"]], pkg, version),
    sprintf("%s/src/contrib/Archive/%s/%s_%s.tar.gz", getOption("repos")[["CRAN"]], pkg, pkg, version)
  )
  for (url in urls) {
    try(suppressWarnings(download.file(url, dest, quiet = TRUE)), silent = TRUE)
    if (file.exists(dest) && file.info(dest)$size > 0) {
      return(dest)
    }
  }
  stop("Failed to download source tarball for ", pkg, " ", version)
}

download_bioc_source <- function(pkg, version) {
  bioc_soft <- unname(BiocManager::repositories()[["BioCsoft"]])
  dest <- file.path(tempdir(), sprintf("%s_%s.tar.gz", pkg, version))
  url <- sprintf("%s/src/contrib/%s_%s.tar.gz", bioc_soft, pkg, version)
  try(suppressWarnings(download.file(url, dest, quiet = TRUE)), silent = TRUE)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(dest)
  }
  stop("Failed to download Bioconductor source tarball for ", pkg, " ", version)
}

install_ggfun_compat <- function(version) {
  message("[INFO] Installing patched ggfun for ggtree compatibility: ", version)
  preinstall_pkgs <- c(
    "ggplotify" = "0.1.3",
    "aplot" = "0.2.9",
    "yulab.utils" = "0.2.4"
  )
  for (dep_pkg in names(preinstall_pkgs)) {
    dep_current <- pkg_version_or_na(dep_pkg)
    dep_version <- preinstall_pkgs[[dep_pkg]]
    if (!is.na(dep_current) && identical(dep_current, dep_version)) {
      message("[OK] ", dep_pkg, " ", dep_version, " already installed")
      next
    }
    message("[INFO] Preinstalling ggfun dependency: ", dep_pkg, " ", dep_version)
    remotes::install_version(
      package = dep_pkg,
      version = dep_version,
      repos = getOption("repos")[["CRAN"]],
      lib = user_lib,
      upgrade = "never",
      dependencies = c("Depends", "Imports", "LinkingTo"),
      quiet = FALSE
    )
  }
  tarball <- download_cran_source("ggfun", version)
  build_dir <- file.path(tempdir(), sprintf("ggfun-%s-src", version))
  unlink(build_dir, recursive = TRUE, force = TRUE)
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(tarball, exdir = build_dir)
  src_dir <- list.dirs(build_dir, recursive = FALSE, full.names = TRUE)[1]
  writeLines(
    "check_linewidth <- function(x = NULL, ...) { x }",
    file.path(src_dir, "R", "check_linewidth.R")
  )
  ns_path <- file.path(src_dir, "NAMESPACE")
  ns <- readLines(ns_path, warn = FALSE)
  if (!any(grepl("^export\\(check_linewidth\\)$", ns))) {
    writeLines(c(ns, "export(check_linewidth)"), ns_path)
  }
  if ("ggfun" %in% loadedNamespaces()) {
    unloadNamespace("ggfun")
  }
  install.packages(src_dir, repos = NULL, type = "source", lib = user_lib, dependencies = FALSE)
  if ("ggfun" %in% loadedNamespaces()) {
    unloadNamespace("ggfun")
  }
  ns_ggfun <- loadNamespace("ggfun")
  if (!exists("check_linewidth", envir = ns_ggfun, inherits = FALSE)) {
    stop("ggfun compatibility patch failed: check_linewidth object was not present in namespace")
  }
  if (!"check_linewidth" %in% getNamespaceExports(ns_ggfun)) {
    stop("ggfun compatibility patch failed: check_linewidth was not exported")
  }
  message("[OK] ggfun compatibility patch verified")
}

install_ggtree_compat <- function(version) {
  ref <- github_overrides[["ggtree"]]
  ggiraph_current <- pkg_version_or_na("ggiraph")
  if (is.na(ggiraph_current)) {
    message("[INFO] Preinstalling ggtree GitHub dependency: ggiraph")
    remotes::install_cran(
      pkgs = "ggiraph",
      lib = user_lib,
      repos = getOption("repos")[["CRAN"]],
      upgrade = "never",
      dependencies = c("Depends", "Imports", "LinkingTo"),
      quiet = FALSE
    )
    ggiraph_current <- pkg_version_or_na("ggiraph")
    if (is.na(ggiraph_current)) {
      stop("ggtree compatibility install failed: ggiraph was not installed")
    }
    message("[OK] ggiraph ", ggiraph_current)
  }
  message("[INFO] Installing ggtree from pinned GitHub ref: ", ref)
  remotes::install_github(
    repo = "YuLab-SMU/ggtree",
    ref = ref,
    lib = user_lib,
    upgrade = "never",
    dependencies = FALSE,
    quiet = FALSE
  )
}

install_cran_exact <- function(pkg, version) {
  effective_version <- get_effective_version(pkg, version)
  if (!identical(effective_version, version)) {
    message("[INFO] Installing CRAN package with compatibility override: ",
            pkg, " ", effective_version, " (manifest recorded ", version, ")")
  } else {
    message("[INFO] Installing CRAN package: ", pkg, " ", version)
  }
  if (identical(pkg, "ggfun")) {
    install_ggfun_compat(version)
    return(invisible(NULL))
  }
  remotes::install_version(
    package = pkg,
    version = effective_version,
    repos = getOption("repos")[["CRAN"]],
    lib = user_lib,
    upgrade = "never",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    quiet = FALSE
  )
}

install_bioc_release <- function(pkg, version) {
  message("[INFO] Installing Bioconductor package: ", pkg, " ", version,
          " (Bioconductor ", expected_bioc, ")")
  if (identical(pkg, "ggtree")) {
    install_ggtree_compat(version)
    return(invisible(NULL))
  }
  BiocManager::install(pkg, ask = FALSE, update = FALSE, force = TRUE, lib = user_lib)
}

for (i in seq_len(nrow(manifest))) {
  pkg <- manifest$package[[i]]
  version <- manifest$version[[i]]
  expected_version <- get_effective_version(pkg, version)
  repository <- manifest$repository[[i]]
  current <- pkg_version_or_na(pkg)
  expected_sha <- if (pkg %in% names(github_overrides)) github_overrides[[pkg]] else NA_character_

  current_sha <- if (!is.na(current) && !is.na(expected_sha)) {
    utils::packageDescription(pkg, fields = "RemoteSha", lib.loc = user_lib)
  } else {
    NA_character_
  }

  if (!is.na(current) && identical(current, expected_version) &&
      (is.na(expected_sha) || identical(current_sha, expected_sha))) {
    if (!identical(expected_version, version)) {
      message("[OK] ", pkg, " ", expected_version,
              " already installed (compatibility override from manifest ", version, ")")
    } else if (!is.na(expected_sha)) {
      message("[OK] ", pkg, " ", current,
              " already installed from GitHub ref ", expected_sha)
    } else {
      message("[OK] ", pkg, " ", version, " already installed")
    }
    next
  }

  if (grepl("^Bioconductor", repository)) {
    install_bioc_release(pkg, version)
  } else {
    install_cran_exact(pkg, version)
  }

  installed_version <- pkg_version_or_na(pkg)
  installed_sha <- if (!is.na(expected_sha) && !is.na(installed_version)) {
    utils::packageDescription(pkg, fields = "RemoteSha", lib.loc = user_lib)
  } else {
    NA_character_
  }
  if (is.na(installed_version) || !identical(installed_version, expected_version)) {
    stop("Version check failed for ", pkg, ": expected ", expected_version,
         " but found ", installed_version %||% "NOT_INSTALLED")
  }
  if (!is.na(expected_sha) && !identical(installed_sha, expected_sha)) {
    stop("Remote SHA check failed for ", pkg, ": expected ", expected_sha,
         " but found ", installed_sha %||% "NOT_INSTALLED")
  }
  message("[OK] ", pkg, " ", installed_version)
}

message("[INFO] Exact dependency installation complete.")
