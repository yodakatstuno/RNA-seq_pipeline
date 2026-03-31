options(repos = c(CRAN = Sys.getenv("CRAN_SNAPSHOT")))

install.packages("renv")
renv::consent(provided = TRUE)

libpath <- renv::paths$library()
dir.create(libpath, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(libpath, .libPaths()))

download_and_install <- function(pkg, version) {
  urls <- c(
    sprintf("%s/src/contrib/%s_%s.tar.gz", Sys.getenv("CRAN_SNAPSHOT"), pkg, version),
    sprintf("https://cran.r-project.org/src/contrib/%s_%s.tar.gz", pkg, version),
    sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", pkg, version)
  )
  dest <- file.path(tempdir(), sprintf("%s_%s.tar.gz", pkg, version))
  ok <- FALSE
  for (u in urls) {
    try(suppressWarnings(download.file(u, dest, quiet = TRUE)), silent = TRUE)
    if (file.exists(dest) && file.info(dest)$size > 0) {
      ok <- TRUE
      break
    }
  }
  if (!ok) {
    stop(sprintf("Failed to download %s %s", pkg, version))
  }
  install.packages(dest, repos = NULL, type = "source")
}

png_urls <- c(
  sprintf("%s/src/contrib/png_0.1-9.tar.gz", Sys.getenv("CRAN_SNAPSHOT")),
  "https://cran.r-project.org/src/contrib/png_0.1-9.tar.gz",
  "https://cran.r-project.org/src/contrib/Archive/png/png_0.1-9.tar.gz"
)

png_dest <- file.path(tempdir(), "png_0.1-9.tar.gz")
ok <- FALSE
for (u in png_urls) {
  try(suppressWarnings(download.file(u, png_dest, quiet = TRUE)), silent = TRUE)
  if (file.exists(png_dest) && file.info(png_dest)$size > 0) {
    ok <- TRUE
    break
  }
}
if (!ok) {
  stop("Failed to download png 0.1-9")
}
install.packages(png_dest, repos = NULL, type = "source")

download_and_install("SQUAREM", "2026.1")
download_and_install("future", "1.70.0")
download_and_install("globals", "0.19.1")
download_and_install("listenv", "0.10.1")
download_and_install("spatstat.utils", "3.2-2")
download_and_install("prodlim", "2026.03.11")

renv::restore(
  lockfile = "/app/renv.lock",
  prompt = FALSE,
  exclude = c(
    "ggtree", "enrichplot", "clusterProfiler", "ReactomePA", "png",
    "SQUAREM", "future", "globals", "listenv", "spatstat.utils", "prodlim"
  )
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(
  version = Sys.getenv("RENV_CONFIG_BIOC_VERSION"),
  ask = FALSE,
  update = FALSE
)

# Patch ggfun for ggtree check_linewidth compatibility

tmp <- tempfile("ggfun")
dir.create(tmp)
download.packages("ggfun", destdir = tmp, type = "source")

tarball <- list.files(tmp, pattern = "^ggfun_.*\\.tar\\.gz$", full.names = TRUE)[1]
untar(tarball, exdir = tmp)
src <- list.files(tmp, pattern = "^ggfun$", full.names = TRUE)[1]

writeLines(
  "check_linewidth <- function(x = NULL, ...) { x }",
  file.path(src, "R", "check_linewidth.R")
)

ns_path <- file.path(src, "NAMESPACE")
ns <- readLines(ns_path)
if (!any(grepl("^export\\(check_linewidth\\)$", ns))) {
  writeLines(c(ns, "export(check_linewidth)"), ns_path)
}
install.packages(src, repos = NULL, type = "source")

install.packages(c("ggplotify", "aplot", "yulab.utils"), repos = Sys.getenv("CRAN_SNAPSHOT"))

BiocManager::install(
  c(
    "tidytree", "treeio", "ggtree",
    "enrichplot", "clusterProfiler", "ReactomePA",
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db"
  ),
  ask = FALSE,
  update = FALSE
)

pkgs <- c(
  "ggfun", "tidytree", "treeio", "ggtree",
  "ggplotify", "aplot", "yulab.utils",
  "clusterProfiler", "enrichplot", "ReactomePA",
  "DESeq2", "edgeR", "limma", "sva", "WGCNA", "SQUAREM",
  "future", "globals", "listenv", "spatstat.utils", "prodlim",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(paste("Missing", p))
  }
  cat("[OK]", p, as.character(packageVersion(p)), "\n")
}
cat("[OK] key packages loaded\n")
