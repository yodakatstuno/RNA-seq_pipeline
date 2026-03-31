#!/usr/bin/env Rscript
# ============================================================
# generate_example.R
# Creates a minimal synthetic RNA-seq dataset for testing
# Output: example/counts.rds, example/metadata.csv
# ============================================================

message("[INFO] Generating example dataset...")

set.seed(42)

# Parameters
n_genes <- 200
n_samples <- 12
conditions <- rep(c("control", "treatment"), each = 6)
dpf_levels <- rep(c(3, 5), times = 6)
batch <- rep(c("A", "B"), each = 6)

# Sample IDs
sample_ids <- paste0("sample_", sprintf("%02d", 1:n_samples))

# Metadata
meta <- data.frame(
  sample_id = sample_ids,
  condition = conditions,
  dpf = dpf_levels,
  batch = batch,
  stringsAsFactors = FALSE
)

# Base expression (Poisson-like counts)
gene_names <- paste0("Gene_", sprintf("%03d", 1:n_genes))
base_means <- runif(n_genes, min = 50, max = 5000)

counts <- matrix(0, nrow = n_genes, ncol = n_samples,
                 dimnames = list(gene_names, sample_ids))

for (j in 1:n_samples) {
  counts[, j] <- rnbinom(n_genes, mu = base_means, size = 10)
}

# Inject DE genes (first 20 genes are differentially expressed)
n_de <- 20
for (j in which(conditions == "treatment")) {
  # Up-regulated in treatment (genes 1-10)
  counts[1:10, j] <- counts[1:10, j] * sample(c(3, 4, 5), 10, replace = TRUE)
  # Down-regulated in treatment (genes 11-20)
  counts[11:n_de, j] <- round(counts[11:n_de, j] / sample(c(3, 4), 10, replace = TRUE))
}

# Ensure integer counts
counts <- round(counts)
storage.mode(counts) <- "integer"

# Save
script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_path)) script_path <- "."
out_dir <- file.path(dirname(script_path), ".")
if (!dir.exists("example")) dir.create("example", recursive = TRUE)
saveRDS(counts, "example/counts.rds")
write.csv(meta, "example/metadata.csv", row.names = FALSE)

message("[INFO] Example dataset created:")
message("[INFO]   Counts: example/counts.rds (", n_genes, " genes x ", n_samples, " samples)")
message("[INFO]   Metadata: example/metadata.csv")
message("[INFO]   DE genes: ", n_de, " (genes 1-10 up, 11-20 down in treatment)")
message("[INFO]   Conditions: ", paste(unique(conditions), collapse = ", "))
message("[INFO]   DPF levels: ", paste(unique(dpf_levels), collapse = ", "))
