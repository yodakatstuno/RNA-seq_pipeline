#!/usr/bin/env Rscript
# ============================================================
# generate_example.R
# Creates a synthetic RNA-seq dataset for testing/debugging
# Output: example/counts.rds, example/metadata.csv
# ============================================================

message("[INFO] Generating example dataset...")

set.seed(42)

n_genes <- 300
genotypes <- c("WT", "MUT")
dpf_levels <- c(3, 5, 7)
replicates <- 2

design <- expand.grid(
  genotype = genotypes,
  dpf = dpf_levels,
  rep = seq_len(replicates),
  stringsAsFactors = FALSE
)
design <- design[order(design$genotype, design$dpf, design$rep), , drop = FALSE]
n_samples <- nrow(design)

sample_ids <- paste0("sample_", sprintf("%02d", seq_len(n_samples)))
condition <- ifelse(design$genotype == "WT", "control", "treatment")
batch <- rep(c("A", "B"), length.out = n_samples)

meta <- data.frame(
  sample_id = sample_ids,
  condition = condition,
  genotype = design$genotype,
  dpf = design$dpf,
  batch = batch,
  replicate = design$rep,
  stringsAsFactors = FALSE
)

gene_names <- paste0("Gene_", sprintf("%03d", seq_len(n_genes)))
base_means <- runif(n_genes, min = 30, max = 3000)

counts <- matrix(
  0,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_names, sample_ids)
)

time_scale <- c(`3` = 0.85, `5` = 1.10, `7` = 1.40)
geno_scale <- c(WT = 1.0, MUT = 1.15)

for (j in seq_len(n_samples)) {
  sample_mu <- base_means * time_scale[as.character(meta$dpf[j])] * geno_scale[meta$genotype[j]]

  if (meta$genotype[j] == "MUT") {
    sample_mu[1:10] <- sample_mu[1:10] * (1.8 + 0.25 * meta$dpf[j])
    sample_mu[11:20] <- sample_mu[11:20] / (1.6 + 0.15 * meta$dpf[j])
  }

  sample_mu[21:30] <- sample_mu[21:30] * c(`3` = 0.7, `5` = 1.0, `7` = 1.5)[as.character(meta$dpf[j])]
  sample_mu[31:40] <- sample_mu[31:40] * c(`3` = 1.4, `5` = 1.0, `7` = 0.65)[as.character(meta$dpf[j])]

  if (meta$genotype[j] == "MUT" && meta$dpf[j] == 7) {
    sample_mu[41:55] <- sample_mu[41:55] * 2.5
  }

  counts[, j] <- rnbinom(n_genes, mu = sample_mu, size = 12)
}

counts <- round(counts)
storage.mode(counts) <- "integer"

if (!dir.exists("example")) {
  dir.create("example", recursive = TRUE)
}
saveRDS(counts, "example/counts.rds")
write.csv(meta, "example/metadata.csv", row.names = FALSE)

message("[INFO] Example dataset created:")
message("[INFO]   Counts: example/counts.rds (", n_genes, " genes x ", n_samples, " samples)")
message("[INFO]   Metadata: example/metadata.csv")
message("[INFO]   Condition levels: ", paste(unique(meta$condition), collapse = ", "))
message("[INFO]   Genotypes: ", paste(unique(meta$genotype), collapse = ", "))
message("[INFO]   DPF levels: ", paste(unique(meta$dpf), collapse = ", "))
