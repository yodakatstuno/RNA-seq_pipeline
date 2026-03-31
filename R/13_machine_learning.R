#!/usr/bin/env Rscript
# ============================================================
# 13_machine_learning.R
# Random Forest / SVM による分類・予測
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
  library(caret)
})

option_list <- c(base_option_list, list(
  make_option("--target_col", type = "character", default = "condition"),
  make_option("--method",     type = "character", default = "both"),
  make_option("--top_n",      type = "integer",   default = 500),
  make_option("--cv_folds",   type = "integer",   default = 5),
  make_option("--test_ratio", type = "double",    default = 0.3)
))

opt <- parse_args(OptionParser(option_list = option_list))
out_sub <- resolve_step_outdir("13_ml", opt$outdir)

counts_raw <- load_counts(opt$counts)
meta_raw   <- load_metadata(opt$metadata)
dat        <- align_data(counts_raw, meta_raw)
counts     <- dat$counts
meta       <- dat$metadata

log_counts <- log2(counts + 1)
filtered <- select_top_variable(log_counts, opt$top_n)

# 特徴量行列 (samples x genes)
X <- t(filtered)
y <- factor(meta[[opt$target_col]])

message("[INFO] サンプル数: ", nrow(X), ", 特徴量数: ", ncol(X))
message("[INFO] クラス: ", paste(levels(y), collapse = ", "))

# Train/Test split
set.seed(42)
train_idx <- createDataPartition(y, p = 1 - opt$test_ratio, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

ctrl <- trainControl(method = "cv", number = opt$cv_folds, classProbs = TRUE,
                      savePredictions = "final")

methods <- if (opt$method == "both") c("rf", "svm") else opt$method
results_summary <- list()

# ---- Random Forest ----
if ("rf" %in% methods) {
  suppressPackageStartupMessages(library(randomForest))
  message("[ML] Random Forest 学習中...")

  set.seed(42)
  rf_model <- train(x = X_train, y = y_train, method = "rf",
                     trControl = ctrl, ntree = 500)

  pred_rf <- predict(rf_model, X_test)
  cm_rf <- confusionMatrix(pred_rf, y_test)

  results_summary[["RF"]] <- cm_rf$overall["Accuracy"]

  # 変数重要度
  imp <- varImp(rf_model)$importance
  imp$gene <- rownames(imp)
  imp <- imp[order(-imp[, 1]), ]
  write.csv(imp, file.path(out_sub, "rf_variable_importance.csv"), row.names = FALSE)

  top_imp <- head(imp, 30)
  top_imp$gene <- factor(top_imp$gene, levels = rev(top_imp$gene))
  p_imp <- ggplot(top_imp, aes(x = gene, y = top_imp[, 1])) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Random Forest - Top 30 Important Genes", x = "", y = "Importance")
  ggsave(file.path(out_sub, "rf_importance.pdf"), p_imp, width = 10, height = 8)
  ggsave(file.path(out_sub, "rf_importance.png"), p_imp, width = 10, height = 8, dpi = 300)

  # Confusion matrix 保存
  sink(file.path(out_sub, "rf_confusion_matrix.txt"))
  print(cm_rf)
  sink()
}

# ---- SVM ----
if ("svm" %in% methods) {
  suppressPackageStartupMessages(library(kernlab))
  message("[ML] SVM 学習中...")

  set.seed(42)
  svm_model <- train(x = X_train, y = y_train, method = "svmRadial",
                      trControl = ctrl, preProcess = c("center", "scale"))

  pred_svm <- predict(svm_model, X_test)
  cm_svm <- confusionMatrix(pred_svm, y_test)

  results_summary[["SVM"]] <- cm_svm$overall["Accuracy"]

  sink(file.path(out_sub, "svm_confusion_matrix.txt"))
  print(cm_svm)
  sink()
}

# ---- サマリー ----
sink(file.path(out_sub, "ml_summary.txt"))
cat("========================================\n")
cat("機械学習解析サマリー\n")
cat("========================================\n")
cat("Target:", opt$target_col, "\n")
cat("Features:", ncol(X), "genes\n")
cat("Samples: Train=", length(y_train), " Test=", length(y_test), "\n")
cat("CV Folds:", opt$cv_folds, "\n\n")
for (m in names(results_summary)) {
  cat(m, "Accuracy:", round(results_summary[[m]], 4), "\n")
}
sink()

message("[完了] 13_machine_learning.R")
