# Check whether the required runtime and packages are available.
required <- c(
  "shiny", "bslib", "DT", "ggplot2", "statmod", "zip", "WGCNA",
  "torch", "data.table", "caret", "pROC", "glmnet", "ranger",
  "mboost", "MASS", "plsRglm", "gbm", "xgboost", "naivebayes",
  "kernlab", "matrixStats", "dplyr", "tidyr", "future",
  "future.apply", "openxlsx", "edgeR", "limma",
  "clusterProfiler", "org.Hs.eg.db"
)

cat("R version: ", R.version.string, "\n", sep = "")
cat("Platform: ", R.version$platform, "\n", sep = "")
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = names(status), installed = unname(status)), row.names = FALSE)
if (any(!status)) {
  stop("Missing required packages: ", paste(names(status)[!status], collapse = ", "))
}
if (!requireNamespace("torch", quietly = TRUE)) stop("torch is missing")
tryCatch({
  torch::torch_tensor(1)
  cat("Torch CPU test: OK\n")
}, error = function(e) {
  stop("Torch runtime is not ready. Run torch::install_torch(). Details: ", conditionMessage(e))
})
cat("Environment check passed.\n")
