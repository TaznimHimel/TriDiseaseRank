# TriDiseaseRank-X Windows setup
# Run from the repository root: Rscript scripts/setup_windows.R

options(timeout = 1200)
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (.Platform$OS.type != "windows") {
  stop("This setup script is intended for Windows.")
}

cran_required <- c(
  "shiny", "bslib", "DT", "ggplot2", "statmod", "zip",
  "WGCNA", "torch", "data.table", "caret", "pROC", "glmnet",
  "ranger", "mboost", "MASS", "plsRglm", "gbm", "xgboost",
  "naivebayes", "kernlab", "matrixStats", "dplyr", "tidyr",
  "future", "future.apply", "openxlsx"
)

install_cran_one <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message("[OK] ", pkg)
    return(TRUE)
  }
  message("[INSTALL] ", pkg)
  tryCatch({
    utils::install.packages(
      pkg,
      type = "binary",
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = 1
    )
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    message("[FAILED] ", pkg, ": ", conditionMessage(e))
    FALSE
  })
}

cran_status <- vapply(cran_required, install_cran_one, logical(1))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  utils::install.packages("BiocManager", type = "binary")
}

bioc_required <- c("edgeR", "limma", "clusterProfiler", "org.Hs.eg.db")

for (pkg in bioc_required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("[BIOC INSTALL] ", pkg)
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE, Ncpus = 1),
      error = function(e) message("[FAILED] ", pkg, ": ", conditionMessage(e))
    )
  } else {
    message("[OK] ", pkg)
  }
}

# Optional pathway databases: the app can run without these.
bioc_optional <- c("ReactomePA", "DOSE", "enrichplot")
for (pkg in bioc_optional) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("[OPTIONAL] Attempting ", pkg)
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE, Ncpus = 1),
      error = function(e) message("[OPTIONAL SKIPPED] ", pkg, ": ", conditionMessage(e))
    )
  }
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  tryCatch(
    utils::install.packages("msigdbr", type = "binary"),
    error = function(e) message("[OPTIONAL SKIPPED] msigdbr: ", conditionMessage(e))
  )
}

# Install the LibTorch runtime used by the GNN module.
if (requireNamespace("torch", quietly = TRUE)) {
  message("[TORCH] Checking LibTorch runtime")
  tryCatch(
    torch::install_torch(),
    error = function(e) message("[TORCH SETUP NEEDED] Run torch::install_torch() manually: ", conditionMessage(e))
  )
}

all_required <- c(cran_required, bioc_required, "BiocManager")
status <- data.frame(
  package = all_required,
  installed = vapply(all_required, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
utils::write.csv(status, "setup_package_report.csv", row.names = FALSE)
print(status, row.names = FALSE)

missing <- status$package[!status$installed]
if (length(missing)) {
  stop(
    "Required packages still missing: ",
    paste(missing, collapse = ", "),
    ". See setup_package_report.csv. Some source packages may require Rtools45."
  )
}

message("\nTriDiseaseRank-X required packages are installed.")
message("Optional pathway packages may be absent without blocking GO/KEGG analysis.")
message("Next: double-click run_app.bat")
