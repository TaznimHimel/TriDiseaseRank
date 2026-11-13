required_files <- c(
  "app.R",
  "R/batch_files.R",
  "R/deg_analysis.R",
  "R/wgcna_analysis.R",
  "R/cross_dataset.R",
  "R/gnn_analysis.R",
  "R/ml_analysis.R",
  "R/step8_analysis.R",
  "original/PPI_GNN_R.R",
  "original/TriDiseaseRankX_ML.R",
  "original/Pathway.R",
  "scripts/setup_windows.R",
  "scripts/check_environment.R",
  "scripts/launch_app.R"
)
missing <- required_files[!file.exists(required_files)]
if (length(missing)) stop("Missing project files: ", paste(missing, collapse = ", "))

r_files <- c(
  "app.R",
  list.files("R", pattern = "\\.R$", full.names = TRUE),
  list.files("original", pattern = "\\.R$", full.names = TRUE),
  list.files("scripts", pattern = "\\.R$", full.names = TRUE)
)
for (f in r_files) {
  parse(file = f)
}
cat("Project layout and R syntax checks passed.\n")
