# Launch from repository root so all relative paths resolve consistently.
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "app.R"))) {
  stop("app.R not found. Run this script from the TriDiseaseRank-X repository folder.")
}
shiny::runApp(appDir = root, launch.browser = TRUE)
