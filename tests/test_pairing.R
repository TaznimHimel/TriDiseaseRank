# Optional lightweight test; run from the project root: source("tests/test_pairing.R")
source(file.path("R", "batch_files.R"))
source(file.path("R", "deg_analysis.R"))
source(file.path("R", "wgcna_analysis.R"))

make_upload <- function(names) {
  data.frame(name = names, datapath = paste0("/temporary/", seq_along(names)),
             stringsAsFactors = FALSE)
}
counts <- make_upload(c("GSE120852_counts_data.csv", "GSE161349_counts_data.csv",
                        "GSE172381_counts_data.csv"))
metadata <- make_upload(c("GSE120852_metadata.csv", "GSE161349_metadata.csv",
                          "GSE172381_metadata.csv"))
pairs <- tdrx_upload_manifest(counts, metadata)
stopifnot(nrow(pairs) == 3L, all(pairs$Status == "Ready"))
stopifnot(identical(sort(pairs$key), sort(c("gse120852", "gse161349", "gse172381"))))

missing_md <- tdrx_upload_manifest(counts, metadata[1:2, , drop = FALSE])
stopifnot(sum(missing_md$Status == "Missing metadata CSV") == 1L)

bad <- tdrx_upload_manifest(make_upload("mydata.csv"), metadata[0, , drop = FALSE])
stopifnot(nrow(bad) == 1L, grepl("Rename file", bad$Status[[1]]))

also_valid <- tdrx_upload_manifest(make_upload("Trial_counts.csv"),
                                   make_upload("Trial_meta.csv"))
stopifnot(also_valid$Status[[1]] == "Ready")

prep <- tdrx_prepare_input(
  data.frame(gene_id = c("ENSG1", "ENSG2", "ENSG3"),
             gene_symbol = c("A", "B", "C"), S1 = c(4L, 5L, 6L),
             S2 = c(5L, 7L, 8L), S3 = c(7L, 4L, 8L), S4 = c(8L, 5L, 9L)),
  data.frame(Samples = c("S1", "S2", "S3", "S4"),
             Condition = c("Control", "Control", "Case", "Case"))
)
stopifnot(ncol(prep$counts) == 4L, nrow(prep$counts) == 3L)
stopifnot(identical(tdrx_minmax01(c(2, 4, 6)), c(0, .5, 1)))
cat("PASS: filename pairing, missing/invalid cases, basic input validation and score scaling.\n")
