# Batch file pairing and CSV loading.
# Input names: DATASET_counts_data.csv and DATASET_metadata.csv
# Alternatives: DATASET_counts.csv / DATASET_raw_counts.csv / DATASET_meta.csv
# No names are used as filesystem output paths without sanitization.

tdrx_dataset_stem <- function(filename, type) {
  name <- basename(as.character(filename))
  if (!grepl("\\.csv$", name, ignore.case = TRUE)) return(NA_character_)
  base <- sub("\\.csv$", "", name, ignore.case = TRUE)
  suffix <- if (identical(type, "counts")) {
    "(_counts_data|_raw_counts|_expression_counts|_counts|_count_matrix)$"
  } else {
    "(_metadata|_meta)$"
  }
  short <- sub(suffix, "", base, ignore.case = TRUE)
  if (identical(short, base) || !nzchar(trimws(short))) return(NA_character_)
  short
}

tdrx_safe_name <- function(x) {
  out <- gsub("[^A-Za-z0-9_-]", "_", as.character(x))
  ifelse(nzchar(out), out, "dataset")
}

tdrx_upload_manifest <- function(counts_upload, metadata_upload) {
  empty <- data.frame(
    key = character(), Dataset = character(), Counts = character(),
    Metadata = character(), counts_path = character(), metadata_path = character(),
    Status = character(), stringsAsFactors = FALSE
  )
  as_entries <- function(u, type) {
    if (is.null(u) || !NROW(u)) {
      return(data.frame(key = character(), Dataset = character(), name = character(),
                        path = character(), stringsAsFactors = FALSE))
    }
    stems <- vapply(u$name, tdrx_dataset_stem, character(1), type = type)
    keys <- tolower(trimws(stems))
    # Every unrecognized upload gets its own visible diagnostic row.
    bad <- is.na(keys) | !nzchar(keys)
    keys[bad] <- paste0("__bad_", type, "_", which(bad))
    stems[bad] <- paste0("Unrecognized: ", basename(u$name[bad]))
    data.frame(key = keys, Dataset = stems,
               name = basename(as.character(u$name)),
               path = as.character(u$datapath), stringsAsFactors = FALSE)
  }
  cts <- as_entries(counts_upload, "counts")
  md <- as_entries(metadata_upload, "metadata")
  keys <- unique(c(cts$key, md$key))
  if (!length(keys)) return(empty)
  rows <- lapply(keys, function(k) {
    x <- cts[cts$key == k, , drop = FALSE]
    y <- md[md$key == k, , drop = FALSE]
    label <- if (nrow(x)) x$Dataset[1] else y$Dataset[1]
    status <- if (grepl("^__bad_", k)) {
      "Rename file using _counts_data.csv or _metadata.csv suffix"
    } else if (nrow(x) > 1L || nrow(y) > 1L) {
      "Duplicate dataset names - rename/remove extra CSVs"
    } else if (!nrow(x)) {
      "Missing counts CSV"
    } else if (!nrow(y)) {
      "Missing metadata CSV"
    } else {
      "Ready"
    }
    data.frame(key = k, Dataset = label,
               Counts = if (nrow(x)) paste(x$name, collapse = ", ") else "-",
               Metadata = if (nrow(y)) paste(y$name, collapse = ", ") else "-",
               counts_path = if (nrow(x) == 1L) x$path[1] else "",
               metadata_path = if (nrow(y) == 1L) y$path[1] else "",
               Status = status, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

tdrx_read_pair <- function(row) {
  if (nrow(row) != 1L || row$Status != "Ready") {
    stop("Selected dataset does not have exactly one matching counts/metadata pair.")
  }
  if (!file.exists(row$counts_path) || !file.exists(row$metadata_path)) {
    stop("Uploaded CSV could not be found; please upload it again.")
  }
  counts <- utils::read.csv(row$counts_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
  metadata <- utils::read.csv(row$metadata_path, check.names = FALSE,
                               stringsAsFactors = FALSE)
  list(counts = counts, metadata = metadata)
}
