# TriDiseaseRank-X: reusable single-dataset DEG module
# Based on the project's Phase 1 approach: edgeR TMM -> limma-voom -> limma.
# This stage intentionally omits SVA, batch/age/sex covariates, and paired designs.

tdrx_prepare_input <- function(counts_df, metadata_df) {
  if (!is.data.frame(counts_df) || ncol(counts_df) < 4L) {
    stop("Counts CSV must contain gene_id, gene_symbol, and at least two sample columns.")
  }
  if (!is.data.frame(metadata_df)) stop("Metadata must be a CSV table.")
  names(metadata_df) <- trimws(names(metadata_df))
  if (!all(c("Samples", "Condition") %in% names(metadata_df))) {
    stop("Metadata must have columns named exactly Samples and Condition.")
  }
  md <- metadata_df
  md$Samples <- trimws(as.character(md$Samples))
  condition <- trimws(as.character(md$Condition))
  if (anyNA(md$Samples) || any(!nzchar(md$Samples))) {
    stop("Metadata contains blank sample names.")
  }
  if (anyDuplicated(md$Samples)) stop("Metadata contains duplicate sample names.")
  if (anyNA(condition) || any(!condition %in% c("Control", "Case"))) {
    stop("Condition must contain only Control and Case labels (case-sensitive).")
  }
  md$Condition <- factor(condition, levels = c("Control", "Case"))
  group_counts <- table(md$Condition)
  if (any(group_counts < 2L)) {
    stop("Please provide at least two independent samples per group for this first version.")
  }
  sample_columns <- names(counts_df)[-(1:2)]
  if (anyDuplicated(sample_columns)) stop("Counts CSV contains duplicate sample column names.")
  missing_samples <- setdiff(md$Samples, sample_columns)
  if (length(missing_samples)) {
    stop("Metadata sample IDs missing from counts: ",
         paste(head(missing_samples, 10L), collapse = ", "))
  }
  selected <- counts_df[, md$Samples, drop = FALSE]
  numeric_columns <- vapply(selected, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop("Non-numeric count columns: ",
         paste(names(selected)[!numeric_columns], collapse = ", "))
  }
  counts <- as.matrix(selected)
  storage.mode(counts) <- "double"
  if (anyNA(counts) || any(!is.finite(counts))) {
    stop("Counts contain missing, NaN, or infinite values.")
  }
  if (any(counts < 0)) stop("Counts must be non-negative.")
  if (any(abs(counts - round(counts)) > 1e-8)) {
    stop("This DEG workflow needs raw integer counts, not TPM/FPKM/normalized values.")
  }
  if (any(colSums(counts) == 0)) stop("One or more samples have zero total counts.")

  gene_ids <- trimws(as.character(counts_df[[1]]))
  gene_symbols <- toupper(trimws(as.character(counts_df[[2]])))
  valid_rows <- !is.na(gene_ids) & nzchar(gene_ids) & rowSums(counts) > 0
  counts <- counts[valid_rows, , drop = FALSE]
  gene_ids <- gene_ids[valid_rows]
  gene_symbols <- gene_symbols[valid_rows]
  if (nrow(counts) < 2L) stop("Not enough expressed genes after removing zero-count rows.")
  rownames(counts) <- gene_ids

  # Mirror the original research script's duplicate gene-ID policy:
  # retain the copy with the largest total count.
  if (anyDuplicated(gene_ids)) {
    ix <- order(gene_ids, -rowSums(counts))
    counts <- counts[ix, , drop = FALSE]
    gene_ids <- gene_ids[ix]
    gene_symbols <- gene_symbols[ix]
    ix_unique <- !duplicated(gene_ids)
    counts <- counts[ix_unique, , drop = FALSE]
    gene_ids <- gene_ids[ix_unique]
    gene_symbols <- gene_symbols[ix_unique]
  }
  invalid_symbols <- is.na(gene_symbols) | gene_symbols %in%
    c("", "NA", "N/A", "---", "--", ".", "NULL")
  gene_symbols[invalid_symbols] <- gene_ids[invalid_symbols]
  rownames(counts) <- gene_ids

  list(counts = counts, metadata = md,
       gene_ids = gene_ids, gene_symbols = gene_symbols,
       group_counts = group_counts)
}

tdrx_deg_analysis <- function(prepared, fdr_cutoff = 0.05, log2fc_cutoff = 1) {
  for (pkg in c("edgeR", "limma", "statmod", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("R package not installed: ", pkg)
    }
  }
  if (!is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff >= 1) {
    stop("FDR cutoff must be strictly between 0 and 1.")
  }
  if (!is.finite(log2fc_cutoff) || log2fc_cutoff < 0) {
    stop("Absolute log2FC cutoff must be non-negative.")
  }

  design <- stats::model.matrix(~ Condition, data = prepared$metadata)
  if (qr(design)$rank != ncol(design)) stop("The experimental design matrix is not full-rank.")
  dge <- edgeR::DGEList(counts = prepared$counts)
  keep <- edgeR::filterByExpr(dge, group = prepared$metadata$Condition, min.count = 10)
  if (sum(keep) < 2L) stop("Not enough genes remain after low-expression filtering.")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
  tab <- limma::topTable(fit, coef = "ConditionCase", number = Inf,
                         adjust.method = "BH", sort.by = "P")
  tab$gene_id <- rownames(tab)
  tab$gene_symbol <- prepared$gene_symbols[
    match(tab$gene_id, prepared$gene_ids)
  ]
  tab$padj <- tab$adj.P.Val
  tab$DEG_raw <- -log10(pmax(tab$padj, 1e-300)) * abs(tab$logFC)
  score <- tab$DEG_raw
  score[!is.finite(score)] <- NA_real_
  limits <- range(score, na.rm = TRUE)
  if (all(is.finite(limits)) && diff(limits) > 0) {
    tab$S_DEG <- (score - limits[1]) / diff(limits)
  } else {
    tab$S_DEG <- rep(0, nrow(tab))
  }
  tab$S_DEG[!is.finite(tab$S_DEG)] <- 0
  # The source project ranks one representative row per gene symbol.
  tab <- tab[!duplicated(tab$gene_symbol), , drop = FALSE]
  tab$Significant <- !is.na(tab$padj) & tab$padj < fdr_cutoff &
    abs(tab$logFC) > log2fc_cutoff
  significant <- tab[tab$Significant, , drop = FALSE]
  tab$minus_log10_FDR <- -log10(pmax(tab$padj, 1e-300))

  volcano <- ggplot2::ggplot(
    tab, ggplot2::aes(x = logFC, y = minus_log10_FDR,
                      color = Significant)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff),
                        linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("FALSE" = "#64748b", "TRUE" = "#e11d48")) +
    ggplot2::labs(title = "Differential gene expression",
                  x = "log2 fold change (Case vs Control)",
                  y = "-log10 adjusted p-value") +
    ggplot2::theme_minimal()

  # Simple exploratory PCA. Unlike the original complete script, this first
  # web version uses the 2,000 most variable genes, without scaling columns.
  expr <- v$E
  gene_variances <- apply(expr, 1L, stats::var)
  variable_genes <- which(is.finite(gene_variances) & gene_variances > 0)
  pca_plot <- NULL
  if (length(variable_genes) >= 2L && ncol(expr) >= 3L) {
    selection <- head(variable_genes[order(gene_variances[variable_genes],
                                           decreasing = TRUE)], 2000L)
    pca <- stats::prcomp(t(expr[selection, , drop = FALSE]), center = TRUE,
                        scale. = FALSE)
    variance_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    pc_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                        Condition = prepared$metadata$Condition,
                        Samples = prepared$metadata$Samples)
    pca_plot <- ggplot2::ggplot(pc_df,
                  ggplot2::aes(x = PC1, y = PC2, color = Condition)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(title = "PCA of voom-normalized gene expression",
                    x = paste0("PC1 (", variance_pct[1], "%)"),
                    y = paste0("PC2 (", variance_pct[2], "%)")) +
      ggplot2::theme_minimal()
  }
  list(all = tab, significant = significant,
       volcano = volcano, pca = pca_plot, expression = expr,
       metadata = prepared$metadata,
       annotation = data.frame(gene_id = prepared$gene_ids,
                               gene_symbol = prepared$gene_symbols,
                               stringsAsFactors = FALSE),
       fdr_cutoff = fdr_cutoff, log2fc_cutoff = log2fc_cutoff,
       tested_genes = nrow(tab), significant_genes = nrow(significant),
       analyzed_samples = ncol(expr))
}
