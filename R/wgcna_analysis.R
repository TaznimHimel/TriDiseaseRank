# TriDiseaseRank-X Step 4: WGCNA module for voom-normalized DEG expression.
# Retains the original project's core choices:
# signed network, Pearson correlation, variance filtering, soft threshold,
# blockwise modules, module-trait association and |MM| * |GS| gene score.
# Does NOT run the downstream shared-candidate, GNN, ML or pathway stages.

tdrx_minmax01 <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  finite <- is.finite(x)
  if (!any(finite)) return(rep(0, length(x)))
  rr <- range(x[finite])
  if (diff(rr) == 0) return(rep(0, length(x)))
  out <- (x - rr[1]) / diff(rr)
  out[!is.finite(out)] <- 0
  pmax(0, pmin(1, out))
}

tdrx_wgcna_analysis <- function(deg_result, top_var_genes = 2000L,
                                min_module_size = 30L, module_p_cutoff = 0.05) {
  for (pkg in c("WGCNA", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required WGCNA package is missing: ", pkg)
    }
  }
  if (!is.list(deg_result) || is.null(deg_result$expression) ||
      is.null(deg_result$metadata) || is.null(deg_result$annotation)) {
    stop("Run DEG first: WGCNA requires its voom expression, metadata and annotation.")
  }
  top_var_genes <- as.integer(top_var_genes)
  min_module_size <- as.integer(min_module_size)
  if (!is.finite(top_var_genes) || top_var_genes < 100L || top_var_genes > 5000L) {
    stop("Top variable genes must be between 100 and 5000.")
  }
  if (!is.finite(min_module_size) || min_module_size < 20L) {
    stop("Minimum module size must be at least 20.")
  }
  if (!is.finite(module_p_cutoff) || module_p_cutoff <= 0 || module_p_cutoff >= 1) {
    stop("Module P-value cutoff must be strictly between 0 and 1.")
  }

  expr <- deg_result$expression
  md <- deg_result$metadata
  if (!is.matrix(expr) || !is.numeric(expr) || !all(is.finite(expr))) {
    stop("The saved voom expression matrix is missing or contains invalid values.")
  }
  if (ncol(expr) < 15L) {
    stop("WGCNA needs at least 15 samples here; fewer are not reliable for this web workflow.")
  }
  if (!all(c("Samples", "Condition") %in% names(md))) {
    stop("DEG metadata is missing Samples or Condition.")
  }
  if (!all(colnames(expr) %in% md$Samples)) {
    stop("Expression sample names do not match the DEG metadata.")
  }
  variances <- apply(expr, 1L, stats::var, na.rm = TRUE)
  eligible <- which(is.finite(variances) & variances > 0)
  if (length(eligible) < min_module_size * 2L) {
    stop("Too few variable genes available to infer modules.")
  }
  eligible <- eligible[order(variances[eligible], decreasing = TRUE)]
  eligible <- head(eligible, top_var_genes)
  datExpr <- t(expr[eligible, , drop = FALSE])
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  if (nrow(datExpr) < 15L || ncol(datExpr) < min_module_size * 2L) {
    stop("WGCNA quality filtering removed too many samples or genes.")
  }
  md <- md[match(rownames(datExpr), as.character(md$Samples)), , drop = FALSE]
  if (anyNA(md$Samples)) stop("Sample mapping failed after WGCNA quality checks.")
  trait <- as.numeric(md$Condition == "Case")
  if (length(unique(trait)) != 2L) stop("WGCNA needs both Case and Control samples.")

  powers <- c(1:10, seq(12, 30, 2))
  sft <- WGCNA::pickSoftThreshold(
    datExpr, powerVector = powers, networkType = "signed",
    corFnc = "cor", corOptions = list(use = "p"), verbose = 0
  )
  fit <- sft$fitIndices
  if (is.null(fit) || !nrow(fit) || ncol(fit) < 3L) {
    stop("WGCNA could not estimate the soft-threshold scale-free topology fit.")
  }
  r2 <- -sign(fit[, 3]) * fit[, 2]
  eligible_power <- which(is.finite(r2) & r2 >= 0.80)
  if (length(eligible_power)) {
    power_idx <- eligible_power[1L]
  } else if (any(is.finite(r2))) {
    power_idx <- which.max(replace(r2, !is.finite(r2), -Inf))
  } else {
    stop("No finite scale-free topology fit; inspect the dataset and sample variation.")
  }
  power <- as.numeric(fit[power_idx, 1L])
  net <- WGCNA::blockwiseModules(
    datExpr, power = power, maxBlockSize = max(5000L, ncol(datExpr)),
    networkType = "signed", TOMType = "signed", corType = "pearson",
    minModuleSize = min_module_size, deepSplit = 2L,
    mergeCutHeight = 0.25, reassignThreshold = 0,
    numericLabels = FALSE, pamRespectsDendro = FALSE,
    saveTOMs = FALSE, verbose = 0
  )
  colors <- as.character(net$colors)
  names(colors) <- colnames(datExpr)
  if (is.null(net$MEs) || ncol(as.matrix(net$MEs)) == 0L) {
    stop("WGCNA did not detect module eigengenes in this dataset.")
  }
  MEs <- WGCNA::orderMEs(net$MEs)
  module_names <- sub("^ME", "", colnames(MEs))
  mt_cor <- as.numeric(WGCNA::cor(MEs, trait, use = "p"))
  mt_p <- as.numeric(WGCNA::corPvalueStudent(matrix(mt_cor, ncol = 1), nrow(datExpr)))
  module_sizes <- table(colors)
  module_table <- data.frame(
    Module = module_names,
    Genes = as.integer(module_sizes[match(module_names, names(module_sizes))]),
    Correlation = mt_cor,
    P_value = mt_p,
    FDR = stats::p.adjust(mt_p, method = "BH"),
    stringsAsFactors = FALSE
  )
  module_table$Genes[is.na(module_table$Genes)] <- 0L
  module_table <- module_table[module_table$Module != "grey", , drop = FALSE]
  module_table$Significant <- !is.na(module_table$P_value) &
    module_table$P_value < module_p_cutoff
  module_table <- module_table[order(module_table$P_value, na.last = TRUE), , drop = FALSE]
  sig_modules <- module_table$Module[module_table$Significant]

  gs <- as.numeric(WGCNA::cor(datExpr, trait, use = "p"))
  names(gs) <- colnames(datExpr)
  mm <- as.matrix(WGCNA::cor(datExpr, MEs, use = "p"))
  rownames(mm) <- colnames(datExpr)
  own_mm <- vapply(colnames(datExpr), function(gene) {
    col <- paste0("ME", colors[[gene]])
    if (col %in% colnames(mm)) as.numeric(mm[gene, col]) else NA_real_
  }, numeric(1))
  genes <- colnames(datExpr)
  annotations <- deg_result$annotation
  symbols <- annotations$gene_symbol[match(genes, annotations$gene_id)]
  symbols[is.na(symbols) | !nzchar(trimws(symbols))] <-
    genes[is.na(symbols) | !nzchar(trimws(symbols))]
  gene_table <- data.frame(
    gene_id = genes,
    gene_symbol = symbols,
    Module = unname(colors[genes]),
    MM = own_mm,
    GS = unname(gs[genes]),
    stringsAsFactors = FALSE
  )
  gene_table$Significant_Module <- gene_table$Module %in% sig_modules
  gene_table$WGCNA_raw <- abs(gene_table$MM) * abs(gene_table$GS)
  gene_table$S_WGCNA <- tdrx_minmax01(gene_table$WGCNA_raw)
  gene_table <- gene_table[order(-gene_table$S_WGCNA), , drop = FALSE]
  gene_table <- gene_table[!duplicated(gene_table$gene_symbol), , drop = FALSE]
  significant_genes <- gene_table[
    gene_table$Significant_Module & gene_table$Module != "grey", , drop = FALSE
  ]

  power_fit <- data.frame(Power = as.numeric(fit[, 1]), Signed_R2 = as.numeric(r2))
  power_plot <- ggplot2::ggplot(power_fit, ggplot2::aes(x = Power, y = Signed_R2)) +
    ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed", color = "#94A3B8") +
    ggplot2::geom_line(color = "#0F766E", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.4, color = "#0F766E", na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = power, linetype = "dotted", color = "#C2410C") +
    ggplot2::labs(title = paste0("Soft-threshold selection (power = ", power, ")"),
                  x = "Candidate power", y = "Signed scale-free topology fit (R2)") +
    ggplot2::theme_minimal(base_size = 12)
  assoc_plot <- if (!nrow(module_table)) {
    NULL
  } else {
    ggplot2::ggplot(module_table,
                    ggplot2::aes(x = stats::reorder(Module, Correlation), y = Correlation,
                                 fill = Correlation > 0)) +
      ggplot2::geom_col(width = 0.68) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#0F766E", "FALSE" = "#EA580C"),
                                 guide = "none") +
      ggplot2::labs(title = "Module associations with Case vs Control",
                    x = "Module", y = "Pearson correlation") +
      ggplot2::theme_minimal(base_size = 12)
  }
  list(modules = module_table, all_genes = gene_table,
       significant_genes = significant_genes,
       soft_power = power, analyzed_genes = ncol(datExpr),
       analyzed_samples = nrow(datExpr),
       power_plot = power_plot, trait_plot = assoc_plot,
       p_cutoff = module_p_cutoff,
       significant_modules = length(sig_modules))
}
