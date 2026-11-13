
# ============================================
# TriDiseaseRank-X STEP 6
# TD-CCHI GNN Batch Analysis Engine
# ============================================

tdrx_gnn_clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "---", "--", ".", "NULL")] <- NA_character_
  x
}

tdrx_gnn_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "_", as.character(x))
  ifelse(nzchar(x), x, "dataset")
}


tdrx_gnn_batch_analysis <- function(
    step5_result,
    deg_results,
    dataset_names = NULL,
    candidate_limit = 40L,
    mode = c("quick", "research"),
    output_root = file.path(
      "TriDiseaseRank_results", "Step6_TD_CCHI"
    ),
    progress = NULL
) {

  mode <- match.arg(mode)

  # 1. Check packages
  needed <- c("torch", "data.table", "ggplot2")

  missing <- needed[
    !vapply(
      needed,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing)) {
    stop(
      "Install required packages: ",
      paste(missing, collapse = ", ")
    )
  }

  # 2. Validate Step 5
  if (
    is.null(step5_result$ranked) ||
    !is.data.frame(step5_result$ranked) ||
    nrow(step5_result$ranked) < 3L
  ) {
    stop(
      "Run Step 5 first. At least 3 ranked candidate genes are required."
    )
  }

  keys <- names(deg_results)

  if (length(keys) < 2L) {
    stop("GNN requires DEG results from at least 2 datasets.")
  }

  candidate_limit <- as.integer(candidate_limit)

  if (
    is.na(candidate_limit) ||
    candidate_limit < 3L ||
    candidate_limit > 100L
  ) {
    stop("Candidate limit must be between 3 and 100.")
  }

  # 3. Load original GNN research code
  gnn_candidates <- c(
    "original/PPI_GNN_R.R",
    "PPI_GNN_R.R",
    "original/TD_CCHI_GNN_R.R",
    "TD_CCHI_GNN_R.R"
  )
  gnn_found <- gnn_candidates[file.exists(gnn_candidates)]

  if (!length(gnn_found)) {
    stop("Original GNN script not found. Expected original/PPI_GNN_R.R.")
  }

  gnn_file <- gnn_found[1L]

  # Private environment prevents function-name conflicts
  core_env <- new.env(parent = globalenv())

  sys.source(
    gnn_file,
    envir = core_env
  )

  run_core <- get(
    "run_td_cchi_r",
    envir = core_env
  )

  # 4. Dataset labels
  labels <- keys

  if (!is.null(dataset_names)) {
    supplied <- as.character(dataset_names[keys])
    ok <- !is.na(supplied) & nzchar(supplied)
    labels[ok] <- supplied[ok]
  }

  # 5. Find genes common to every dataset
  symbols <- setNames(
    vector("list", length(keys)),
    keys
  )

  for (key in keys) {

    d <- deg_results[[key]]

    if (
      is.null(d$expression) ||
      is.null(d$metadata) ||
      is.null(d$annotation)
    ) {
      stop("Missing DEG data for ", key)
    }

    expr <- as.matrix(d$expression)
    ann <- d$annotation

    if (
      is.null(rownames(expr)) ||
      !all(c("gene_id", "gene_symbol") %in% names(ann))
    ) {
      stop("Invalid gene annotation for ", key)
    }

    ix <- match(
      rownames(expr),
      as.character(ann$gene_id)
    )

    sym <- tdrx_gnn_clean_gene(
      ann$gene_symbol[ix]
    )

    miss <- is.na(sym)

    sym[miss] <- tdrx_gnn_clean_gene(
      rownames(expr)[miss]
    )

    symbols[[key]] <- unique(
      sym[!is.na(sym)]
    )
  }

  common <- Reduce(intersect, symbols)

  ranked_genes <- tdrx_gnn_clean_gene(
    step5_result$ranked$Gene
  )

  ranked_genes <- unique(
    ranked_genes[!is.na(ranked_genes)]
  )

  genes <- head(
    ranked_genes[ranked_genes %in% common],
    candidate_limit
  )

  if (length(genes) < 3L) {
    stop(
      "Only ", length(genes),
      " Step 5 genes are present in every dataset. Need at least 3."
    )
  }

  # 6. Create temporary input files
  tmp_in <- tempfile("tdrx_gnn_input_")
  dir.create(tmp_in, recursive = TRUE)

  on.exit(
    unlink(tmp_in, recursive = TRUE, force = TRUE),
    add = TRUE
  )

  expr_files <- setNames(
    vector("list", length(keys)),
    keys
  )

  meta_files <- setNames(
    vector("list", length(keys)),
    keys
  )

  for (i in seq_along(keys)) {

    key <- keys[i]
    d <- deg_results[[key]]

    expr <- as.matrix(d$expression)
    storage.mode(expr) <- "double"

    if (any(!is.finite(expr))) {
      stop("Non-finite expression values in ", labels[i])
    }

    ann <- d$annotation

    ix <- match(
      rownames(expr),
      as.character(ann$gene_id)
    )

    sym <- tdrx_gnn_clean_gene(
      ann$gene_symbol[ix]
    )

    miss <- is.na(sym)

    sym[miss] <- tdrx_gnn_clean_gene(
      rownames(expr)[miss]
    )

    expr_df <- data.frame(
      probe_id = rownames(expr),
      gene_symbol = sym,
      as.data.frame(expr, check.names = FALSE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    md <- d$metadata

    if (!all(c("Samples", "Condition") %in% names(md))) {
      stop("Metadata must contain Samples and Condition for ", labels[i])
    }

    md_df <- data.frame(
      Samples = trimws(as.character(md$Samples)),
      Condition = trimws(as.character(md$Condition)),
      stringsAsFactors = FALSE
    )

    if (
      anyDuplicated(md_df$Samples) ||
      any(!md_df$Condition %in% c("Control", "Case"))
    ) {
      stop("Invalid sample IDs or Condition values for ", labels[i])
    }

    nm <- tdrx_gnn_safe_name(key)

    expr_files[[key]] <- file.path(
      tmp_in, paste0(nm, "_expression.csv")
    )

    meta_files[[key]] <- file.path(
      tmp_in, paste0(nm, "_metadata.csv")
    )

    utils::write.csv(
      expr_df,
      expr_files[[key]],
      row.names = FALSE,
      na = ""
    )

    utils::write.csv(
      md_df,
      meta_files[[key]],
      row.names = FALSE,
      na = ""
    )
  }

  candidate_file <- file.path(
    tmp_in, "TD_CCHI_candidates.csv"
  )

  utils::write.csv(
    data.frame(Gene = genes),
    candidate_file,
    row.names = FALSE
  )

  # 7. Pairwise dataset combinations
  pairs <- utils::combn(
    keys,
    2L,
    simplify = FALSE
  )

  n_pairs <- length(pairs)

  # 8. Easy analysis presets
  if (mode == "quick") {

    cfg <- list(
      bootstraps = 5L,
      warmup = 30L,
      cluster = 30L,
      patience = 10L
    )

  } else {

    cfg <- list(
      bootstraps = 100L,
      warmup = 250L,
      cluster = 200L,
      patience = 60L
    )
  }

  # 9. Persistent output directory
  run_dir <- file.path(
    output_root,
    paste0(
      "run_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  )

  dir.create(
    run_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::write.csv(
    data.frame(Gene = genes),
    file.path(run_dir, "TD_CCHI_candidates_used.csv"),
    row.names = FALSE
  )

  pair_tables <- list()
  errors <- list()

  # 10. Run original TD-CCHI model
  for (i in seq_along(pairs)) {

    a <- pairs[[i]][1L]
    b <- pairs[[i]][2L]

    la <- labels[match(a, keys)]
    lb <- labels[match(b, keys)]

    pair_name <- paste0(la, "_vs_", lb)

    pair_dir <- file.path(
      run_dir,
      tdrx_gnn_safe_name(pair_name)
    )

    dir.create(
      pair_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (is.function(progress)) {
      progress(
        (i - 1L) / n_pairs,
        paste("Running", i, "of", n_pairs, "·", pair_name)
      )
    }

    fit <- tryCatch(

      run_core(
        cohort1_expression = expr_files[[a]],
        cohort1_metadata = meta_files[[a]],

        cohort2_expression = expr_files[[b]],
        cohort2_metadata = meta_files[[b]],

        candidate_gene_file = candidate_file,

        cohort1_name = tdrx_gnn_safe_name(la),
        cohort2_name = tdrx_gnn_safe_name(lb),

        output_dir = pair_dir,

        bootstraps = cfg$bootstraps,
        top_fraction = 0.10,

        warmup_epochs = cfg$warmup,
        cluster_epochs = cfg$cluster,
        patience = cfg$patience,

        k_min = 2L,
        k_max = 8L,

        seed = 20260908L + i * 1000L,

        device = "auto",
        verbose = FALSE
      ),

      error = function(e) {
        errors[[pair_name]] <<- conditionMessage(e)
        NULL
      }
    )

    if (
      !is.null(fit) &&
      is.data.frame(fit$results) &&
      nrow(fit$results)
    ) {

      z <- fit$results

      z$Pair <- pair_name

      z$Pair_Top <- z$Rank <= max(
        1L,
        ceiling(nrow(z) * 0.10)
      )

      pair_tables[[pair_name]] <- z
    }

    if (is.function(progress)) {
      progress(
        i / n_pairs,
        paste("Finished", i, "of", n_pairs, "·", pair_name)
      )
    }
  }

  if (!length(pair_tables)) {

    details <- paste(
      paste(names(errors), unlist(errors), sep = ": "),
      collapse = " | "
    )

    stop("All GNN pair runs failed. ", details)
  }

  # 11. Combine pairwise GNN results
  pairwise <- do.call(rbind, pair_tables)
  rownames(pairwise) <- NULL

  groups <- split(pairwise, pairwise$Gene)

  agg <- do.call(
    rbind,
    lapply(groups, function(z) {

      avg <- function(x) {
        x <- x[is.finite(x)]
        if (length(x)) mean(x) else 0
      }

      spread <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) > 1L) stats::sd(x) else 0
      }

      data.frame(
        Gene = z$Gene[1L],

        S_GNN = avg(z$CCHS),
        GNN_SD = spread(z$CCHS),

        GNN_PairTopSupport = sum(
          z$Pair_Top,
          na.rm = TRUE
        ),

        GNN_PairCount = nrow(z),

        GNN_C = avg(z$Centrality_C),
        GNN_I = avg(z$InterCluster_I),
        GNN_S = avg(z$Stability_S),
        GNN_R = avg(z$Reproducibility_R),

        stringsAsFactors = FALSE
      )
    })
  )

  rownames(agg) <- NULL

  agg <- agg[
    order(
      -agg$S_GNN,
      -agg$GNN_PairTopSupport,
      agg$Gene
    ),
    ,
    drop = FALSE
  ]

  agg$GNN_Rank <- seq_len(nrow(agg))

  # 12. Integrate with Step 5 ranking
  ranked <- merge(
    step5_result$ranked,
    agg,
    by = "Gene",
    all = FALSE,
    sort = FALSE
  )

  if (!nrow(ranked)) {
    stop("No Step 5 candidates matched the GNN results.")
  }

  for (nm in c("S_DEG", "S_WGCNA", "S_GNN")) {

    if (!nm %in% names(ranked)) {
      ranked[[nm]] <- NA_real_
    }

    ranked[[nm]] <- as.numeric(ranked[[nm]])

    ranked[[nm]][!is.finite(ranked[[nm]])] <- NA_real_

    ranked[[nm]] <- pmax(
      0,
      pmin(1, ranked[[nm]])
    )
  }

  # Intermediate three-layer score, NOT final five-layer score
  ranked$Network_Combined_Score <- rowMeans(
    as.matrix(
      ranked[, c("S_DEG", "S_WGCNA", "S_GNN"), drop = FALSE]
    ),
    na.rm = TRUE
  )

  ranked$Network_Combined_Score[
    !is.finite(ranked$Network_Combined_Score)
  ] <- NA_real_

  ranked$Network_Rank <- rank(
    -ranked$Network_Combined_Score,
    ties.method = "min",
    na.last = "keep"
  )

  ranked <- ranked[
    order(
      ranked$Network_Rank,
      -ranked$S_GNN,
      ranked$Gene,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(ranked) <- NULL

  # 13. Top 20 genes for the later ML stage
  hubs_for_ml <- data.frame(
    Gene = head(
      agg$Gene,
      min(20L, nrow(agg))
    )
  )

  # 14. Save outputs
  utils::write.csv(
    pairwise,
    file.path(run_dir, "TD_CCHI_pairwise_long.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    agg,
    file.path(run_dir, "TD_CCHI_aggregate.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    ranked,
    file.path(run_dir, "TD_CCHI_network_ranked.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    hubs_for_ml,
    file.path(run_dir, "TD_CCHI_hubs_for_113ML.csv"),
    row.names = FALSE
  )

  audit <- list(
    network_source = "Absolute Spearman correlation from normalized gene expression",
    curated_ppi_database_used = FALSE,
    datasets = keys,
    candidate_genes = genes,
    attempted_pairs = n_pairs,
    successful_pairs = length(pair_tables),
    errors = errors,
    mode = mode,
    configuration = cfg,
    output_directory = normalizePath(
      run_dir,
      winslash = "/",
      mustWork = FALSE
    )
  )

  saveRDS(
    audit,
    file.path(run_dir, "TD_CCHI_analysis_audit.rds")
  )

  # 15. Top hub score plot
  top <- head(agg, 20L)

  score_plot <- ggplot2::ggplot(
    top,
    ggplot2::aes(
      x = reorder(Gene, S_GNN),
      y = S_GNN
    )
  ) +
    ggplot2::geom_col(
      fill = "#0F766E",
      width = 0.68
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Top TD-CCHI network hubs",
      x = NULL,
      y = "Mean Cross-Cluster Hub Score"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  # 16. GNN evidence component heatmap
  component_plot <- NULL

  if (nrow(top)) {

    comp <- data.frame(
      Gene = rep(top$Gene, 4L),

      Component = rep(
        c(
          "Centrality",
          "Inter-cluster",
          "Stability",
          "Reproducibility"
        ),
        each = nrow(top)
      ),

      Score = c(
        top$GNN_C,
        top$GNN_I,
        top$GNN_S,
        top$GNN_R
      )
    )

    comp$Gene <- factor(
      comp$Gene,
      levels = rev(top$Gene)
    )

    component_plot <- ggplot2::ggplot(
      comp,
      ggplot2::aes(
        x = Gene,
        y = Component,
        fill = Score
      )
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(
        low = "#E2E8F0",
        high = "#0F766E",
        limits = c(0, 1)
      ) +
      ggplot2::labs(
        title = "GNN evidence components",
        x = NULL,
        y = NULL,
        fill = "Score"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 35,
          hjust = 1
        )
      )
  }

  # 17. Return results to Shiny
  list(
    ranked = ranked,
    aggregate = agg,
    pairwise_long = pairwise,
    candidate_genes = genes,
    pair_count = n_pairs,
    successful_pairs = length(pair_tables),
    errors = errors,
    mode = mode,
    configuration = cfg,
    output_dir = run_dir,
    score_plot = score_plot,
    component_plot = component_plot,
    audit = audit
  )
}
