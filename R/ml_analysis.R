
# ============================================================
# TriDiseaseRank-X STEP 7
# Corrected GNN Top-20 + Original 113 ML Integration
# ============================================================

tdrx_ml_clean_gene <- function(x) {

  x <- trimws(as.character(x))

  x[x %in% c(
    "", "NA", "N/A", "NaN",
    "---", "--", "NULL", "."
  )] <- NA_character_

  x <- sub(
    "[[:space:]]*(///|//|;|\\|).*$",
    "",
    x
  )

  toupper(trimws(x))
}


tdrx_ml_safe_name <- function(x) {

  x <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    as.character(x)
  )

  ifelse(nzchar(x), x, "dataset")
}


# ------------------------------------------------------------
# 1. Find original ML script
# ------------------------------------------------------------

tdrx_find_original_ml <- function() {

  candidates <- c(
    "original/TriDiseaseRankX_ML.R",
    "original/TriDiseaseRankX_113ML_Internal.R",
    "TriDiseaseRankX_ML.R",
    "TriDiseaseRankX_113ML_Internal.R"
  )

  found <- candidates[file.exists(candidates)]

  if (!length(found)) {

    stop(
      "Original ML script not found. ",
      "Expected original/TriDiseaseRankX_ML.R"
    )
  }

  normalizePath(
    found[1L],
    winslash = "/",
    mustWork = TRUE
  )
}


# ------------------------------------------------------------
# 2. Main ML Analysis Function
# ------------------------------------------------------------

tdrx_ml_batch_analysis <- function(

    gnn_result,

    deg_results,

    dataset_names = NULL,

    mode = c("quick", "full"),

    output_root = file.path(
      "TriDiseaseRank_results",
      "Step7_113ML"
    ),

    progress = NULL,

    seed = 20260908L

) {

  mode <- match.arg(mode)

  # ----------------------------------------------------------
  # A. Check dependencies
  # ----------------------------------------------------------

  required_packages <- c(
    "caret", "pROC", "glmnet", "ranger",
    "mboost", "MASS", "plsRglm", "gbm",
    "xgboost", "naivebayes", "kernlab",
    "data.table", "matrixStats", "dplyr",
    "tidyr", "future", "future.apply",
    "ggplot2"
  )

  missing <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing)) {

    stop(
      "Missing ML packages: ",
      paste(missing, collapse = ", "),
      ". Run install_ml_packages.R first."
    )
  }

  # ----------------------------------------------------------
  # B. Validate GNN result
  # ----------------------------------------------------------

  if (
    is.null(gnn_result$aggregate) ||
    !is.data.frame(gnn_result$aggregate)
  ) {

    stop(
      "GNN aggregate result missing. ",
      "Complete Step 6 first."
    )
  }

  gnn <- as.data.frame(
    gnn_result$aggregate,
    stringsAsFactors = FALSE
  )

  if (!all(c("Gene", "S_GNN") %in% names(gnn))) {

    stop(
      "GNN aggregate requires Gene and S_GNN columns."
    )
  }

  # ----------------------------------------------------------
  # C. CORRECTED GNN TOP-20 SORTING
  # ----------------------------------------------------------

  gnn$Gene <- tdrx_ml_clean_gene(
    gnn$Gene
  )

  gnn$S_GNN <- suppressWarnings(
    as.numeric(gnn$S_GNN)
  )

  if ("GNN_PairTopSupport" %in% names(gnn)) {

    gnn$GNN_PairTopSupport <- suppressWarnings(
      as.numeric(gnn$GNN_PairTopSupport)
    )

  } else {

    gnn$GNN_PairTopSupport <- rep(
      0,
      nrow(gnn)
    )
  }

  gnn$GNN_PairTopSupport[
    !is.finite(gnn$GNN_PairTopSupport)
  ] <- 0

  gnn <- gnn[
    !is.na(gnn$Gene) &
      nzchar(gnn$Gene) &
      is.finite(gnn$S_GNN),
    ,
    drop = FALSE
  ]

  if (!nrow(gnn)) {

    stop("No valid GNN candidate genes found.")
  }

  # Primary: highest S_GNN
  # Secondary: highest PairTopSupport
  # Tertiary: alphabetical Gene symbol

  gnn <- gnn[
    order(
      -gnn$S_GNN,
      -gnn$GNN_PairTopSupport,
      gnn$Gene
    ),
    ,
    drop = FALSE
  ]

  # Remove duplicate genes AFTER sorting
  gnn <- gnn[
    !duplicated(gnn$Gene),
    ,
    drop = FALSE
  ]

  # Exactly the top 20 hub genes
  candidate_genes <- head(
    gnn$Gene,
    20L
  )

  if (length(candidate_genes) < 2L) {

    stop(
      "At least two valid GNN hub genes are required."
    )
  }

  # Save GNN ranking used for ML
  gnn_top20 <- gnn[
    match(candidate_genes, gnn$Gene),
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # D. Validate DEG datasets
  # ----------------------------------------------------------

  keys <- names(deg_results)

  if (
    length(keys) < 2L ||
    is.null(keys) ||
    anyNA(keys) ||
    any(!nzchar(keys)) ||
    anyDuplicated(keys)
  ) {

    stop(
      "ML requires at least two named DEG datasets."
    )
  }

  labels <- keys

  if (!is.null(dataset_names)) {

    supplied <- as.character(
      dataset_names[keys]
    )

    good <- !is.na(supplied) & nzchar(supplied)

    labels[good] <- supplied[good]
  }

  # ----------------------------------------------------------
  # E. Load original 113 ML engine
  # ----------------------------------------------------------

  ml_script <- tdrx_find_original_ml()

  ml_env <- new.env(
    parent = globalenv()
  )

  sys.source(
    ml_script,
    envir = ml_env
  )

  if (!exists(
    "model_definitions",
    envir = ml_env,
    inherits = FALSE
  )) {

    stop("Original ML model definitions not found.")
  }

  original_models <- get(
    "model_definitions",
    envir = ml_env,
    inherits = FALSE
  )

  if (length(original_models) != 113L) {

    stop(
      "Original ML engine must contain 113 models. Found: ",
      length(original_models)
    )
  }

  if (!exists(
    "run_113ml_internal",
    envir = ml_env,
    inherits = FALSE
  )) {

    stop("run_113ml_internal() not found in original script.")
  }

  # ----------------------------------------------------------
  # F. Quick mode uses 8 ORIGINAL model definitions
  # ----------------------------------------------------------

  quick_models <- c(
    "SVM",
    "Ridge",
    "Lasso",
    "RF",
    "GBM",
    "LDA",
    "XGBoost",
    "NaiveBayes"
  )

  if (!all(quick_models %in% original_models)) {

    stop(
      "Quick model names do not match original model definitions."
    )
  }

  selected_models <- if (mode == "quick") {
    quick_models
  } else {
    original_models
  }

  # Override only the in-memory model list.
  # Original source file remains unchanged.

  ml_env$model_definitions <- selected_models

  # ----------------------------------------------------------
  # G. Create reproducible run directory
  # ----------------------------------------------------------

  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  run_dir <- tempfile(
    pattern = paste0(
      "run_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_"
    ),
    tmpdir = output_root
  )

  dir.create(
    run_dir,
    recursive = TRUE
  )

  candidate_file <- file.path(
    run_dir,
    "GNN_Top20_Hub_Genes.csv"
  )

  utils::write.csv(
    gnn_top20,
    candidate_file,
    row.names = FALSE,
    na = ""
  )

  # ----------------------------------------------------------
  # H. Preserve and restore Future plan
  # ----------------------------------------------------------

  old_plan <- future::plan()

  on.exit(
    try(
      future::plan(old_plan),
      silent = TRUE
    ),
    add = TRUE
  )

  # ----------------------------------------------------------
  # I. Run original ML engine per dataset
  # ----------------------------------------------------------

  results <- list()
  errors <- list()

  total <- length(keys)

  for (i in seq_along(keys)) {

    key <- keys[i]
    label <- labels[i]

    if (is.function(progress)) {

      progress(
        (i - 1L) / total,
        paste("Preparing:", label)
      )
    }

    dataset_dir <- file.path(
      run_dir,
      paste0(
        sprintf("%02d", i),
        "_",
        tdrx_ml_safe_name(key)
      )
    )

    dir.create(
      dataset_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    # ----------------------------------------
    # Dataset-specific validation and execution
    # ----------------------------------------

    res <- tryCatch({

      d <- deg_results[[key]]

      if (
        is.null(d$expression) ||
        is.null(d$metadata) ||
        is.null(d$annotation)
      ) {

        stop(
          "DEG expression, metadata or annotation missing."
        )
      }

      expr <- as.matrix(
        d$expression
      )

      if (!is.numeric(expr)) {
        storage.mode(expr) <- "double"
      }

      if (
        is.null(rownames(expr)) ||
        is.null(colnames(expr)) ||
        any(!is.finite(expr))
      ) {

        stop("Invalid normalized expression matrix.")
      }

      md <- as.data.frame(
        d$metadata,
        stringsAsFactors = FALSE
      )

      if (!all(
        c("Samples", "Condition") %in% names(md)
      )) {

        stop(
          "Metadata requires Samples and Condition."
        )
      }

      sample_ids <- trimws(
        as.character(md$Samples)
      )

      condition <- trimws(
        as.character(md$Condition)
      )

      if (
        anyNA(sample_ids) ||
        any(!nzchar(sample_ids)) ||
        anyDuplicated(sample_ids)
      ) {

        stop("Invalid or duplicated sample IDs.")
      }

      if (!all(
        condition %in% c("Control", "Case")
      )) {

        stop(
          "Condition must contain exactly Control and Case."
        )
      }

      group_counts <- table(
        factor(
          condition,
          levels = c("Control", "Case")
        )
      )

      if (any(group_counts < 2L)) {

        stop(
          "At least two samples per class are required."
        )
      }

      missing_samples <- setdiff(
        sample_ids,
        colnames(expr)
      )

      if (length(missing_samples)) {

        stop(
          "Missing expression samples: ",
          paste(head(missing_samples, 5L), collapse = ", ")
        )
      }

      # Preserve metadata order
      expr <- expr[
        ,
        sample_ids,
        drop = FALSE
      ]

      # ----------------------------------------
      # Gene annotation
      # ----------------------------------------

      ann <- as.data.frame(
        d$annotation,
        stringsAsFactors = FALSE
      )

      if (!all(
        c("gene_id", "gene_symbol") %in% names(ann)
      )) {

        stop(
          "Annotation requires gene_id and gene_symbol."
        )
      }

      gene_ids <- rownames(expr)

      ann_index <- match(
        gene_ids,
        as.character(ann$gene_id)
      )

      symbols <- tdrx_ml_clean_gene(
        ann$gene_symbol[ann_index]
      )

      missing_symbols <- is.na(symbols)

      symbols[missing_symbols] <- tdrx_ml_clean_gene(
        gene_ids[missing_symbols]
      )

      # ----------------------------------------
      # Match GNN Top20 to dataset genes
      # ----------------------------------------

      available_symbols <- unique(
        symbols[!is.na(symbols)]
      )

      dataset_candidates <- intersect(
        candidate_genes,
        available_symbols
      )

      if (length(dataset_candidates) < 2L) {

        stop(
          "Fewer than two GNN Top20 genes are available in this dataset."
        )
      }

      dataset_candidate_file <- file.path(
        dataset_dir,
        "GNN_candidates_available.csv"
      )

      utils::write.csv(
        data.frame(Gene = dataset_candidates),
        dataset_candidate_file,
        row.names = FALSE
      )

      # ----------------------------------------
      # Prepare original ML input format
      # ----------------------------------------

      expression_df <- as.data.frame(
        expr,
        check.names = FALSE
      )

      expression_export <- data.frame(
        probe_id = gene_ids,
        gene_symbol = symbols,
        expression_df,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      expression_file <- file.path(
        dataset_dir,
        "normalized_expression.csv"
      )

      utils::write.csv(
        expression_export,
        expression_file,
        row.names = FALSE,
        na = ""
      )

      metadata_export <- data.frame(
        Samples = sample_ids,
        Condition = condition,
        stringsAsFactors = FALSE
      )

      # Preserve paired biological sample grouping
      if ("Pair_ID" %in% names(md)) {

        metadata_export$Pair_ID <- as.character(
          md$Pair_ID
        )
      }

      metadata_file <- file.path(
        dataset_dir,
        "metadata.csv"
      )

      utils::write.csv(
        metadata_export,
        metadata_file,
        row.names = FALSE,
        na = ""
      )

      # ----------------------------------------
      # Run ORIGINAL run_113ml_internal()
      # ----------------------------------------

      ml_output_dir <- file.path(
        dataset_dir,
        "ML_results"
      )

      future::plan(
        future::sequential
      )

      ml_env$run_113ml_internal(

        expression_file = expression_file,

        metadata_file = metadata_file,

        candidate_gene_file = dataset_candidate_file,

        dataset_name = tdrx_ml_safe_name(label),

        output_dir = ml_output_dir,

        outer_folds = 5L,

        inner_folds = 4L,

        tuning_len = 3L,

        workers = 1L,

        consensus_top_models = 20L,

        seed = as.integer(seed + i * 10000L)
      )

    }, error = function(e) {

      structure(
        list(message = conditionMessage(e)),
        class = "tdrx_ml_error"
      )
    })

    if (inherits(res, "tdrx_ml_error")) {

      errors[[key]] <- res$message

    } else {

      if (
        is.null(res$performance) ||
        is.null(res$gene_scores) ||
        is.null(res$best_model)
      ) {

        errors[[key]] <- "Original ML engine returned incomplete results."

      } else {

        results[[key]] <- res
      }
    }

    # Ensure future execution does not remain parallel
    future::plan(
      future::sequential
    )

    if (is.function(progress)) {

      progress(
        i / total,
        paste(
          "Completed:",
          length(results),
          "of",
          total,
          "datasets"
        )
      )
    }
  }

  # ----------------------------------------------------------
  # J. Stop if all datasets failed
  # ----------------------------------------------------------

  if (!length(results)) {

    error_text <- paste(
      paste(
        names(errors),
        unlist(errors),
        sep = ": "
      ),
      collapse = " | "
    )

    stop(
      "All ML datasets failed. ",
      error_text
    )
  }

  # ----------------------------------------------------------
  # K. Combine original model performance
  # ----------------------------------------------------------

  performance_list <- lapply(
    names(results),
    function(key) {

      z <- as.data.frame(
        results[[key]]$performance
      )

      z$Dataset <- labels[match(key, keys)]

      z
    }
  )

  performance_all <- do.call(
    rbind,
    performance_list
  )

  rownames(performance_all) <- NULL

  # ----------------------------------------------------------
  # L. Combine gene selection scores
  # ----------------------------------------------------------

  gene_score_list <- lapply(
    names(results),
    function(key) {

      z <- as.data.frame(
        results[[key]]$gene_scores
      )

      z$Dataset <- labels[match(key, keys)]

      z
    }
  )

  gene_scores_long <- do.call(
    rbind,
    gene_score_list
  )

  rownames(gene_scores_long) <- NULL

  # ----------------------------------------------------------
  # M. Combine best models
  # ----------------------------------------------------------

  best_model_list <- lapply(
    names(results),
    function(key) {

      z <- as.data.frame(
        results[[key]]$best_model
      )

      z$Dataset <- labels[match(key, keys)]

      z
    }
  )

  best_models <- do.call(
    rbind,
    best_model_list
  )

  rownames(best_models) <- NULL

  # ----------------------------------------------------------
  # N. Aggregate ML scores across datasets
  # ----------------------------------------------------------

  aggregate <- data.frame(
    Gene = candidate_genes,
    stringsAsFactors = FALSE
  )

  score_columns <- character(0)

  for (i in seq_along(keys)) {

    key <- keys[i]

    label <- labels[i]

    col_name <- paste0(
      "S_ML_",
      tdrx_ml_safe_name(label)
    )

    if (key %in% names(results)) {

      gs <- results[[key]]$gene_scores

      gs_gene <- tdrx_ml_clean_gene(
        gs$Gene
      )

      idx <- match(
        candidate_genes,
        gs_gene
      )

      values <- as.numeric(
        gs$S_ML_Disease[idx]
      )

    } else {

      values <- rep(
        NA_real_,
        length(candidate_genes)
      )
    }

    col_name <- make.unique(
      c(score_columns, col_name)
    )[length(score_columns) + 1L]

    aggregate[[col_name]] <- values

    score_columns <- c(
      score_columns,
      col_name
    )
  }

  score_matrix <- as.matrix(
    aggregate[
      ,
      score_columns,
      drop = FALSE
    ]
  )

  # Average only across datasets with available scores
  aggregate$S_ML <- rowMeans(
    score_matrix,
    na.rm = TRUE
  )

  aggregate$S_ML[
    !is.finite(aggregate$S_ML)
  ] <- NA_real_

  aggregate$ML_DatasetSupport <- rowSums(
    score_matrix > 0,
    na.rm = TRUE
  )

  aggregate$ML_DatasetCoverage <- rowSums(
    !is.na(score_matrix)
  )

  # Add original GNN evidence
  gnn_index <- match(
    aggregate$Gene,
    gnn$Gene
  )

  aggregate$S_GNN <- gnn$S_GNN[gnn_index]

  aggregate$GNN_PairTopSupport <-
    gnn$GNN_PairTopSupport[gnn_index]

  aggregate <- aggregate[
    order(
      -aggregate$S_ML,
      -aggregate$ML_DatasetSupport,
      -aggregate$S_GNN,
      aggregate$Gene,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  aggregate$ML_Rank <- seq_len(
    nrow(aggregate)
  )

  aggregate <- aggregate[
    ,
    c(
      "ML_Rank",
      setdiff(names(aggregate), "ML_Rank")
    ),
    drop = FALSE
  ]

  rownames(aggregate) <- NULL

  # ----------------------------------------------------------
  # O. Save all output files
  # ----------------------------------------------------------

  utils::write.csv(
    gnn_top20,
    file.path(
      run_dir,
      "GNN_Top20_used_for_ML.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    performance_all,
    file.path(
      run_dir,
      "113ML_model_performance_all_datasets.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    gene_scores_long,
    file.path(
      run_dir,
      "113ML_gene_scores_long.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    aggregate,
    file.path(
      run_dir,
      "113ML_aggregate_gene_scores.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    best_models,
    file.path(
      run_dir,
      "113ML_best_models_by_dataset.csv"
    ),
    row.names = FALSE
  )

  # Error report
  error_table <- if (length(errors)) {

    data.frame(
      Dataset = names(errors),
      Error = unlist(errors),
      stringsAsFactors = FALSE
    )

  } else {

    data.frame(
      Dataset = character(),
      Error = character()
    )
  }

  utils::write.csv(
    error_table,
    file.path(
      run_dir,
      "113ML_dataset_errors.csv"
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # P. Model performance visualization
  # ----------------------------------------------------------

  valid_perf <- performance_all[
    performance_all$Status == "OK" &
      is.finite(performance_all$Internal_AUC),
    ,
    drop = FALSE
  ]

  model_plot <- NULL

  if (nrow(valid_perf)) {

    split_perf <- split(
      valid_perf,
      valid_perf$Dataset
    )

    top_perf <- do.call(
      rbind,
      lapply(
        split_perf,
        function(z) {

          z <- z[
            order(
              -z$Internal_AUC,
              z$Model
            ),
            ,
            drop = FALSE
          ]

          head(z, 5L)
        }
      )
    )

    top_perf$Model_Label <- paste(
      top_perf$Dataset,
      top_perf$Model,
      sep = " | "
    )

    model_plot <- ggplot2::ggplot(
      top_perf,
      ggplot2::aes(
        x = stats::reorder(
          Model_Label,
          Internal_AUC
        ),
        y = Internal_AUC,
        fill = Dataset
      )
    ) +
      ggplot2::geom_col(
        width = 0.7
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Top ML model performance",
        x = NULL,
        y = "Internal out-of-fold AUC"
      ) +
      ggplot2::theme_minimal(
        base_size = 12
      )
  }

  # ----------------------------------------------------------
  # Q. ML-supported gene visualization
  # ----------------------------------------------------------

  plot_data <- aggregate[
    is.finite(aggregate$S_ML),
    ,
    drop = FALSE
  ]

  plot_data <- head(
    plot_data,
    20L
  )

  gene_plot <- NULL

  if (nrow(plot_data)) {

    gene_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = stats::reorder(
          Gene,
          S_ML
        ),
        y = S_ML
      )
    ) +
      ggplot2::geom_col(
        fill = "#0F766E",
        width = 0.7
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "GNN hub genes ranked by ML evidence",
        x = NULL,
        y = "Mean ML selection-stability score"
      ) +
      ggplot2::theme_minimal(
        base_size = 12
      )
  }

  # ----------------------------------------------------------
  # R. Experiment audit
  # ----------------------------------------------------------

  audit <- list(

    original_ml_script = ml_script,

    candidate_genes = candidate_genes,

    candidate_gene_count = length(candidate_genes),

    mode = mode,

    models_requested = length(selected_models),

    original_models_available = length(original_models),

    datasets_requested = keys,

    datasets_completed = names(results),

    dataset_errors = errors,

    outer_folds = 5L,

    inner_folds = 4L,

    tuning_length = 3L,

    workers = 1L,

    validation = "Internal nested cross-validation only",

    feature_selection_note =
      "GNN Top20 was selected before ML cross-validation.",

    output_directory = normalizePath(
      run_dir,
      winslash = "/",
      mustWork = FALSE
    )
  )

  saveRDS(
    audit,
    file.path(
      run_dir,
      "Step7_ML_audit.rds"
    )
  )

  # ----------------------------------------------------------
  # S. Return results to Shiny
  # ----------------------------------------------------------

  list(

    candidate_genes = candidate_genes,

    gnn_top20 = gnn_top20,

    mode = mode,

    models_requested = length(selected_models),

    models_available = length(original_models),

    successful_datasets = length(results),

    requested_datasets = length(keys),

    dataset_results = results,

    dataset_errors = errors,

    model_performance = performance_all,

    gene_scores_long = gene_scores_long,

    gene_scores = aggregate,

    best_models = best_models,

    model_plot = model_plot,

    gene_plot = gene_plot,

    output_dir = run_dir,

    audit = audit
  )
}
