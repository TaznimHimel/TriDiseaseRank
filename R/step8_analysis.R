
# ============================================================
# TriDiseaseRank-X STEP 8
# Original Pathway Engine + Five-Layer Evidence Fusion
# ============================================================

tdrx_step8_clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NAN", "NULL", "---", "--", ".")] <- NA_character_
  x
}

tdrx_step8_clamp01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- 0
  pmax(0, pmin(1, x))
}

tdrx_step8_find_pathway_script <- function() {

  candidates <- c(
    "original/Pathway.R",
    "Pathway.R",
    "original/TriDiseaseRankX_Pathway.R",
    "TriDiseaseRankX_Pathway.R"
  )

  found <- candidates[file.exists(candidates)]

  if (!length(found)) {
    stop(
      "Original Pathway.R not found. ",
      "Check original/Pathway.R"
    )
  }

  normalizePath(
    found[1L],
    winslash = "/",
    mustWork = TRUE
  )
}


# ============================================================
# MAIN STEP 8 FUNCTION
# ============================================================

tdrx_step8_analysis <- function(

    deg_results,
    wgcna_results,
    cross_result,
    gnn_result,
    ml_result,

    dataset_names = NULL,
    top_n = 10L,

    output_root = file.path(
      "TriDiseaseRank_results",
      "Step8_Final_Ranking"
    ),

    progress = NULL

) {

  # ----------------------------------------------------------
  # 1. Required packages
  # ----------------------------------------------------------

  required <- c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "dplyr",
    "tidyr",
    "data.table",
    "ggplot2",
    "openxlsx"
  )

  missing <- required[
    !vapply(
      required,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing)) {
    stop(
      "Missing core packages: ",
      paste(missing, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 2. Validate previous stages
  # ----------------------------------------------------------

  if (
    is.null(cross_result$ranked) ||
    is.null(gnn_result$aggregate) ||
    is.null(ml_result$gene_scores)
  ) {
    stop(
      "Complete Step 5, Step 6 and Step 7 before Step 8."
    )
  }

  deg_keys <- names(deg_results)

  if (
    is.null(deg_keys) ||
    length(deg_keys) < 2L
  ) {
    stop("At least two DEG datasets are required.")
  }

  pathway_keys <- deg_keys[
    deg_keys %in% names(wgcna_results)
  ]

  if (length(pathway_keys) < 2L) {
    stop(
      "At least two completed WGCNA datasets are required."
    )
  }

  labels <- deg_keys

  if (!is.null(dataset_names)) {

    supplied <- as.character(
      dataset_names[deg_keys]
    )

    good <- !is.na(supplied) & nzchar(supplied)

    labels[good] <- supplied[good]
  }

  # ----------------------------------------------------------
  # 3. Step 5 candidate genes
  # ----------------------------------------------------------

  candidate_table <- as.data.frame(
    cross_result$ranked,
    stringsAsFactors = FALSE
  )

  required_candidate_columns <- c(
    "Gene", "S_DEG", "S_WGCNA"
  )

  if (!all(
    required_candidate_columns %in%
      names(candidate_table)
  )) {
    stop(
      "Step 5 table requires Gene, S_DEG and S_WGCNA."
    )
  }

  candidate_table$Gene <- tdrx_step8_clean_symbol(
    candidate_table$Gene
  )

  candidate_table <- candidate_table[
    !is.na(candidate_table$Gene),
    ,
    drop = FALSE
  ]

  candidate_table <- candidate_table[
    !duplicated(candidate_table$Gene),
    ,
    drop = FALSE
  ]

  candidate_genes <- unique(
    candidate_table$Gene
  )

  if (!length(candidate_genes)) {
    stop("No candidate genes available from Step 5.")
  }

  # ----------------------------------------------------------
  # 4. Prepare disease gene lists and backgrounds
  # ----------------------------------------------------------

  disease_gene_lists <- setNames(
    vector("list", length(pathway_keys)),
    pathway_keys
  )

  universe_gene_lists <- setNames(
    vector("list", length(pathway_keys)),
    pathway_keys
  )

  pathway_input_counts <- vector(
    "list",
    length(pathway_keys)
  )

  names(pathway_input_counts) <- pathway_keys

  for (key in pathway_keys) {

    de <- deg_results[[key]]
    wg <- wgcna_results[[key]]

    # Significant DEG
    de_sig <- de$significant

    if (
      is.null(de_sig) &&
      !is.null(de$all) &&
      "Significant" %in% names(de$all)
    ) {

      de_sig <- de$all[
        !is.na(de$all$Significant) &
          de$all$Significant,
        ,
        drop = FALSE
      ]
    }

    if (
      is.null(de_sig) ||
      !"gene_symbol" %in% names(de_sig)
    ) {
      stop(
        "Significant DEG gene symbols missing: ",
        key
      )
    }

    # Significant WGCNA module genes
    wg_sig <- wg$significant_genes

    if (is.null(wg_sig)) {
      wg_sig <- wg$sig_module_genes
    }

    if (
      is.null(wg_sig) ||
      !"gene_symbol" %in% names(wg_sig)
    ) {
      stop(
        "WGCNA significant genes missing: ",
        key
      )
    }

    de_symbols <- unique(
      tdrx_step8_clean_symbol(
        de_sig$gene_symbol
      )
    )

    wg_symbols <- unique(
      tdrx_step8_clean_symbol(
        wg_sig$gene_symbol
      )
    )

    de_symbols <- de_symbols[
      !is.na(de_symbols)
    ]

    wg_symbols <- wg_symbols[
      !is.na(wg_symbols)
    ]

    # Original research strategy:
    # Significant DEG intersect significant WGCNA genes
    disease_gene_lists[[key]] <- intersect(
      de_symbols,
      wg_symbols
    )

    # Dataset-specific tested gene universe
    tested <- if (
      !is.null(de$all) &&
      "gene_symbol" %in% names(de$all)
    ) {
      de$all$gene_symbol
    } else {
      de$annotation$gene_symbol
    }

    tested <- unique(
      tdrx_step8_clean_symbol(tested)
    )

    tested <- tested[
      !is.na(tested)
    ]

    if (!length(tested)) {
      stop("No tested gene universe found: ", key)
    }

    universe_gene_lists[[key]] <- tested

    pathway_input_counts[[key]] <- data.frame(
      Dataset = labels[match(key, deg_keys)],
      Significant_DEG = length(de_symbols),
      Significant_WGCNA_genes = length(wg_symbols),
      DEG_WGCNA_intersection =
        length(disease_gene_lists[[key]]),
      Tested_gene_universe = length(tested),
      stringsAsFactors = FALSE
    )
  }

  pathway_input_counts <- do.call(
    rbind,
    pathway_input_counts
  )

  rownames(pathway_input_counts) <- NULL

  if (is.function(progress)) {
    progress(0.08, "Preparing pathway gene sets")
  }

  # ----------------------------------------------------------
  # 5. Load ORIGINAL Pathway.R
  # ----------------------------------------------------------

  suppressPackageStartupMessages(
    library(dplyr)
  )

  pathway_script <- tdrx_step8_find_pathway_script()

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
    recursive = TRUE,
    showWarnings = FALSE
  )

  pathway_dir <- file.path(
    run_dir,
    "Pathway"
  )

  dir.create(
    pathway_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ----------------------------------------------------------
  # 6. Safe compatibility wrapper for original Pathway.R
  # ----------------------------------------------------------

  # Empty enrichment tables should not break CSV export.
  pathway_env <- new.env(
    parent = globalenv()
  )

  pathway_env$tdrx_safe_fwrite <- function(x, file) {

    if (
      is.data.frame(x) &&
      ncol(x) > 0L
    ) {

      data.table::fwrite(
        x,
        file
      )

    } else {

      writeLines(
        "No significant enrichment rows.",
        file
      )
    }
  }

  pathway_source <- readLines(
    pathway_script,
    warn = FALSE,
    encoding = "UTF-8"
  )

  # Original script maps SYMBOL to Entrez IDs.
  # Use NCBI gene IDs for KEGG input.
  pathway_source <- gsub(
    'keyType = "kegg"',
    'keyType = "ncbi-geneid"',
    pathway_source,
    fixed = TRUE
  )

  # Guard empty GO/KEGG/optional enrichment tables.
  pathway_source <- gsub(
    "data.table::fwrite(tab, file.path",
    "tdrx_safe_fwrite(tab, file.path",
    pathway_source,
    fixed = TRUE
  )

  patched_pathway <- file.path(
    pathway_dir,
    "Pathway_web_compatible.R"
  )

  writeLines(
    pathway_source,
    patched_pathway,
    useBytes = TRUE
  )

  invisible(
    parse(file = patched_pathway)
  )

  sys.source(
    patched_pathway,
    envir = pathway_env
  )

  if (!exists(
    "run_tdrx_pathway",
    envir = pathway_env,
    inherits = FALSE
  )) {
    stop("run_tdrx_pathway() not found.")
  }

  # ----------------------------------------------------------
  # 7. Run original Pathway Analysis
  # ----------------------------------------------------------

  if (is.function(progress)) {
    progress(
      0.15,
      "Running GO, KEGG and available optional databases"
    )
  }

  pathway_result <- pathway_env$run_tdrx_pathway(

    disease_gene_lists = disease_gene_lists,

    universe_gene_lists = universe_gene_lists,

    candidate_genes = candidate_genes,

    output_dir = pathway_dir,

    padj_cutoff = 0.05,

    pvalue_cutoff = 0.05,

    qvalue_cutoff = 0.20,

    min_gs_size = 5L,

    max_gs_size = 500L,

    run_reactome = requireNamespace(
      "ReactomePA",
      quietly = TRUE
    ),

    run_wikipathways = requireNamespace(
      "msigdbr",
      quietly = TRUE
    ),

    run_disease_ontology = requireNamespace(
      "DOSE",
      quietly = TRUE
    ),

    top_terms_plot = 20L
  )

  if (is.function(progress)) {
    progress(
      0.78,
      "Combining five evidence layers"
    )
  }

  # ----------------------------------------------------------
  # 8. Prepare DEG + WGCNA evidence
  # ----------------------------------------------------------

  support_col <- if (
    "WGCNA_Support" %in% names(candidate_table)
  ) {
    "WGCNA_Support"
  } else if (
    "WGCNA_support_count" %in% names(candidate_table)
  ) {
    "WGCNA_support_count"
  } else {
    NA_character_
  }

  candidate_cols <- c(
    "Gene", "S_DEG", "S_WGCNA"
  )

  if (!is.na(support_col)) {
    candidate_cols <- c(
      candidate_cols,
      support_col
    )
  }

  candidate_layer <- candidate_table[
    ,
    candidate_cols,
    drop = FALSE
  ]

  if (
    !is.na(support_col) &&
    support_col != "WGCNA_support_count"
  ) {

    names(candidate_layer)[
      names(candidate_layer) == support_col
    ] <- "WGCNA_support_count"
  }

  # ----------------------------------------------------------
  # 9. Prepare GNN evidence
  # ----------------------------------------------------------

  gnn <- as.data.frame(
    gnn_result$aggregate,
    stringsAsFactors = FALSE
  )

  if (!all(
    c("Gene", "S_GNN") %in% names(gnn)
  )) {
    stop("GNN aggregate requires Gene and S_GNN.")
  }

  gnn$Gene <- tdrx_step8_clean_symbol(
    gnn$Gene
  )

  gnn <- gnn[
    !is.na(gnn$Gene) &
      !duplicated(gnn$Gene),
    ,
    drop = FALSE
  ]

  gnn_cols <- intersect(
    c(
      "Gene", "S_GNN", "GNN_SD",
      "GNN_PairTopSupport", "GNN_PairCount",
      "GNN_C", "GNN_I", "GNN_S", "GNN_R"
    ),
    names(gnn)
  )

  gnn_layer <- gnn[
    ,
    gnn_cols,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 10. Prepare ML evidence
  # ----------------------------------------------------------

  ml <- as.data.frame(
    ml_result$gene_scores,
    stringsAsFactors = FALSE
  )

  if (!all(
    c("Gene", "S_ML") %in% names(ml)
  )) {
    stop(
      "Step 7 ML aggregate must contain Gene and S_ML."
    )
  }

  ml$Gene <- tdrx_step8_clean_symbol(
    ml$Gene
  )

  ml <- ml[
    !is.na(ml$Gene) &
      !duplicated(ml$Gene),
    ,
    drop = FALSE
  ]

  ml_cols <- intersect(
    c(
      "Gene",
      "S_ML",
      "ML_DatasetSupport",
      "ML_DatasetCoverage"
    ),
    names(ml)
  )

  ml_layer <- ml[
    ,
    ml_cols,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 11. Prepare Pathway evidence
  # ----------------------------------------------------------

  pathway_scores <- as.data.frame(
    pathway_result$gene_scores,
    stringsAsFactors = FALSE
  )

  if (!all(
    c("Gene", "S_Pathway") %in%
      names(pathway_scores)
  )) {
    stop("Pathway engine did not return S_Pathway.")
  }

  pathway_scores$Gene <- tdrx_step8_clean_symbol(
    pathway_scores$Gene
  )

  pathway_scores <- pathway_scores[
    !is.na(pathway_scores$Gene) &
      !duplicated(pathway_scores$Gene),
    ,
    drop = FALSE
  ]

  pathway_cols <- intersect(
    c(
      "Gene",
      "S_Pathway",
      "PathwayCount",
      "DatasetPathwaySupport",
      "DatabaseSupport",
      "ConvergentPathwayCount",
      "SharedAllDatasetPathwayCount",
      "BestAdjustedP",
      "SupportedDatasets",
      "SupportedDatabases"
    ),
    names(pathway_scores)
  )

  pathway_layer <- pathway_scores[
    ,
    pathway_cols,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 12. Five-layer integration
  # ----------------------------------------------------------

  base <- merge(
    candidate_layer,
    gnn_layer,
    by = "Gene",
    all = FALSE,
    sort = FALSE
  )

  base <- merge(
    base,
    ml_layer,
    by = "Gene",
    all = FALSE,
    sort = FALSE
  )

  base <- merge(
    base,
    pathway_layer,
    by = "Gene",
    all = FALSE,
    sort = FALSE
  )

  if (!nrow(base)) {
    stop(
      "No genes matched across all five evidence layers."
    )
  }

  layers <- c(
    "DEG",
    "WGCNA",
    "GNN",
    "ML",
    "Pathway"
  )

  score_columns <- paste0(
    "S_",
    layers
  )

  for (nm in score_columns) {

    if (!nm %in% names(base)) {
      stop("Missing evidence layer: ", nm)
    }

    base[[nm]] <- tdrx_step8_clamp01(
      base[[nm]]
    )
  }

  # Original research default:
  # 20% weight per evidence layer.

  primary_weights <- c(
    DEG = 0.20,
    WGCNA = 0.20,
    GNN = 0.20,
    ML = 0.20,
    Pathway = 0.20
  )

  for (layer in layers) {

    base[[paste0(
      "Contribution_",
      layer
    )]] <-
      primary_weights[layer] *
      base[[paste0("S_", layer)]]
  }

  contribution_cols <- paste0(
    "Contribution_",
    layers
  )

  base$TriDiseaseRankX_Score <- rowSums(
    base[, contribution_cols, drop = FALSE]
  )

  # ----------------------------------------------------------
  # 13. Rank fusion sensitivity
  # ----------------------------------------------------------

  rank_pct <- function(x) {

    if (length(x) <= 1L) {
      return(rep(1, length(x)))
    }

    1 - (
      rank(-x, ties.method = "average") - 1
    ) / (length(x) - 1)
  }

  rank_matrix <- do.call(
    cbind,
    lapply(
      score_columns,
      function(nm) rank_pct(base[[nm]])
    )
  )

  colnames(rank_matrix) <- layers

  base$RankFusion_Sensitivity <- rowMeans(
    rank_matrix
  )

  score_matrix <- as.matrix(
    base[, score_columns, drop = FALSE]
  )

  base$Strongest_Evidence <- layers[
    max.col(
      score_matrix,
      ties.method = "first"
    )
  ]

  base$Weakest_Evidence <- layers[
    max.col(
      -score_matrix,
      ties.method = "first"
    )
  ]

  # ----------------------------------------------------------
  # 14. Research validation status
  # ----------------------------------------------------------

  ml_mode <- if (!is.null(ml_result$mode)) {
    as.character(ml_result$mode[1L])
  } else {
    "unknown"
  }

  gnn_complete <- !is.null(gnn_result$pair_count) &&
    !is.null(gnn_result$successful_pairs) &&
    gnn_result$successful_pairs >= gnn_result$pair_count

  ml_complete <- FALSE

  if (
    identical(ml_mode, "full") &&
    identical(as.integer(ml_result$models_requested), 113L) &&
    ml_result$successful_datasets >= length(deg_keys)
  ) {

    perf <- ml_result$model_performance

    if (
      is.data.frame(perf) &&
      all(c("Dataset", "Status") %in% names(perf))
    ) {

      ml_complete <- all(vapply(
        labels,
        function(label) {

          z <- perf[
            perf$Dataset == label,
            ,
            drop = FALSE
          ]

          nrow(z) == 113L &&
            all(z$Status == "OK")
        },
        logical(1)
      ))
    }
  }

  all_pathway_datasets <- length(pathway_keys) ==
    length(deg_keys)

  ranking_status <- if (
    ml_complete &&
    gnn_complete &&
    all_pathway_datasets
  ) {

    "FULL_113ML_INTERNAL_VALIDATION"

  } else if (
    identical(ml_mode, "quick")
  ) {

    "PRELIMINARY_QUICK_ML"

  } else {

    "PRELIMINARY_PARTIAL_COVERAGE"
  }

  base$Ranking_Status <- ranking_status

  base$External_Validation_Status <- "NOT_PERFORMED"

  # ----------------------------------------------------------
  # 15. Final ranking
  # ----------------------------------------------------------

  base <- base[
    order(
      -base$TriDiseaseRankX_Score,
      -base$RankFusion_Sensitivity,
      base$Gene
    ),
    ,
    drop = FALSE
  ]

  rownames(base) <- NULL

  base$Final_Rank <- seq_len(
    nrow(base)
  )

  base <- base[
    ,
    c(
      "Final_Rank",
      setdiff(names(base), "Final_Rank")
    ),
    drop = FALSE
  ]

  top_n <- as.integer(top_n)

  if (
    length(top_n) != 1L ||
    is.na(top_n) ||
    top_n < 1L ||
    top_n > 100L
  ) {
    stop("Top biomarker count must be between 1 and 100.")
  }

  final_panel <- utils::head(
    base,
    min(top_n, nrow(base))
  )

  # ----------------------------------------------------------
  # 16. Save all result tables
  # ----------------------------------------------------------

  pathway_enrichment <- pathway_result$enrichment
  pathway_convergence <- pathway_result$convergence

  if (!is.data.frame(pathway_enrichment)) {
    pathway_enrichment <- data.frame()
  }

  if (!is.data.frame(pathway_convergence)) {
    pathway_convergence <- data.frame()
  }

  ranking_file <- file.path(
    run_dir,
    "TriDiseaseRankX_five_layer_ranking.csv"
  )

  top_file <- file.path(
    run_dir,
    "TriDiseaseRankX_top_biomarkers.csv"
  )

  pathway_score_file <- file.path(
    run_dir,
    "TriDiseaseRankX_pathway_gene_scores.csv"
  )

  enrichment_file <- file.path(
    run_dir,
    "TriDiseaseRankX_pathway_enrichment.csv"
  )

  convergence_file <- file.path(
    run_dir,
    "TriDiseaseRankX_pathway_convergence.csv"
  )

  weights_file <- file.path(
    run_dir,
    "TriDiseaseRankX_weights.csv"
  )

  panel_file <- file.path(
    run_dir,
    "TriDiseaseRankX_panel_for_external_validation.csv"
  )

  utils::write.csv(base, ranking_file, row.names = FALSE, na = "")
  utils::write.csv(final_panel, top_file, row.names = FALSE, na = "")
  utils::write.csv(pathway_scores, pathway_score_file, row.names = FALSE, na = "")
  utils::write.csv(pathway_enrichment, enrichment_file, row.names = FALSE, na = "")
  utils::write.csv(pathway_convergence, convergence_file, row.names = FALSE, na = "")

  utils::write.csv(
    data.frame(
      Layer = names(primary_weights),
      Weight = as.numeric(primary_weights)
    ),
    weights_file,
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(Gene = final_panel$Gene),
    panel_file,
    row.names = FALSE
  )

  utils::write.csv(
    pathway_input_counts,
    file.path(run_dir, "Pathway_input_gene_counts.csv"),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 17. Final ranking plot
  # ----------------------------------------------------------

  final_score_plot <- ggplot2::ggplot(
    utils::head(base, 20L),
    ggplot2::aes(
      x = stats::reorder(
        Gene,
        TriDiseaseRankX_Score
      ),
      y = TriDiseaseRankX_Score
    )
  ) +
    ggplot2::geom_col(
      fill = "#0F766E",
      width = 0.68
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Top TriDiseaseRank-X biomarkers",
      x = NULL,
      y = "Five-layer weighted evidence score"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  # ----------------------------------------------------------
  # 18. Five-layer evidence heatmap
  # ----------------------------------------------------------

  top_evidence <- utils::head(base, 15L)

  evidence_long <- data.frame(
    Gene = rep(
      top_evidence$Gene,
      times = length(layers)
    ),
    Layer = rep(
      layers,
      each = nrow(top_evidence)
    ),
    Score = unlist(
      lapply(
        score_columns,
        function(nm) top_evidence[[nm]]
      )
    ),
    stringsAsFactors = FALSE
  )

  evidence_long$Gene <- factor(
    evidence_long$Gene,
    levels = rev(top_evidence$Gene)
  )

  evidence_long$Layer <- factor(
    evidence_long$Layer,
    levels = layers
  )

  evidence_plot <- ggplot2::ggplot(
    evidence_long,
    ggplot2::aes(
      x = Layer,
      y = Gene,
      fill = Score
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_gradient(
      low = "#E2E8F0",
      high = "#0F766E",
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Five-layer evidence profile",
      x = NULL,
      y = NULL,
      fill = "Score"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  # ----------------------------------------------------------
  # 19. Pathway convergence plot
  # ----------------------------------------------------------

  pathway_plot <- NULL

  if (
    nrow(pathway_enrichment) > 0L &&
    all(
      c(
        "Pathway_Key",
        "Pathway_Description",
        "Database",
        "Dataset",
        "p.adjust"
      ) %in% names(pathway_enrichment)
    )
  ) {

    pathway_keys_unique <- unique(
      pathway_enrichment$Pathway_Key
    )

    pathway_plot_data <- do.call(
      rbind,
      lapply(
        pathway_keys_unique,
        function(k) {

          z <- pathway_enrichment[
            pathway_enrichment$Pathway_Key == k,
            ,
            drop = FALSE
          ]

          data.frame(
            Pathway_Key = k,
            Pathway_Description = as.character(
              z$Pathway_Description[1L]
            ),
            Database = as.character(
              z$Database[1L]
            ),
            DatasetSupport = length(
              unique(z$Dataset)
            ),
            MinPadj = min(
              as.numeric(z$p.adjust),
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )

    pathway_plot_data <- pathway_plot_data[
      is.finite(pathway_plot_data$MinPadj),
      ,
      drop = FALSE
    ]

    pathway_plot_data <- pathway_plot_data[
      order(
        -pathway_plot_data$DatasetSupport,
        pathway_plot_data$MinPadj
      ),
      ,
      drop = FALSE
    ]

    pathway_plot_data <- utils::head(
      pathway_plot_data,
      15L
    )

    if (nrow(pathway_plot_data)) {

      pathway_plot_data$Label <- factor(
        pathway_plot_data$Pathway_Description,
        levels = rev(unique(
          pathway_plot_data$Pathway_Description
        ))
      )

      pathway_plot <- ggplot2::ggplot(
        pathway_plot_data,
        ggplot2::aes(
          x = -log10(pmax(MinPadj, 1e-300)),
          y = Label
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            size = DatasetSupport,
            shape = Database
          ),
          alpha = 0.85
        ) +
        ggplot2::labs(
          title = "Top convergent pathways",
          x = "-log10 adjusted P-value",
          y = NULL,
          size = "Dataset support",
          shape = "Database"
        ) +
        ggplot2::theme_minimal(base_size = 11)
    }
  }

  # ----------------------------------------------------------
  # 20. Save graphs
  # ----------------------------------------------------------

  ggplot2::ggsave(
    file.path(run_dir, "TriDiseaseRankX_final_ranking.png"),
    final_score_plot,
    width = 9,
    height = 6,
    dpi = 220
  )

  ggplot2::ggsave(
    file.path(run_dir, "TriDiseaseRankX_evidence_heatmap.png"),
    evidence_plot,
    width = 9,
    height = 7,
    dpi = 220
  )

  if (!is.null(pathway_plot)) {

    ggplot2::ggsave(
      file.path(run_dir, "TriDiseaseRankX_pathway_convergence.png"),
      pathway_plot,
      width = 10,
      height = 7,
      dpi = 220
    )
  }

  # ----------------------------------------------------------
  # 21. Excel report
  # ----------------------------------------------------------

  safe_sheet <- function(x, message) {

    if (
      is.data.frame(x) &&
      ncol(x) > 0L
    ) {
      return(x)
    }

    data.frame(
      Message = message,
      stringsAsFactors = FALSE
    )
  }

  ml_success <- ml_result$successful_datasets
  ml_requested <- ml_result$models_requested

  run_metadata <- data.frame(
    Metric = c(
      "Ranking status",
      "ML mode",
      "ML models requested per dataset",
      "ML datasets completed",
      "DEG datasets",
      "Pathway datasets",
      "Candidate genes",
      "Final ranked genes",
      "Significant pathway rows",
      "External validation"
    ),
    Value = c(
      ranking_status,
      ml_mode,
      as.character(ml_requested),
      as.character(ml_success),
      as.character(length(deg_keys)),
      as.character(length(pathway_keys)),
      as.character(length(candidate_genes)),
      as.character(nrow(base)),
      as.character(nrow(pathway_enrichment)),
      "Not performed"
    ),
    stringsAsFactors = FALSE
  )

  workbook_file <- file.path(
    run_dir,
    "TriDiseaseRankX_results.xlsx"
  )

  openxlsx::write.xlsx(
    list(
      Final_Ranking = base,
      Top_Biomarkers = final_panel,
      DEG_WGCNA_Candidates = candidate_layer,
      GNN_Scores = gnn_layer,
      ML_Scores = ml_layer,
      Pathway_Gene_Scores = pathway_scores,
      Pathway_Convergence = safe_sheet(
        pathway_convergence,
        "No significant pathway convergence."
      ),
      Pathway_Enrichment = safe_sheet(
        pathway_enrichment,
        "No significant enrichment terms."
      ),
      Pathway_Input_Counts = pathway_input_counts,
      Evidence_Weights = data.frame(
        Layer = names(primary_weights),
        Weight = as.numeric(primary_weights)
      ),
      Run_Metadata = run_metadata
    ),
    file = workbook_file,
    overwrite = TRUE
  )

  capture.output(
    utils::sessionInfo(),
    file = file.path(run_dir, "sessionInfo.txt")
  )

  saveRDS(
    list(
      pathway_result = pathway_result,
      final_ranking = base,
      audit = run_metadata
    ),
    file.path(run_dir, "Step8_complete_results.rds")
  )

  if (is.function(progress)) {
    progress(1, "Five-layer ranking completed")
  }

  # ----------------------------------------------------------
  # 22. Return to Shiny
  # ----------------------------------------------------------

  list(
    final_ranking = base,
    top_biomarkers = final_panel,

    pathway_result = pathway_result,
    pathway_input_counts = pathway_input_counts,

    weights = primary_weights,
    ranking_status = ranking_status,
    ml_mode = ml_mode,

    final_score_plot = final_score_plot,
    evidence_plot = evidence_plot,
    pathway_plot = pathway_plot,

    output_dir = run_dir,

    ranking_file = ranking_file,
    top_file = top_file,
    pathway_score_file = pathway_score_file,
    enrichment_file = enrichment_file,
    convergence_file = convergence_file,
    weights_file = weights_file,
    panel_file = panel_file,
    workbook_file = workbook_file
  )
}
