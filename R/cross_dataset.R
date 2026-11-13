
# TriDiseaseRank-X Step 5
# Cross-Dataset DEG + WGCNA Integration
# Uses Step 4 analysis results without rerunning models.

tdrx_cross_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))

  x[x %in% c(
    "", "NA", "N/A", "--", "---", ".", "NULL"
  )] <- NA_character_

  x
}

tdrx_cross_analysis <- function(
    deg_results,
    wgcna_results = list(),
    dataset_names = NULL,
    min_deg_support = 2L,
    min_wgcna_support = 0L,
    require_same_direction = FALSE,
    deg_weight = 0.50
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install the 'ggplot2' package first.")
  }

  keys <- names(deg_results)

  if (
    length(keys) < 2L ||
    anyNA(keys) ||
    any(!nzchar(keys)) ||
    anyDuplicated(keys)
  ) {
    stop(
      "Step 5 requires DEG results from at least TWO named datasets."
    )
  }

  n <- length(keys)

  valid_wg <- intersect(keys, names(wgcna_results))

  valid_wg <- valid_wg[
    vapply(valid_wg, function(k) {
      is.list(wgcna_results[[k]]) &&
        is.data.frame(
          wgcna_results[[k]]$significant_genes
        )
    }, logical(1))
  ]

  m <- length(valid_wg)

  min_deg_support <- as.integer(min_deg_support)
  min_wgcna_support <- as.integer(min_wgcna_support)

  if (
    length(min_deg_support) != 1L ||
    is.na(min_deg_support) ||
    min_deg_support < 2L ||
    min_deg_support > n
  ) {
    stop(
      "DEG support must be between 2 and the number of completed DEG datasets (",
      n, ")."
    )
  }

  if (
    length(min_wgcna_support) != 1L ||
    is.na(min_wgcna_support) ||
    min_wgcna_support < 0L ||
    min_wgcna_support > m
  ) {
    stop(
      "WGCNA support must be between 0 and the number of completed WGCNA datasets (",
      m, ")."
    )
  }

  if (
    length(deg_weight) != 1L ||
    !is.finite(deg_weight) ||
    deg_weight < 0 ||
    deg_weight > 1
  ) {
    stop("DEG score weight must be between 0 and 1.")
  }

  if (
    !is.logical(require_same_direction) ||
    length(require_same_direction) != 1L ||
    is.na(require_same_direction)
  ) {
    stop("The same-direction option must be TRUE or FALSE.")
  }

  # Human-readable dataset names
  labels <- keys

  if (!is.null(dataset_names)) {
    v <- as.character(dataset_names[keys])
    good <- !is.na(v) & nzchar(v)
    labels[good] <- v[good]
  }

  suffixes <- make.unique(
    gsub("[^A-Za-z0-9_]", "_", labels)
  )

  de_tables <- vector("list", n)
  names(de_tables) <- keys

  all_genes <- character(0)

  # -----------------------------------------
  # 1. Collect DEG results from all datasets
  # -----------------------------------------

  for (i in seq_along(keys)) {

    d <- deg_results[[keys[i]]]$all

    need <- c(
      "gene_symbol",
      "logFC",
      "S_DEG",
      "Significant",
      "padj"
    )

    if (
      !is.data.frame(d) ||
      !all(need %in% names(d))
    ) {
      stop(
        "DEG result for ", labels[i],
        " does not contain: ",
        paste(need, collapse = ", ")
      )
    }

    z <- data.frame(
      Gene = tdrx_cross_gene(d$gene_symbol),
      logFC = as.numeric(d$logFC),
      FDR = as.numeric(d$padj),
      DEG_Score = as.numeric(d$S_DEG),
      DEG_Significant =
        !is.na(d$Significant) &
        as.logical(d$Significant),
      stringsAsFactors = FALSE
    )

    z <- z[!is.na(z$Gene), , drop = FALSE]

    z$DEG_Score[
      !is.finite(z$DEG_Score)
    ] <- 0

    z$DEG_Score <- pmax(
      0, pmin(1, z$DEG_Score)
    )

    z <- z[
      order(
        z$FDR,
        -abs(z$logFC),
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    z <- z[
      !duplicated(z$Gene),
      ,
      drop = FALSE
    ]

    if (!nrow(z)) {
      stop("No usable gene symbols found for ", labels[i])
    }

    de_tables[[keys[i]]] <- z

    all_genes <- union(
      all_genes, z$Gene
    )
  }

  # -----------------------------------------
  # 2. Join DEG tables by gene symbol
  # -----------------------------------------

  out <- data.frame(
    Gene = sort(all_genes),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(keys)) {

    z <- de_tables[[keys[i]]]

    names(z)[-1] <- paste0(
      names(z)[-1], "_", suffixes[i]
    )

    out <- merge(
      out, z,
      by = "Gene",
      all.x = TRUE,
      sort = FALSE
    )
  }

  # -----------------------------------------
  # 3. Add significant WGCNA module evidence
  # -----------------------------------------

  for (i in seq_along(keys)) {

    key <- keys[i]

    score_col <- paste0(
      "WGCNA_Score_", suffixes[i]
    )

    module_col <- paste0(
      "WGCNA_Module_", suffixes[i]
    )

    out[[score_col]] <- NA_real_
    out[[module_col]] <- NA_character_

    if (!key %in% valid_wg) next

    w <- wgcna_results[[key]]$significant_genes

    if (!all(
      c("gene_symbol", "S_WGCNA", "Module") %in%
        names(w)
    )) {
      stop(
        "WGCNA result for ",
        labels[i],
        " has missing columns."
      )
    }

    if (!nrow(w)) next

    wtab <- data.frame(
      Gene = tdrx_cross_gene(w$gene_symbol),
      Score = as.numeric(w$S_WGCNA),
      Module = as.character(w$Module),
      stringsAsFactors = FALSE
    )

    wtab$Score[
      !is.finite(wtab$Score)
    ] <- NA_real_

    wtab <- wtab[
      !is.na(wtab$Gene),
      ,
      drop = FALSE
    ]

    wtab <- wtab[
      order(
        -wtab$Score,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    wtab <- wtab[
      !duplicated(wtab$Gene),
      ,
      drop = FALSE
    ]

    pos <- match(out$Gene, wtab$Gene)

    out[[score_col]] <- wtab$Score[pos]
    out[[module_col]] <- wtab$Module[pos]
  }

  # -----------------------------------------
  # 4. Calculate dataset support
  # -----------------------------------------

  sig_columns <- paste0(
    "DEG_Significant_", suffixes
  )

  lfc_columns <- paste0(
    "logFC_", suffixes
  )

  ds_columns <- paste0(
    "DEG_Score_", suffixes
  )

  ws_columns <- paste0(
    "WGCNA_Score_", suffixes
  )

  sig_matrix <- as.matrix(
    out[, sig_columns, drop = FALSE]
  )

  sig_matrix[
    is.na(sig_matrix)
  ] <- FALSE

  storage.mode(sig_matrix) <- "logical"

  out$DEG_Support <- rowSums(sig_matrix)

  wg_matrix <- as.matrix(
    out[, ws_columns, drop = FALSE]
  )

  out$WGCNA_Support <- rowSums(
    !is.na(wg_matrix)
  )

  wg_matrix[
    !is.finite(wg_matrix)
  ] <- 0

  wg_matrix[wg_matrix < 0] <- 0
  wg_matrix[wg_matrix > 1] <- 1

  # -----------------------------------------
  # 5. Mean normalized evidence scores
  # -----------------------------------------

  deg_matrix <- as.matrix(
    out[, ds_columns, drop = FALSE]
  )

  deg_matrix[
    !is.finite(deg_matrix)
  ] <- 0

  out$S_DEG <- rowMeans(deg_matrix)

  if (m > 0L) {

    out$S_WGCNA <- rowMeans(
      wg_matrix[
        ,
        keys %in% valid_wg,
        drop = FALSE
      ]
    )

  } else {

    out$S_WGCNA <- rep(
      NA_real_, nrow(out)
    )
  }

  # Direction must agree across all
  # completed DEG datasets.

  lfc_matrix <- as.matrix(
    out[, lfc_columns, drop = FALSE]
  )

  out$Same_direction <- apply(
    lfc_matrix, 1L,
    function(z) {

      if (
        any(!is.finite(z)) ||
        any(z == 0)
      ) {
        return(FALSE)
      }

      all(z > 0) || all(z < 0)
    }
  )

  # -----------------------------------------
  # 6. Combined gene ranking
  # -----------------------------------------

  effective_deg_weight <- if (m > 0L) {
    deg_weight
  } else {
    1
  }

  effective_wg_weight <- 1 -
    effective_deg_weight

  out$Combined_Score <-
    effective_deg_weight * out$S_DEG +
    effective_wg_weight * (
      if (m > 0L) out$S_WGCNA else 0
    )

  out$Combined_Score[
    !is.finite(out$Combined_Score)
  ] <- 0

  ranked <- out[
    out$DEG_Support >= min_deg_support &
      out$WGCNA_Support >= min_wgcna_support &
      (
        !require_same_direction |
          out$Same_direction
      ),
    ,
    drop = FALSE
  ]

  ranked <- ranked[
    order(
      -ranked$Combined_Score,
      -ranked$DEG_Support,
      -ranked$WGCNA_Support,
      ranked$Gene
    ),
    ,
    drop = FALSE
  ]

  ranked$Rank <- seq_len(nrow(ranked))

  ranked <- ranked[
    ,
    c(
      "Rank",
      setdiff(names(ranked), "Rank")
    ),
    drop = FALSE
  ]

  # -----------------------------------------
  # 7. Pairwise DEG overlaps and Jaccard
  # -----------------------------------------

  pairs <- utils::combn(
    seq_along(keys), 2L
  )

  pairwise <- do.call(
    rbind,
    lapply(
      seq_len(ncol(pairs)),
      function(j) {

        a <- pairs[1, j]
        b <- pairs[2, j]

        ga <- de_tables[[keys[a]]]$Gene[
          de_tables[[keys[a]]]$DEG_Significant
        ]

        gb <- de_tables[[keys[b]]]$Gene[
          de_tables[[keys[b]]]$DEG_Significant
        ]

        un <- union(ga, gb)

        overlap <- length(
          intersect(ga, gb)
        )

        data.frame(
          Dataset_A = labels[a],
          Dataset_B = labels[b],
          Significant_A = length(ga),
          Significant_B = length(gb),
          Shared_DEGs = overlap,
          Jaccard = if (length(un)) {
            overlap / length(un)
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
    )
  )

  # -----------------------------------------
  # 8. DEG support distribution plot
  # -----------------------------------------

  freq <- as.data.frame(
    table(
      factor(
        out$DEG_Support[
          out$DEG_Support > 0
        ],
        levels = seq_len(n)
      )
    )
  )

  names(freq) <- c("Datasets", "Genes")

  freq$Datasets <- as.integer(
    as.character(freq$Datasets)
  )

  support_plot <- ggplot2::ggplot(
    freq,
    ggplot2::aes(
      x = factor(Datasets),
      y = Genes
    )
  ) +
    ggplot2::geom_col(
      fill = "#0F766E",
      width = 0.65
    ) +
    ggplot2::labs(
      title = "Significant DEG support across datasets",
      x = "Number of datasets",
      y = "Gene count"
    ) +
    ggplot2::theme_minimal(
      base_size = 12
    )

  # -----------------------------------------
  # 9. Top ranked gene heatmap
  # -----------------------------------------

  heatmap_plot <- NULL

  if (nrow(ranked)) {

    top <- utils::head(
      ranked, 20L
    )

    heat <- do.call(
      rbind,
      lapply(
        seq_along(keys),
        function(i) {

          data.frame(
            Gene = top$Gene,
            Dataset = labels[i],
            logFC = top[[
              lfc_columns[i]
            ]],
            stringsAsFactors = FALSE
          )
        }
      )
    )

    heat$Gene <- factor(
      heat$Gene,
      levels = rev(top$Gene)
    )

    heat$Dataset <- factor(
      heat$Dataset,
      levels = labels
    )

    heatmap_plot <- ggplot2::ggplot(
      heat,
      ggplot2::aes(
        x = Dataset,
        y = Gene,
        fill = logFC
      )
    ) +
      ggplot2::geom_tile(
        color = "white",
        linewidth = 0.5
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2563EB",
        mid = "#F8FAFC",
        high = "#DC2626",
        midpoint = 0,
        na.value = "#E2E8F0"
      ) +
      ggplot2::labs(
        title = "Top ranked genes: log2 fold change",
        x = "Dataset",
        y = NULL,
        fill = "log2FC"
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 25,
          hjust = 1
        )
      )
  }

  # -----------------------------------------
  # 10. Return all results
  # -----------------------------------------

  list(
    ranked = ranked,
    universe = out,
    pairwise = pairwise,
    dataset_count = n,
    wgcna_count = m,

    shared_all = sum(
      out$DEG_Support == n
    ),

    shared_at_least = sum(
      out$DEG_Support >= min_deg_support
    ),

    deg_weight = effective_deg_weight,
    wgcna_weight = effective_wg_weight,

    min_deg_support = min_deg_support,
    min_wgcna_support = min_wgcna_support,

    require_same_direction =
      require_same_direction,

    support_plot = support_plot,
    heatmap_plot = heatmap_plot
  )
}
