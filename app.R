# TriDiseaseRank-X | Step 4: multi-dataset DEG + WGCNA research workspace
# Original research scripts remain untouched; these are reusable web modules.
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
})
source(file.path("R", "batch_files.R"), local = TRUE)
source(file.path("R", "deg_analysis.R"), local = TRUE)
source(file.path("R", "wgcna_analysis.R"), local = TRUE)
# STEP5_INTEGRATION_V1
source(file.path("R", "cross_dataset.R"), local = TRUE)
# STEP6_GNN_INTEGRATION_V1
source(file.path("R", "gnn_analysis.R"), local = TRUE)
# STEP7_ML_INTEGRATION_V1
source(file.path("R", "ml_analysis.R"), local = TRUE)
# STEP8_PATHWAY_FUSION_V1
source(file.path("R", "step8_analysis.R"), local = TRUE)

# These helpers are loaded from the local R modules above at app startup.
utils::globalVariables(c(
  "tdrx_read_pair", "tdrx_prepare_input", "tdrx_deg_analysis",
  "tdrx_wgcna_analysis", "tdrx_cross_analysis", "tdrx_safe_name"
))

# WGCNA_COR_CONFLICT_FIX
suppressPackageStartupMessages(library(WGCNA))
cor <- WGCNA::cor

options(shiny.maxRequestSize = 250 * 1024^2)

metric_tile <- function(label, value, caption, icon) {
  tags$div(class = "metric-card",
    tags$div(class = "metric-head", tags$span(label), tags$span(class = "metric-icon", icon)),
    tags$div(class = "metric-value", value),
    tags$div(class = "metric-caption", caption)
  )
}

section_title <- function(eyebrow, title, caption = NULL) {
  tags$div(class = "section-heading",
    tags$div(class = "section-eyebrow", eyebrow),
    tags$h3(title),
    if (!is.null(caption)) tags$p(caption)
  )
}

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bg = "#F6F8FC", fg = "#172237", primary = "#0F766E"),
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  tags$div(class = "tdrx-layout",
    tags$aside(class = "tdrx-sidebar",
      tags$div(class = "logo-wrap",
        tags$div(class = "logo-mark", "T", tags$span("X")),
        tags$div(tags$strong("TriDiseaseRank-X"), tags$small("RESEARCH WORKSPACE"))
      ),
      tags$div(class = "sidebar-section-label", "01 · INPUT DATA"),
      tags$div(class = "side-panel",
        fileInput("counts_files", "Counts matrices · CSV",
                  multiple = TRUE, accept = ".csv", buttonLabel = "Choose files",
                  placeholder = "Select one or more CSVs"),
        fileInput("metadata_files", "Matching metadata · CSV",
                  multiple = TRUE, accept = ".csv", buttonLabel = "Choose files",
                  placeholder = "Select one or more CSVs"),
        tags$p(class = "hint", "Match by prefix: GSE120852_counts_data.csv + GSE120852_metadata.csv"),
        actionButton("validate_pairs", "Validate all datasets", icon = icon("check-circle"),
                     class = "btn-secondary action-wide")
      ),
      tags$div(class = "sidebar-section-label", "02 · DEG SETTINGS"),
      tags$div(class = "side-panel",
        numericInput("fdr", "DEG adjusted P (FDR)", value = 0.05,
                     min = 0.0001, max = 0.99, step = 0.01),
        numericInput("logfc", "Minimum absolute log2FC", value = 1,
                     min = 0, max = 10, step = 0.25),
        actionButton("run_deg", "Run DEG · all datasets", icon = icon("play"),
                     class = "btn-primary action-wide")
      ),
      tags$div(class = "sidebar-section-label", "03 · WGCNA SETTINGS"),
      tags$div(class = "side-panel",
        numericInput("top_genes", "Top variable genes", value = 2000L,
                     min = 100L, max = 5000L, step = 500L),
        numericInput("min_module", "Minimum module size", value = 30L,
                     min = 20L, max = 200L, step = 10L),
        numericInput("module_p", "Module association P-value", value = 0.05,
                     min = 0.0001, max = 0.99, step = 0.01),
        actionButton("run_wgcna", "Run WGCNA · DEG datasets", icon = icon("network-wired"),
                     class = "btn-primary action-wide")
      ),
      tags$div(class = "sidebar-section-label", "04 · CROSS-DATASET SETTINGS"),
      tags$div(class = "side-panel",
        numericInput("min_deg_support", "Minimum datasets with significant DEG",
                     value = 2L, min = 2L, max = 20L, step = 1L),
        numericInput("min_wgcna_support", "Minimum WGCNA-supporting datasets",
                     value = 0L, min = 0L, max = 20L, step = 1L),
        checkboxInput("same_direction", "Require the same log2FC direction", FALSE),
        sliderInput("deg_weight", "DEG contribution to ranking",
                    min = 0, max = 1, value = 0.5, step = 0.05),
        tags$p(class = "hint",
          "DEG support = 3 requires significant genes in all three datasets. WGCNA support = 3 requires support in three WGCNA datasets."),
        actionButton("run_cross", "Run Step 5 · Cross-Dataset",
                     icon = icon("layer-group"), class = "btn-primary action-wide")
      ),
      tags$div(class = "sidebar-section-label", "05 · NETWORK ANALYSIS"),
      tags$div(class = "side-panel",
        selectInput("gnn_mode", "Analysis mode",
          choices = c("Quick preview · recommended" = "quick",
                     "Full research run · slower" = "research"),
          selected = "quick", width = "100%"),
        numericInput("gnn_candidates", "Top Step 5 genes",
          value = 40L, min = 3L, max = 100L, step = 5L),
        tags$p(class = "hint",
          "Uses genes present in every completed DEG dataset. More genes and full mode need more CPU/RAM."),
        actionButton("run_gnn", "Run GNN · all dataset pairs",
          icon = icon("network-wired"), class = "btn-primary action-wide")
      ),
      tags$div(class = "sidebar-section-label", "07 · MACHINE LEARNING"),
      tags$div(class = "side-panel",
        selectInput(
          "ml_mode",
          "Analysis mode",
          choices = c(
            "Quick test · 8 models" = "quick",
            "Full research · 113 models" = "full"
          ),
          selected = "quick",
          width = "100%"
        ),
        tags$p(
          class = "hint",
          "The Top 20 hub genes from GNN are selected automatically."
        ),
        actionButton(
          "run_ml",
          "Run ML · GNN Top 20",
          icon = icon("play"),
          class = "btn-primary action-wide"
        )
      ),
      tags$div(
        class = "sidebar-section-label",
        "08 · FINAL DISCOVERY"
      ),

      tags$div(
        class = "side-panel",

        numericInput(
          "step8_top_n",
          "Top biomarkers to report",
          value = 10L,
          min = 5L,
          max = 50L,
          step = 5L
        ),

        tags$p(
          class = "hint",
          "Combines DEG, WGCNA, GNN, ML and Pathway evidence."
        ),

        actionButton(
          "run_step8",
          "Run Step 8 · Final ranking",
          icon = icon("award"),
          class = "btn-primary action-wide"
        )
      ),
      tags$div(class = "sidebar-footer",
        tags$span(class = "status-dot"), " Local analysis · no cloud upload"
      )
    ),
    tags$main(class = "tdrx-content",
      tags$div(class = "topline",
        tags$span(class = "breadcrumbs", "Projects", tags$span(class = "sep", "/"), "TriDiseaseRank-X"),
        tags$span(class = "version-pill", "BUILD 08 · FIVE-LAYER RANKING")
      ),
      tags$div(class = "hero",
        tags$div(class = "hero-copy",
          tags$span(class = "hero-kicker", "TRANSCRIPTOMICS ANALYSIS"),
          tags$h1("From raw counts to", tags$span(" meaningful signals.")),
          tags$p("Upload multiple RNA-seq datasets, compare differential expression, and explore gene co-expression modules from one workspace.")
        ),
        tags$div(class = "hero-art", tags$div(class = "art-circle a1"),
                 tags$div(class = "art-circle a2"), tags$div(class = "art-circle a3"))
      ),
      uiOutput("overview_cards"),
      tags$div(class = "workspace-controls",
        tags$div(class = "dataset-select",
          selectInput("dataset_id", "CURRENT DATASET", choices = character(0), width = "100%")
        ),
        tags$div(class = "control-note",
          tags$span(class = "status-dot teal"),
          "Batch runs process matching datasets independently."
        )
      ),
      tags$div(class = "workspace-tabs",
        tabsetPanel(id = "result_tabs", type = "tabs",
          tabPanel(title = "Datasets", value = "datasets",
            section_title("INPUT CHECKS", "Your dataset library",
              "Counts and metadata are paired using the dataset name before their filename suffix."),
            DTOutput("paired_table"),
            tags$div(class = "content-grid two",
              tags$div(class = "content-card",
                tags$h4("Counts preview"), DTOutput("counts_preview")
              ),
              tags$div(class = "content-card",
                tags$h4("Metadata preview"), DTOutput("metadata_preview")
              )
            ),
            tags$div(class = "content-card",
              tags$h4("Batch execution status"), DTOutput("execution_table")
            )
          ),
          tabPanel(title = "DEG Analysis", value = "deg",
            section_title("STEP 03", "Differential expression",
              "TMM normalization, limma-voom and moderated statistics (Case vs Control)."),
            uiOutput("deg_status"),
            uiOutput("deg_metrics"),
            tags$div(class = "content-card",
              tags$h4("Gene-level results"), DTOutput("deg_table")
            ),
            tags$div(class = "content-grid two",
              tags$div(class = "content-card", tags$h4("Volcano plot"),
                       plotOutput("volcano_plot", height = "360px")),
              tags$div(class = "content-card", tags$h4("PCA · voom expression"),
                       plotOutput("pca_plot", height = "360px"))
            )
          ),
          tabPanel(title = "WGCNA · Step 4", value = "wgcna",
            section_title("STEP 04", "Gene co-expression networks",
              "Signed WGCNA network, selected power, module assignment and module–trait relationships."),
            uiOutput("wgcna_status"),
            uiOutput("wgcna_metrics"),
            tags$div(class = "content-grid two",
              tags$div(class = "content-card", tags$h4("Soft-threshold fit"),
                       plotOutput("power_plot", height = "360px")),
              tags$div(class = "content-card", tags$h4("Module–trait associations"),
                       plotOutput("trait_plot", height = "360px"))
            ),
            tags$div(class = "content-card", tags$h4("Module statistics"),
                     DTOutput("modules_table")),
            tags$div(class = "content-card", tags$h4("WGCNA gene ranking"),
                     DTOutput("wgcna_gene_table"))
          ),
          tabPanel(title = "Cross-Dataset · Step 5", value = "cross",
            section_title("STEP 05", "Shared DEG + WGCNA ranking",
              "Compare shared DEGs, inspect overlaps, and rank combined evidence."),
            uiOutput("cross_status"),
            uiOutput("cross_metrics"),
            tags$div(class = "content-card",
              tags$h4("Ranked candidates (first 2,000 displayed)"),
              DTOutput("cross_table")
            ),
            tags$div(class = "content-grid two",
              tags$div(class = "content-card",
                tags$h4("DEG support distribution"),
                plotOutput("cross_support_plot", height = "350px")),
              tags$div(class = "content-card",
                tags$h4("Top gene log2FC heatmap"),
                plotOutput("cross_heatmap", height = "450px"))
            ),
            tags$div(class = "content-card",
              tags$h4("Pairwise significant DEG overlap"),
              DTOutput("cross_pairwise")
            ),
            tags$div(class = "cross-actions",
              downloadButton("download_cross", "Download ranked genes",
                             class = "btn-primary"),
              downloadButton("download_cross_all", "Download complete gene evidence",
                             class = "btn-outline-primary"),
              downloadButton("download_cross_pairs", "Download pairwise overlap",
                             class = "btn-outline-primary")
            )
          ),
          tabPanel(title = "Network Ranking · Step 6", value = "gnn",
            section_title("STEP 06", "Network-based gene ranking",
              "TD-CCHI runs pairwise across datasets and combines network evidence with Step 5 scores."),
            uiOutput("gnn_status"),
            uiOutput("gnn_metrics"),

            tags$div(class = "content-card",
              tags$h4("Network-ranked candidate genes · first 1,000 displayed"),
              DTOutput("gnn_table")
            ),

            tags$div(class = "content-grid two",
              tags$div(class = "content-card",
                tags$h4("Top network hubs"),
                plotOutput("gnn_score_plot", height = "390px")),
              tags$div(class = "content-card",
                tags$h4("GNN evidence components"),
                plotOutput("gnn_component_plot", height = "390px"))
            ),

            tags$div(class = "content-card",
              tags$h4("Pairwise GNN results"),
              DTOutput("gnn_pairwise_table")
            ),

            tags$div(class = "cross-actions",
              downloadButton("download_gnn", "Download network ranking",
                class = "btn-primary"),
              downloadButton("download_gnn_aggregate", "Download GNN scores",
                class = "btn-outline-primary"),
              downloadButton("download_gnn_pairs", "Download pairwise scores",
                class = "btn-outline-primary")
            ),

            tags$div(class = "research-note",
              tags$strong("Network source. "),
              "The supplied PPI_GNN_R.R builds gene-gene relationships from expression Spearman correlations. It does not currently load a curated STRING, BioGRID or IntAct PPI database."
            )
          ),
          tabPanel(
            title = "Machine Learning · Step 7",
            value = "ml",

            section_title(
              "STEP 07",
              "GNN-guided machine learning",
              "Compare 113 ML pipelines and identify stable candidate genes."
            ),

            uiOutput("ml_status"),
            uiOutput("ml_metrics"),

            tags$div(
              class = "content-card",
              tags$h4("Best model by dataset"),
              DTOutput("ml_best_table")
            ),

            tags$div(
              class = "content-grid two",

              tags$div(
                class = "content-card",
                tags$h4("Model performance"),
                plotOutput("ml_auc_plot", height = "430px")
              ),

              tags$div(
                class = "content-card",
                tags$h4("ML-supported hub genes"),
                plotOutput("ml_gene_plot", height = "430px")
              )
            ),

            tags$div(
              class = "content-card",
              tags$h4("All model performance results"),
              DTOutput("ml_model_table")
            ),

            tags$div(
              class = "content-card",
              tags$h4("Combined ML gene ranking"),
              DTOutput("ml_gene_table")
            ),

            tags$div(
              class = "cross-actions",

              downloadButton(
                "download_ml_performance",
                "Download model performance",
                class = "btn-primary"
              ),

              downloadButton(
                "download_ml_genes",
                "Download ML gene scores",
                class = "btn-outline-primary"
              ),

              downloadButton(
                "download_ml_best",
                "Download best models",
                class = "btn-outline-primary"
              )
            ),

            tags$div(
              class = "research-note",
              tags$strong("Validation note. "),
              "AUC represents internal out-of-fold performance. ",
              "GNN candidate selection was performed before ML cross-validation; ",
              "independent external validation is still required."
            )
          ),
          tabPanel(
            title = "Final Ranking · Step 8",
            value = "step8",

            section_title(
              "STEP 08",
              "Final biomarker discovery",
              "Pathway enrichment and five-layer evidence integration."
            ),

            uiOutput("step8_status"),
            uiOutput("step8_metrics"),

            tags$div(
              class = "content-card",
              tags$h4("Top prioritized biomarkers"),
              DTOutput("step8_top_table")
            ),

            tags$div(
              class = "content-grid two",

              tags$div(
                class = "content-card",
                tags$h4("Final ranking"),
                plotOutput("step8_final_plot", height = "390px")
              ),

              tags$div(
                class = "content-card",
                tags$h4("Five-layer evidence profile"),
                plotOutput("step8_evidence_plot", height = "430px")
              )
            ),

            tags$div(
              class = "content-card",
              tags$h4("Pathway convergence"),
              plotOutput("step8_pathway_plot", height = "420px")
            ),

            tags$div(
              class = "content-card",
              tags$h4("Complete five-layer ranking"),
              DTOutput("step8_final_table")
            ),

            tags$div(
              class = "content-card",
              tags$h4("Pathway enrichment"),
              DTOutput("step8_pathway_table")
            ),

            tags$div(
              class = "content-card",
              tags$h4("Cross-dataset pathway convergence"),
              DTOutput("step8_convergence_table")
            ),

            tags$div(
              class = "cross-actions",

              downloadButton(
                "download_step8_final",
                "Download final ranking",
                class = "btn-primary"
              ),

              downloadButton(
                "download_step8_top",
                "Download top biomarkers",
                class = "btn-outline-primary"
              ),

              downloadButton(
                "download_step8_pathway",
                "Download pathway results",
                class = "btn-outline-primary"
              ),

              downloadButton(
                "download_step8_workbook",
                "Download Excel report",
                class = "btn-outline-primary"
              )
            ),

            tags$div(
              class = "research-note",
              tags$strong("Research note. "),
              "Quick ML produces a preliminary ranking. ",
              "Full 113-model validation and independent external ",
              "validation are needed for publication-level claims."
            )
          ),
          tabPanel(title = "Exports", value = "exports",
            section_title("RESULTS LIBRARY", "Export your analysis",
              "Download selected dataset tables or package all completed analyses into one ZIP."),
            tags$div(class = "download-grid",
              tags$div(class = "download-card", tags$div(class = "download-symbol", "CSV"),
                       tags$h4("Selected DEG · All genes"),
                       tags$p("Statistics, fold changes and DEG scores."),
                       downloadButton("download_deg", "Download CSV", class = "btn-outline-primary")),
              tags$div(class = "download-card", tags$div(class = "download-symbol", "CSV"),
                       tags$h4("Selected DEG · Significant"),
                       tags$p("Genes passing selected FDR and log2FC cutoffs."),
                       downloadButton("download_deg_sig", "Download CSV", class = "btn-outline-primary")),
              tags$div(class = "download-card", tags$div(class = "download-symbol", "CSV"),
                       tags$h4("Selected WGCNA · Genes"),
                       tags$p("Module membership, gene significance and score."),
                       downloadButton("download_wgcna", "Download CSV", class = "btn-outline-primary")),
              tags$div(class = "download-card", tags$div(class = "download-symbol", "ZIP"),
                       tags$h4("All completed results"),
                       tags$p("A folder per dataset with tables and scientific plots."),
                       downloadButton("download_zip", "Download ZIP", class = "btn-primary"))
            )
          )
        )
      ),
      tags$div(class = "research-note",
        tags$strong("Research note. "),
        "WGCNA module significance uses nominal P < cutoff (the original project convention); FDR is also displayed. ",
        "Paired tissue, covariates, SVA and cohort validation are not adjusted in this prototype. ",
        "Step 5 is DEG + WGCNA evidence only; Step 6 adds TD-CCHI GNN. ML and pathway stages are not included yet."
      )
    )
  )
)

server <- function(input, output, session) {
  state <- reactiveValues(
    step8 = NULL,
    step8_error = NULL,
    ml = NULL,
    ml_error = NULL,
    gnn = NULL,
    gnn_error = NULL,
    cross = NULL,
    checks = list(), deg = list(), wgcna = list(),
    deg_errors = list(), wgcna_errors = list()
  )

  pairs <- reactive({
    tdrx_upload_manifest(input$counts_files, input$metadata_files)
  })
  ready_pairs <- reactive({
    x <- pairs()
    x[x$Status == "Ready", , drop = FALSE]
  })

  observeEvent(list(input$counts_files, input$metadata_files), {
    state$cross <- NULL
    state$checks <- list()
    state$deg <- list()
    state$wgcna <- list()
    state$deg_errors <- list()
    state$wgcna_errors <- list()
  }, ignoreInit = TRUE)

  observeEvent(list(input$fdr, input$logfc), {
    state$cross <- NULL
    state$deg <- list()
    state$wgcna <- list()
    state$deg_errors <- list()
    state$wgcna_errors <- list()
  }, ignoreInit = TRUE)

  observeEvent(list(input$top_genes, input$min_module, input$module_p), {
    state$cross <- NULL
    state$wgcna <- list()
    state$wgcna_errors <- list()
  }, ignoreInit = TRUE)

  observeEvent(list(input$min_deg_support, input$min_wgcna_support,
                    input$same_direction, input$deg_weight), {
    state$cross <- NULL
  }, ignoreInit = TRUE)

  observeEvent(list(input$counts_files, input$metadata_files), {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(list(input$fdr, input$logfc), {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(list(input$top_genes, input$min_module, input$module_p), {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(list(input$min_deg_support, input$min_wgcna_support,
                    input$same_direction, input$deg_weight), {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(list(input$gnn_mode, input$gnn_candidates), {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_deg, {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_wgcna, {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_cross, {
    state$gnn <- NULL
    state$gnn_error <- NULL
  }, ignoreInit = TRUE)
  observeEvent(
    list(input$counts_files, input$metadata_files),
    {
      state$ml <- NULL
      state$ml_error <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(
    list(input$fdr, input$logfc),
    {
      state$ml <- NULL
      state$ml_error <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(
    list(input$top_genes, input$min_module, input$module_p),
    {
      state$ml <- NULL
      state$ml_error <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(
    list(
      input$min_deg_support,
      input$min_wgcna_support,
      input$same_direction,
      input$deg_weight
    ),
    {
      state$ml <- NULL
      state$ml_error <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(
    list(input$gnn_mode, input$gnn_candidates),
    {
      state$ml <- NULL
      state$ml_error <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(input$run_deg, {
    state$ml <- NULL
    state$ml_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_wgcna, {
    state$ml <- NULL
    state$ml_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_cross, {
    state$ml <- NULL
    state$ml_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_gnn, {
    state$ml <- NULL
    state$ml_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$ml_mode, {
    state$ml <- NULL
    state$ml_error <- NULL
  }, ignoreInit = TRUE)
  observeEvent(
    list(
      input$counts_files,
      input$metadata_files,
      input$fdr,
      input$logfc,
      input$top_genes,
      input$min_module,
      input$module_p,
      input$min_deg_support,
      input$min_wgcna_support,
      input$same_direction,
      input$deg_weight,
      input$gnn_mode,
      input$gnn_candidates,
      input$ml_mode
    ),
    {
      state$step8 <- NULL
      state$step8_error <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(input$run_deg, {
    state$step8 <- NULL
    state$step8_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_wgcna, {
    state$step8 <- NULL
    state$step8_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_cross, {
    state$step8 <- NULL
    state$step8_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_gnn, {
    state$step8 <- NULL
    state$step8_error <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$run_ml, {
    state$step8 <- NULL
    state$step8_error <- NULL
  }, ignoreInit = TRUE)
  observeEvent(ready_pairs(), {
    x <- ready_pairs()
    selected <- isolate(input$dataset_id)
    if (!length(selected) || !selected %in% x$key) {
      selected <- if (nrow(x)) x$key[1L] else character(0)
    }
    choices <- stats::setNames(x$key, x$Dataset)
    updateSelectInput(session, "dataset_id", choices = choices, selected = selected)
  }, ignoreInit = FALSE)

  active_pair <- reactive({
    req(input$dataset_id)
    x <- ready_pairs()
    row <- x[x$key == input$dataset_id, , drop = FALSE]
    req(nrow(row) == 1L)
    row
  })
  active_files <- reactive({
    tdrx_read_pair(active_pair())
  })
  current_deg <- reactive({
    req(input$dataset_id)
    x <- state$deg[[input$dataset_id]]
    req(!is.null(x))
    x
  })
  current_wgcna <- reactive({
    req(input$dataset_id)
    x <- state$wgcna[[input$dataset_id]]
    req(!is.null(x))
    x
  })

  output$overview_cards <- renderUI({
    x <- ready_pairs()
    completed_deg <- sum(x$key %in% names(state$deg))
    completed_wgcna <- sum(x$key %in% names(state$wgcna))
    counts <- state$checks
    samples <- sum(vapply(counts, function(a) if (!is.null(a$samples)) a$samples else 0L,
                          numeric(1)), na.rm = TRUE)
    tags$div(class = "metrics-grid",
      metric_tile("Paired datasets", nrow(x), "Ready for validation", "▣"),
      metric_tile("Validated samples", if (length(counts)) samples else "—",
                  "Across checked datasets", "◎"),
      metric_tile("DEG completed", completed_deg, "Separate gene-level results", "↗"),
      metric_tile("WGCNA completed", completed_wgcna, "Network modules & traits", "◉")
    )
  })

  output$paired_table <- DT::renderDT({
    x <- pairs()
    if (!nrow(x)) {
      return(DT::datatable(data.frame(Message = "Select counts CSVs and metadata CSVs in the left panel."),
                           options = list(dom = "t"), rownames = FALSE))
    }
    x$Genes <- "—"
    x$Samples <- "—"
    x$Case <- "—"
    x$Control <- "—"
    x$Note <- ""
    for (i in seq_len(nrow(x))) {
      c <- state$checks[[x$key[i]]]
      if (!is.null(c)) {
        x$Status[i] <- c$status
        x$Genes[i] <- if (!is.null(c$genes)) as.character(c$genes) else "—"
        x$Samples[i] <- if (!is.null(c$samples)) as.character(c$samples) else "—"
        x$Case[i] <- if (!is.null(c$case)) as.character(c$case) else "—"
        x$Control[i] <- if (!is.null(c$control)) as.character(c$control) else "—"
        x$Note[i] <- if (!is.null(c$note)) c$note else ""
      }
    }
    shown <- x[, c("Dataset", "Counts", "Metadata", "Status", "Genes", "Samples",
                    "Case", "Control", "Note"), drop = FALSE]
    DT::datatable(shown, options = list(pageLength = 8, scrollX = TRUE, autoWidth = TRUE),
                  rownames = FALSE)
  }, server = FALSE)

  output$counts_preview <- DT::renderDT({
    x <- active_files()$counts
    cols <- seq_len(min(7L, ncol(x)))
    DT::datatable(utils::head(x[, cols, drop = FALSE], 8L),
                  options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  }, server = FALSE)
  output$metadata_preview <- DT::renderDT({
    x <- active_files()$metadata
    DT::datatable(utils::head(x, 10L),
                  options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  }, server = FALSE)

  observeEvent(input$validate_pairs, {
    x <- ready_pairs()
    if (!nrow(x)) {
      showNotification("Upload at least one matching counts + metadata pair.", type = "warning")
      return()
    }
    results <- list()
    withProgress(message = "Validating all uploaded datasets", value = 0, {
      for (i in seq_len(nrow(x))) {
        setProgress(value = (i - 1) / nrow(x), detail = x$Dataset[i])
        key <- x$key[i]
        results[[key]] <- tryCatch({
          csv <- tdrx_read_pair(x[i, , drop = FALSE])
          p <- tdrx_prepare_input(csv$counts, csv$metadata)
          list(status = "Validated", genes = as.numeric(length(p$gene_ids)),
               samples = as.numeric(ncol(p$counts)),
               case = as.numeric(p$group_counts[["Case"]]),
               control = as.numeric(p$group_counts[["Control"]]),
               note = "Counts and metadata match")
        }, error = function(e) {
          list(status = "Validation error", note = conditionMessage(e))
        })
        setProgress(value = i / nrow(x))
      }
    })
    state$checks <- results
    good <- sum(vapply(results, function(x) identical(x$status, "Validated"), logical(1)))
    showNotification(paste0(good, "/", nrow(x), " dataset pairs validated."),
                     type = if (good == nrow(x)) "message" else "warning")
  })

  observeEvent(input$run_deg, {
    state$cross <- NULL
    x <- ready_pairs()
    if (!nrow(x)) {
      showNotification("Upload matching counts and metadata CSVs first.", type = "error")
      return()
    }
    state$deg <- list()
    state$wgcna <- list()
    completed <- list()
    errors <- list()
    withProgress(message = "Running batch DEG analysis", value = 0, {
      for (i in seq_len(nrow(x))) {
        key <- x$key[i]
        setProgress(value = (i - 1) / nrow(x), detail = paste0("DEG · ", x$Dataset[i]))
        result <- tryCatch({
          data <- tdrx_read_pair(x[i, , drop = FALSE])
          p <- tdrx_prepare_input(data$counts, data$metadata)
          tdrx_deg_analysis(p, fdr_cutoff = input$fdr,
                            log2fc_cutoff = input$logfc)
        }, error = function(e) {
          errors[[key]] <<- conditionMessage(e)
          NULL
        })
        if (!is.null(result)) completed[[key]] <- result
        setProgress(value = i / nrow(x))
      }
    })
    state$deg <- completed
    state$deg_errors <- errors
    state$wgcna_errors <- list()
    updateTabsetPanel(session, "result_tabs", selected = "deg")
    showNotification(paste0("DEG completed: ", length(completed), "/", nrow(x),
                            " datasets. See Datasets > Batch execution status."),
                     type = if (length(completed) == nrow(x)) "message" else "warning",
                     duration = 10)
  })

  observeEvent(input$run_wgcna, {
    state$cross <- NULL
    x <- ready_pairs()
    available <- x[x$key %in% names(state$deg), , drop = FALSE]
    if (!nrow(available)) {
      showNotification("Run DEG successfully before starting WGCNA.", type = "warning")
      return()
    }
    state$wgcna <- list()
    completed <- list()
    errors <- list()
    withProgress(message = "Inferring co-expression networks", value = 0, {
      for (i in seq_len(nrow(available))) {
        key <- available$key[i]
        setProgress(value = (i - 1) / nrow(available),
                    detail = paste0("WGCNA · ", available$Dataset[i]))
        result <- tryCatch({
          tdrx_wgcna_analysis(state$deg[[key]],
                              top_var_genes = input$top_genes,
                              min_module_size = input$min_module,
                              module_p_cutoff = input$module_p)
        }, error = function(e) {
          errors[[key]] <<- conditionMessage(e)
          NULL
        })
        if (!is.null(result)) completed[[key]] <- result
        setProgress(value = i / nrow(available))
      }
    })
    state$wgcna <- completed
    state$wgcna_errors <- errors
    updateTabsetPanel(session, "result_tabs", selected = "wgcna")
    showNotification(paste0("WGCNA completed: ", length(completed), "/", nrow(available),
                            " DEG datasets. Errors, if any, appear under Datasets."),
                     type = if (length(completed) == nrow(available)) "message" else "warning",
                     duration = 10)
  })

  observeEvent(input$run_cross, {
    if (length(state$deg) < 2L) {
      showNotification("Complete DEG on at least 2 datasets first.", type = "warning")
      return()
    }
    pairs_now <- ready_pairs()
    labels <- stats::setNames(pairs_now$Dataset, pairs_now$key)
    result <- tryCatch(
      tdrx_cross_analysis(
        deg_results = state$deg,
        wgcna_results = state$wgcna,
        dataset_names = labels,
        min_deg_support = input$min_deg_support,
        min_wgcna_support = input$min_wgcna_support,
        require_same_direction = input$same_direction,
        deg_weight = input$deg_weight
      ),
      error = function(e) {
        showNotification(paste("Step 5:", conditionMessage(e)),
                         type = "error", duration = 12)
        NULL
      }
    )
    state$cross <- result
    if (!is.null(result)) {
      updateTabsetPanel(session, "result_tabs", selected = "cross")
      showNotification(
        paste("Step 5 complete:", nrow(result$ranked), "ranked candidates."),
        type = "message", duration = 8
      )
    }
  })

  output$cross_status <- renderUI({
    x <- state$cross
    if (is.null(x)) {
      return(tags$div(class = "callout",
        "Complete DEG on 2+ datasets, then click Run Step 5 in the sidebar."))
    }
    msg <- if (x$wgcna_count == 0L) {
      "DEG-only comparison: no completed WGCNA analysis yet. DEG weight automatically becomes 100%."
    } else if (x$wgcna_count < x$dataset_count) {
      "Partial WGCNA coverage: only completed WGCNA datasets contribute WGCNA evidence."
    } else {
      "All completed DEG datasets have corresponding WGCNA results."
    }
    if (!nrow(x$ranked)) {
      msg <- paste(msg, "No genes pass current support/direction filters. Relax settings and rerun.")
    }
    tags$div(class = "callout", msg,
      tags$div(class = "cross-tip",
        "Exploratory rankings only. Case-Control contrasts may differ biologically among datasets."))
  })

  output$cross_metrics <- renderUI({
    x <- state$cross
    req(!is.null(x))
    tags$div(class = "mini-metrics",
      tags$div(tags$span("DEG datasets"), tags$strong(x$dataset_count)),
      tags$div(tags$span("Shared in all DEG datasets"), tags$strong(x$shared_all)),
      tags$div(tags$span("WGCNA datasets"), tags$strong(x$wgcna_count)),
      tags$div(tags$span("Ranked candidates"), tags$strong(nrow(x$ranked)))
    )
  })

  output$cross_table <- DT::renderDT({
    x <- state$cross
    req(!is.null(x))
    fields <- c("Rank", "Gene", "DEG_Support", "WGCNA_Support",
                "Same_direction", "S_DEG", "S_WGCNA", "Combined_Score")
    fields <- c(fields, grep("^(logFC|FDR)_", names(x$ranked), value = TRUE))
    if (!nrow(x$ranked)) {
      return(DT::datatable(data.frame(Message = "No candidate genes pass these settings."),
                           options = list(dom = "t"), rownames = FALSE))
    }
    DT::datatable(utils::head(x$ranked[, fields, drop = FALSE], 2000L),
                  options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  }, server = FALSE)

  output$cross_pairwise <- DT::renderDT({
    x <- state$cross
    req(!is.null(x))
    DT::datatable(x$pairwise, options = list(dom = "t", scrollX = TRUE),
                  rownames = FALSE)
  }, server = FALSE)

  output$cross_support_plot <- renderPlot({
    x <- state$cross
    req(!is.null(x))
    print(x$support_plot)
  }, res = 110)

  output$cross_heatmap <- renderPlot({
    x <- state$cross
    req(!is.null(x))
    if (is.null(x$heatmap_plot)) {
      plot.new()
      text(0.5, 0.5, "No genes pass the selected filters")
    } else {
      print(x$heatmap_plot)
    }
  }, res = 110)

  observeEvent(input$run_gnn, {

    state$gnn <- NULL
    state$gnn_error <- NULL

    if (is.null(state$cross)) {
      showNotification(
        "Run Step 5 · Cross-Dataset first.",
        type = "warning"
      )
      return()
    }

    if (length(state$deg) < 2L) {
      showNotification(
        "GNN needs DEG results from at least two datasets.",
        type = "warning"
      )
      return()
    }

    if (!requireNamespace("torch", quietly = TRUE)) {
      showNotification(
        "R torch is not installed. Install it, restart R, then try again.",
        type = "error",
        duration = 12
      )
      return()
    }

    available <- ready_pairs()

    labels <- stats::setNames(
      available$Dataset,
      available$key
    )

    result <- tryCatch({

      withProgress(
        message = "Running TD-CCHI GNN across dataset pairs",
        value = 0,
        {

          tdrx_gnn_batch_analysis(
            step5_result = state$cross,
            deg_results = state$deg,
            dataset_names = labels,
            candidate_limit = input$gnn_candidates,
            mode = input$gnn_mode,

            output_root = file.path(
              "TriDiseaseRank_results",
              "Step6_TD_CCHI"
            ),

            progress = function(value, detail) {
              setProgress(
                value = value,
                detail = detail
              )
            }
          )
        }
      )

    }, error = function(e) {

      state$gnn_error <- conditionMessage(e)

      showNotification(
        paste("GNN failed:", conditionMessage(e)),
        type = "error",
        duration = 15
      )

      NULL
    })

    state$gnn <- result

    if (!is.null(result)) {

      updateTabsetPanel(
        session,
        "result_tabs",
        selected = "gnn"
      )

      showNotification(
        paste0(
          "GNN completed: ",
          result$successful_pairs, "/",
          result$pair_count,
          " pair runs; ",
          length(result$candidate_genes),
          " common candidate genes."
        ),
        type = if (length(result$errors)) "warning" else "message",
        duration = 10
      )
    }
  })

  # ------------------------------------------
  # 7. GNN status
  # ------------------------------------------

  output$gnn_status <- renderUI({

    x <- state$gnn

    if (is.null(x)) {

      msg <- if (!is.null(state$gnn_error)) {
        paste("GNN failed:", state$gnn_error)
      } else {
        "Complete DEG and Step 5 first. Choose Quick Preview or Full Research, then click Run GNN."
      }

      return(tags$div(class = "callout", msg))
    }

    msg <- paste0(
      "Completed ", x$successful_pairs, " of ",
      x$pair_count, " dataset pairs using ",
      length(x$candidate_genes), " genes. Mode: ",
      if (x$mode == "quick") "Quick preview" else "Full research",
      "."
    )

    if (length(x$errors)) {
      msg <- paste0(
        msg,
        " Failed pairs: ",
        paste(names(x$errors), collapse = ", "),
        "."
      )
    }

    tags$div(
      class = "callout",
      msg,
      tags$div(
        class = "cross-tip",
        "This is an expression-derived gene network, not a curated protein-protein interaction network."
      )
    )
  })

  # ------------------------------------------
  # 8. Summary cards
  # ------------------------------------------

  output$gnn_metrics <- renderUI({

    x <- state$gnn
    req(!is.null(x))

    top_gene <- if (nrow(x$ranked)) {
      x$ranked$Gene[1L]
    } else {
      "—"
    }

    tags$div(
      class = "mini-metrics",

      tags$div(
        tags$span("Candidate genes"),
        tags$strong(length(x$candidate_genes))
      ),

      tags$div(
        tags$span("Successful pair runs"),
        tags$strong(paste0(
          x$successful_pairs, "/", x$pair_count
        ))
      ),

      tags$div(
        tags$span("Top network gene"),
        tags$strong(top_gene)
      ),

      tags$div(
        tags$span("Ranked genes"),
        tags$strong(nrow(x$ranked))
      )
    )
  })

  # ------------------------------------------
  # 9. Network ranking table
  # ------------------------------------------

  output$gnn_table <- DT::renderDT({

    x <- state$gnn
    req(!is.null(x))

    fields <- c(
      "Network_Rank",
      "Gene",
      "Network_Combined_Score",
      "S_GNN",
      "GNN_Rank",
      "GNN_PairTopSupport",
      "GNN_PairCount",
      "S_DEG",
      "S_WGCNA",
      "GNN_C",
      "GNN_I",
      "GNN_S",
      "GNN_R"
    )

    fields <- intersect(fields, names(x$ranked))

    DT::datatable(
      utils::head(
        x$ranked[, fields, drop = FALSE],
        1000L
      ),
      options = list(
        pageLength = 12,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  # ------------------------------------------
  # 10. Pairwise results
  # ------------------------------------------

  output$gnn_pairwise_table <- DT::renderDT({

    x <- state$gnn
    req(!is.null(x))

    fields <- intersect(
      c(
        "Pair",
        "Gene",
        "Rank",
        "CCHS",
        "Centrality_C",
        "InterCluster_I",
        "Stability_S",
        "Reproducibility_R",
        "Pair_Top"
      ),
      names(x$pairwise_long)
    )

    DT::datatable(
      utils::head(
        x$pairwise_long[, fields, drop = FALSE],
        1000L
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  # ------------------------------------------
  # 11. Graphs
  # ------------------------------------------

  output$gnn_score_plot <- renderPlot({

    x <- state$gnn
    req(!is.null(x))

    print(x$score_plot)

  }, res = 110)

  output$gnn_component_plot <- renderPlot({

    x <- state$gnn
    req(!is.null(x))

    if (is.null(x$component_plot)) {
      plot.new()
      text(0.5, 0.5, "No GNN component scores available")
    } else {
      print(x$component_plot)
    }

  }, res = 110)
  observeEvent(input$run_ml, {

    state$ml <- NULL
    state$ml_error <- NULL

    if (
      is.null(state$gnn) ||
      is.null(state$gnn$aggregate)
    ) {
      showNotification(
        "Complete Step 6 GNN analysis first.",
        type = "warning"
      )
      return()
    }

    if (length(state$deg) < 2L) {
      showNotification(
        "At least two completed DEG datasets are required.",
        type = "warning"
      )
      return()
    }

    available <- ready_pairs()

    labels <- stats::setNames(
      available$Dataset,
      available$key
    )

    current_mode <- input$ml_mode

    if (is.null(current_mode)) {
      current_mode <- "quick"
    }

    error_message <- NULL

    result <- tryCatch({

      withProgress(
        message = if (current_mode == "full") {
          "Running all 113 ML pipelines"
        } else {
          "Testing 8 representative ML pipelines"
        },
        value = 0,
        {

          tdrx_ml_batch_analysis(
            gnn_result = state$gnn,
            deg_results = state$deg,
            dataset_names = labels,
            mode = current_mode,

            output_root = file.path(
              "TriDiseaseRank_results",
              "Step7_113ML"
            ),

            progress = function(value, detail) {
              setProgress(
                value = value,
                detail = detail
              )
            }
          )
        }
      )

    }, error = function(e) {

      error_message <<- conditionMessage(e)

      showNotification(
        paste("ML Analysis failed:", conditionMessage(e)),
        type = "error",
        duration = 15
      )

      NULL
    })

    state$ml_error <- error_message
    state$ml <- result

    if (!is.null(result)) {

      updateTabsetPanel(
        session,
        "result_tabs",
        selected = "ml"
      )

      showNotification(
        paste0(
          "ML completed: ",
          result$successful_datasets,
          "/",
          result$requested_datasets,
          " datasets. ",
          result$models_requested,
          " pipelines requested per dataset."
        ),
        type = if (
          length(result$dataset_errors)
        ) "warning" else "message",
        duration = 12
      )
    }
  })

  # -------------------------------------------------------
  # 7. ML status
  # -------------------------------------------------------

  output$ml_status <- renderUI({

    x <- state$ml

    if (is.null(x)) {

      msg <- if (!is.null(state$ml_error)) {
        paste("ML failed:", state$ml_error)
      } else {
        "Complete DEG, Step 5 and Step 6 first. Then select an ML mode and click Run ML."
      }

      return(
        tags$div(class = "callout", msg)
      )
    }

    msg <- paste0(
      x$models_requested,
      " model pipelines requested per dataset. ",
      x$successful_datasets,
      " of ",
      x$requested_datasets,
      " datasets completed."
    )

    if (length(x$dataset_errors)) {
      msg <- paste0(
        msg,
        " Failed datasets: ",
        paste(names(x$dataset_errors), collapse = ", "),
        "."
      )
    }

    tags$div(
      class = "callout",
      msg,
      tags$div(
        class = "cross-tip",
        "Internal validation only; upstream GNN selection is not nested inside these CV folds."
      )
    )
  })

  # -------------------------------------------------------
  # 8. ML summary cards
  # -------------------------------------------------------

  output$ml_metrics <- renderUI({

    x <- state$ml
    req(!is.null(x))

    successful_models <- sum(
      x$model_performance$Status == "OK",
      na.rm = TRUE
    )

    top_gene <- if (nrow(x$gene_scores)) {
      x$gene_scores$Gene[1L]
    } else {
      "—"
    }

    tags$div(
      class = "mini-metrics",

      tags$div(
        tags$span("GNN candidate genes"),
        tags$strong(length(x$candidate_genes))
      ),

      tags$div(
        tags$span("Completed datasets"),
        tags$strong(paste0(
          x$successful_datasets,
          "/",
          x$requested_datasets
        ))
      ),

      tags$div(
        tags$span("Successful model runs"),
        tags$strong(successful_models)
      ),

      tags$div(
        tags$span("Top ML-supported gene"),
        tags$strong(top_gene)
      )
    )
  })

  # -------------------------------------------------------
  # 9. Best model table
  # -------------------------------------------------------

  output$ml_best_table <- DT::renderDT({

    x <- state$ml
    req(!is.null(x))

    fields <- intersect(
      c(
        "Dataset",
        "Model_ID",
        "Model",
        "Selector",
        "Classifier",
        "Internal_AUC",
        "Selected_Genes_n",
        "Outer_Folds"
      ),
      names(x$best_models)
    )

    DT::datatable(
      x$best_models[, fields, drop = FALSE],
      options = list(
        dom = "t",
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  # -------------------------------------------------------
  # 10. All 113 model performance table
  # -------------------------------------------------------

  output$ml_model_table <- DT::renderDT({

    x <- state$ml
    req(!is.null(x))

    z <- x$model_performance

    z <- z[
      order(
        z$Dataset,
        z$Status != "OK",
        -z$Internal_AUC,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    DT::datatable(
      utils::head(z, 500L),
      options = list(
        pageLength = 15,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  # -------------------------------------------------------
  # 11. Combined gene ranking table
  # -------------------------------------------------------

  output$ml_gene_table <- DT::renderDT({

    x <- state$ml
    req(!is.null(x))

    DT::datatable(
      x$gene_scores,
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  # -------------------------------------------------------
  # 12. Model performance plot
  # -------------------------------------------------------

  output$ml_auc_plot <- renderPlot({

    x <- state$ml
    req(!is.null(x))

    if (is.null(x$model_plot)) {
      plot.new()
      text(0.5, 0.5, "No successful model AUC available")
    } else {
      print(x$model_plot)
    }

  }, res = 110)

  # -------------------------------------------------------
  # 13. ML gene score plot
  # -------------------------------------------------------

  output$ml_gene_plot <- renderPlot({

    x <- state$ml
    req(!is.null(x))

    if (is.null(x$gene_plot)) {
      plot.new()
      text(0.5, 0.5, "No ML gene scores available")
    } else {
      print(x$gene_plot)
    }

  }, res = 110)
  observeEvent(input$run_step8, {

    state$step8 <- NULL
    state$step8_error <- NULL

    if (is.null(state$cross)) {
      showNotification(
        "Complete Step 5 first.",
        type = "warning"
      )
      return()
    }

    if (is.null(state$gnn)) {
      showNotification(
        "Complete Step 6 GNN first.",
        type = "warning"
      )
      return()
    }

    if (is.null(state$ml)) {
      showNotification(
        "Complete Step 7 ML first.",
        type = "warning"
      )
      return()
    }

    available <- ready_pairs()

    labels <- stats::setNames(
      available$Dataset,
      available$key
    )

    result <- tryCatch({

      withProgress(
        message = "Running Pathway and Final Ranking",
        value = 0,
        {

          tdrx_step8_analysis(

            deg_results = state$deg,

            wgcna_results = state$wgcna,

            cross_result = state$cross,

            gnn_result = state$gnn,

            ml_result = state$ml,

            dataset_names = labels,

            top_n = input$step8_top_n,

            output_root = file.path(
              "TriDiseaseRank_results",
              "Step8_Final_Ranking"
            ),

            progress = function(value, detail) {
              setProgress(
                value = value,
                detail = detail
              )
            }
          )
        }
      )

    }, error = function(e) {

      state$step8_error <- conditionMessage(e)

      showNotification(
        paste(
          "Step 8 failed:",
          conditionMessage(e)
        ),
        type = "error",
        duration = 18
      )

      NULL
    })

    state$step8 <- result

    if (!is.null(result)) {

      updateTabsetPanel(
        session,
        "result_tabs",
        selected = "step8"
      )

      showNotification(
        paste(
          "Step 8 completed.",
          nrow(result$final_ranking),
          "genes ranked.",
          result$ranking_status
        ),
        type = "message",
        duration = 12
      )
    }
  })

  # ----------------------------------------------------------
  # Step 8 status
  # ----------------------------------------------------------

  output$step8_status <- renderUI({

    x <- state$step8

    if (is.null(x)) {

      msg <- if (!is.null(state$step8_error)) {
        paste("Step 8 failed:", state$step8_error)
      } else {
        "Complete DEG, WGCNA, Step 5, GNN and ML first."
      }

      return(
        tags$div(class = "callout", msg)
      )
    }

    msg <- if (
      x$ranking_status == "PRELIMINARY_QUICK_ML"
    ) {

      "PRELIMINARY: Quick ML was used. Run Full Research · 113 models and rerun Step 8."

    } else if (
      x$ranking_status == "FULL_113ML_INTERNAL_VALIDATION"
    ) {

      "Full 113-model internal validation completed. Independent external validation is still required."

    } else {

      "PRELIMINARY: Some model runs or evidence layers have incomplete coverage."
    }

    databases <- x$pathway_result$databases

    if (is.null(databases) || !length(databases)) {
      databases <- "No significant enrichment terms returned"
    } else {
      databases <- paste(databases, collapse = ", ")
    }

    tags$div(
      class = "callout",
      msg,
      tags$div(
        class = "cross-tip",
        paste("Pathway databases:", databases)
      )
    )
  })

  # ----------------------------------------------------------
  # Step 8 summary cards
  # ----------------------------------------------------------

  output$step8_metrics <- renderUI({

    x <- state$step8
    req(!is.null(x))

    pathway_rows <- nrow(
      x$pathway_result$enrichment
    )

    top_gene <- if (nrow(x$final_ranking)) {
      x$final_ranking$Gene[1L]
    } else {
      "—"
    }

    tags$div(
      class = "mini-metrics",

      tags$div(
        tags$span("Final candidates"),
        tags$strong(nrow(x$final_ranking))
      ),

      tags$div(
        tags$span("Top biomarker"),
        tags$strong(top_gene)
      ),

      tags$div(
        tags$span("Enrichment rows"),
        tags$strong(pathway_rows)
      ),

      tags$div(
        tags$span("Pathway datasets"),
        tags$strong(nrow(x$pathway_input_counts))
      )
    )
  })

  # ----------------------------------------------------------
  # Step 8 tables
  # ----------------------------------------------------------

  output$step8_top_table <- DT::renderDT({

    x <- state$step8
    req(!is.null(x))

    DT::datatable(
      x$top_biomarkers,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  output$step8_final_table <- DT::renderDT({

    x <- state$step8
    req(!is.null(x))

    DT::datatable(
      utils::head(
        x$final_ranking,
        1000L
      ),
      options = list(
        pageLength = 15,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  output$step8_pathway_table <- DT::renderDT({

    x <- state$step8
    req(!is.null(x))

    z <- x$pathway_result$enrichment

    if (!is.data.frame(z) || !nrow(z)) {

      return(
        DT::datatable(
          data.frame(
            Message = "No significant pathway terms returned."
          ),
          options = list(dom = "t"),
          rownames = FALSE
        )
      )
    }

    if ("p.adjust" %in% names(z)) {
      z <- z[
        order(z$p.adjust, na.last = TRUE),
        ,
        drop = FALSE
      ]
    }

    DT::datatable(
      utils::head(z, 1000L),
      options = list(
        pageLength = 12,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  output$step8_convergence_table <- DT::renderDT({

    x <- state$step8
    req(!is.null(x))

    z <- x$pathway_result$convergence

    if (!is.data.frame(z) || !nrow(z)) {

      return(
        DT::datatable(
          data.frame(
            Message = "No significant pathway convergence."
          ),
          options = list(dom = "t"),
          rownames = FALSE
        )
      )
    }

    DT::datatable(
      utils::head(z, 1000L),
      options = list(
        pageLength = 12,
        scrollX = TRUE
      ),
      rownames = FALSE
    )

  }, server = FALSE)

  # ----------------------------------------------------------
  # Step 8 plots
  # ----------------------------------------------------------

  output$step8_final_plot <- renderPlot({

    x <- state$step8
    req(!is.null(x))

    print(x$final_score_plot)

  }, res = 110)

  output$step8_evidence_plot <- renderPlot({

    x <- state$step8
    req(!is.null(x))

    print(x$evidence_plot)

  }, res = 110)

  output$step8_pathway_plot <- renderPlot({

    x <- state$step8
    req(!is.null(x))

    if (is.null(x$pathway_plot)) {

      plot.new()

      text(
        0.5, 0.5,
        "No significant pathway terms available"
      )

    } else {

      print(x$pathway_plot)
    }

  }, res = 110)
  output$execution_table <- DT::renderDT({
    x <- ready_pairs()
    if (!nrow(x)) return(DT::datatable(data.frame(Message = "No valid pairs yet."),
                                           options = list(dom = "t"), rownames = FALSE))
    rows <- lapply(seq_len(nrow(x)), function(i) {
      k <- x$key[i]
      deg_error <- state$deg_errors[[k]]
      wg_error <- state$wgcna_errors[[k]]
      data.frame(
        Dataset = x$Dataset[i],
        DEG = if (k %in% names(state$deg)) "Complete" else
          if (!is.null(deg_error)) "Failed" else "Not run",
        WGCNA = if (k %in% names(state$wgcna)) "Complete" else
          if (!is.null(wg_error)) "Failed" else "Not run",
        Error_details = paste(c(deg_error, wg_error), collapse = " | "),
        stringsAsFactors = FALSE
      )
    })
    DT::datatable(do.call(rbind, rows), options = list(dom = "t", scrollX = TRUE),
                  rownames = FALSE)
  }, server = FALSE)

  output$deg_status <- renderUI({
    key <- input$dataset_id
    if (is.null(key) || !key %in% names(state$deg)) {
      msg <- if (!is.null(key) && !is.null(state$deg_errors[[key]])) {
        paste0("DEG failed: ", state$deg_errors[[key]])
      } else "No DEG result for the selected dataset yet. Click 'Run DEG · all datasets'."
      tags$div(class = "callout", msg)
    }
  })
  output$deg_metrics <- renderUI({
    x <- current_deg()
    tags$div(class = "mini-metrics",
      tags$div(tags$span("Samples"), tags$strong(x$analyzed_samples)),
      tags$div(tags$span("Tested genes"), tags$strong(format(x$tested_genes, big.mark = ","))),
      tags$div(tags$span("Significant DEGs"), tags$strong(format(x$significant_genes, big.mark = ",")))
    )
  })
  output$deg_table <- DT::renderDT({
    x <- current_deg()$all
    columns <- intersect(c("gene_id", "gene_symbol", "logFC", "P.Value", "adj.P.Val",
                           "S_DEG", "Significant"), names(x))
    DT::datatable(utils::head(x[, columns, drop = FALSE], 2000L), options = list(pageLength = 12,
                  scrollX = TRUE), rownames = FALSE)
  }, server = FALSE)
  output$volcano_plot <- renderPlot({
    print(current_deg()$volcano)
  }, res = 100)
  output$pca_plot <- renderPlot({
    plot <- current_deg()$pca
    if (is.null(plot)) {
      plot.new()
      text(0.5, 0.5, "Not enough sample variation for PCA")
    } else print(plot)
  }, res = 100)

  output$wgcna_status <- renderUI({
    key <- input$dataset_id
    if (is.null(key) || !key %in% names(state$wgcna)) {
      msg <- if (!is.null(key) && !is.null(state$wgcna_errors[[key]])) {
        paste0("WGCNA failed: ", state$wgcna_errors[[key]])
      } else "WGCNA has not been run for this dataset. Complete DEG, then run WGCNA."
      tags$div(class = "callout", msg)
    }
  })
  output$wgcna_metrics <- renderUI({
    x <- current_wgcna()
    tags$div(class = "mini-metrics",
      tags$div(tags$span("Soft power"), tags$strong(x$soft_power)),
      tags$div(tags$span("Detected modules"), tags$strong(nrow(x$modules))),
      tags$div(tags$span("Nominally significant modules"),
               tags$strong(x$significant_modules)),
      tags$div(tags$span("Module-supported genes"),
               tags$strong(nrow(x$significant_genes)))
    )
  })
  output$power_plot <- renderPlot({
    print(current_wgcna()$power_plot)
  }, res = 100)
  output$trait_plot <- renderPlot({
    plot <- current_wgcna()$trait_plot
    if (is.null(plot)) {
      plot.new()
      text(0.5, 0.5, "No non-grey co-expression modules detected")
    } else print(plot)
  }, res = 100)
  output$modules_table <- DT::renderDT({
    DT::datatable(current_wgcna()$modules, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  }, server = FALSE)
  output$wgcna_gene_table <- DT::renderDT({
    DT::datatable(utils::head(current_wgcna()$all_genes, 2000L), rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  }, server = FALSE)

  output$download_cross <- downloadHandler(
    filename = function() "Step5_DEG_WGCNA_ranked_candidates.csv",
    content = function(file) {
      req(!is.null(state$cross))
      utils::write.csv(state$cross$ranked, file, row.names = FALSE, na = "")
    }
  )
  output$download_cross_all <- downloadHandler(
    filename = function() "Step5_all_gene_evidence.csv",
    content = function(file) {
      req(!is.null(state$cross))
      utils::write.csv(state$cross$universe, file, row.names = FALSE, na = "")
    }
  )
  output$download_cross_pairs <- downloadHandler(
    filename = function() "Step5_pairwise_DEG_overlap.csv",
    content = function(file) {
      req(!is.null(state$cross))
      utils::write.csv(state$cross$pairwise, file, row.names = FALSE, na = "")
    }
  )

  output$download_gnn <- downloadHandler(
    filename = function() "Step6_TD_CCHI_network_ranked.csv",
    content = function(file) {
      req(!is.null(state$gnn))
      utils::write.csv(
        state$gnn$ranked,
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )

  output$download_gnn_aggregate <- downloadHandler(
    filename = function() "Step6_TD_CCHI_aggregate_scores.csv",
    content = function(file) {
      req(!is.null(state$gnn))
      utils::write.csv(
        state$gnn$aggregate,
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )

  output$download_gnn_pairs <- downloadHandler(
    filename = function() "Step6_TD_CCHI_pairwise_scores.csv",
    content = function(file) {
      req(!is.null(state$gnn))
      utils::write.csv(
        state$gnn$pairwise_long,
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )
  output$download_ml_performance <- downloadHandler(

    filename = function() {
      "Step7_113ML_model_performance.csv"
    },

    content = function(file) {
      req(!is.null(state$ml))

      utils::write.csv(
        state$ml$model_performance,
        file,
        row.names = FALSE
      )
    }
  )

  output$download_ml_genes <- downloadHandler(

    filename = function() {
      "Step7_113ML_aggregate_gene_scores.csv"
    },

    content = function(file) {
      req(!is.null(state$ml))

      utils::write.csv(
        state$ml$gene_scores,
        file,
        row.names = FALSE
      )
    }
  )

  output$download_ml_best <- downloadHandler(

    filename = function() {
      "Step7_113ML_best_models.csv"
    },

    content = function(file) {
      req(!is.null(state$ml))

      utils::write.csv(
        state$ml$best_models,
        file,
        row.names = FALSE
      )
    }
  )
  output$download_step8_final <- downloadHandler(

    filename = function() {
      "TriDiseaseRankX_five_layer_ranking.csv"
    },

    content = function(file) {
      req(!is.null(state$step8))

      file.copy(
        state$step8$ranking_file,
        file,
        overwrite = TRUE
      )
    }
  )

  output$download_step8_top <- downloadHandler(

    filename = function() {
      "TriDiseaseRankX_top_biomarkers.csv"
    },

    content = function(file) {
      req(!is.null(state$step8))

      file.copy(
        state$step8$top_file,
        file,
        overwrite = TRUE
      )
    }
  )

  output$download_step8_pathway <- downloadHandler(

    filename = function() {
      "TriDiseaseRankX_pathway_enrichment.csv"
    },

    content = function(file) {
      req(!is.null(state$step8))

      file.copy(
        state$step8$enrichment_file,
        file,
        overwrite = TRUE
      )
    }
  )

  output$download_step8_workbook <- downloadHandler(

    filename = function() {
      "TriDiseaseRankX_results.xlsx"
    },

    content = function(file) {
      req(!is.null(state$step8))

      file.copy(
        state$step8$workbook_file,
        file,
        overwrite = TRUE
      )
    }
  )
  output$download_deg <- downloadHandler(
    filename = function() paste0(tdrx_safe_name(input$dataset_id), "_DEG_all.csv"),
    content = function(file) {
      utils::write.csv(current_deg()$all, file, row.names = FALSE, na = "")
    }
  )
  output$download_deg_sig <- downloadHandler(
    filename = function() paste0(tdrx_safe_name(input$dataset_id), "_DEG_significant.csv"),
    content = function(file) {
      utils::write.csv(current_deg()$significant, file, row.names = FALSE, na = "")
    }
  )
  output$download_wgcna <- downloadHandler(
    filename = function() paste0(
      tdrx_safe_name(input$dataset_id),
      "_WGCNA_all_genes.csv"
    ),
    content = function(file) {
      utils::write.csv(current_wgcna()$all_genes, file, row.names = FALSE, na = "")
    }
  )
  output$download_zip <- downloadHandler(
    filename = function() paste0("TriDiseaseRankX_Batch_", Sys.Date(), ".zip"),
    content = function(file) {
      req(length(state$deg) + length(state$wgcna) > 0L)
      if (!requireNamespace("zip", quietly = TRUE)) {
        stop("Install the 'zip' package to download batch results.")
      }
      export_root <- tempfile("tdrx_results_")
      dir.create(export_root)
      on.exit(unlink(export_root, recursive = TRUE), add = TRUE)
      keys <- unique(c(names(state$deg), names(state$wgcna)))
      for (i in seq_along(keys)) {
        key <- keys[i]
        folder <- file.path(export_root, sprintf("%02d_%s", i, tdrx_safe_name(key)))
        dir.create(folder)
        de <- state$deg[[key]]
        if (!is.null(de)) {
          utils::write.csv(de$all, file.path(folder, "DEG_all.csv"), row.names = FALSE)
          utils::write.csv(de$significant, file.path(folder, "DEG_significant.csv"), row.names = FALSE)
          ggplot2::ggsave(file.path(folder, "Volcano.png"), plot = de$volcano,
                          width = 7, height = 5.5, dpi = 180)
          if (!is.null(de$pca)) {
            ggplot2::ggsave(file.path(folder, "PCA.png"), plot = de$pca,
                            width = 7, height = 5.5, dpi = 180)
          }
        }
        wg <- state$wgcna[[key]]
        if (!is.null(wg)) {
          utils::write.csv(wg$modules, file.path(folder, "module_trait.csv"), row.names = FALSE)
          utils::write.csv(wg$all_genes, file.path(folder, "WGCNA_all_genes.csv"), row.names = FALSE)
          utils::write.csv(wg$significant_genes,
                           file.path(folder, "WGCNA_significant_module_genes.csv"), row.names = FALSE)
          ggplot2::ggsave(file.path(folder, "Soft_Power.png"), plot = wg$power_plot,
                          width = 7, height = 5.5, dpi = 180)
          if (!is.null(wg$trait_plot)) {
            ggplot2::ggsave(file.path(folder, "Module_Trait.png"), plot = wg$trait_plot,
                            width = 7, height = 5.5, dpi = 180)
          }
        }
      }
      if (!is.null(state$cross)) {
        utils::write.csv(state$cross$ranked,
          file.path(export_root, "Step5_DEG_WGCNA_ranked.csv"), row.names = FALSE)
        utils::write.csv(state$cross$universe,
          file.path(export_root, "Step5_all_gene_evidence.csv"), row.names = FALSE)
        utils::write.csv(state$cross$pairwise,
          file.path(export_root, "Step5_pairwise_overlap.csv"), row.names = FALSE)
      }
      if (!is.null(state$gnn)) {

        utils::write.csv(
          state$gnn$ranked,
          file.path(export_root, "Step6_TD_CCHI_network_ranked.csv"),
          row.names = FALSE
        )

        utils::write.csv(
          state$gnn$aggregate,
          file.path(export_root, "Step6_TD_CCHI_aggregate.csv"),
          row.names = FALSE
        )

        utils::write.csv(
          state$gnn$pairwise_long,
          file.path(export_root, "Step6_TD_CCHI_pairwise_long.csv"),
          row.names = FALSE
        )
      }
      if (!is.null(state$ml)) {

        utils::write.csv(
          state$ml$model_performance,
          file.path(
            export_root,
            "Step7_113ML_model_performance.csv"
          ),
          row.names = FALSE
        )

        utils::write.csv(
          state$ml$gene_scores,
          file.path(
            export_root,
            "Step7_113ML_gene_scores.csv"
          ),
          row.names = FALSE
        )

        utils::write.csv(
          state$ml$best_models,
          file.path(
            export_root,
            "Step7_113ML_best_models.csv"
          ),
          row.names = FALSE
        )
      }
      manifest <- data.frame(Dataset = keys,
        DEG = keys %in% names(state$deg),
        WGCNA = keys %in% names(state$wgcna), stringsAsFactors = FALSE)
      utils::write.csv(manifest, file.path(export_root, "manifest.csv"), row.names = FALSE)
      zip::zipr(zipfile = file,
                files = list.files(export_root, recursive = TRUE, full.names = FALSE),
                root = export_root)
    }
  )
}

shinyApp(ui, server)
