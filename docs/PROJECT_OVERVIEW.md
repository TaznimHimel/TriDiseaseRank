# Project overview

TriDiseaseRank-X is a local Shiny application that wraps the existing R research workflow.

## Analysis flow

1. Counts/metadata pairing and validation (`R/batch_files.R`)
2. Differential expression (`R/deg_analysis.R`)
3. Co-expression modules (`R/wgcna_analysis.R`)
4. Shared DEG and cross-dataset evidence (`R/cross_dataset.R`)
5. TD-CCHI expression-correlation graph analysis (`R/gnn_analysis.R` + `original/PPI_GNN_R.R`)
6. GNN-guided ML (`R/ml_analysis.R` + `original/TriDiseaseRankX_ML.R`)
7. Pathway and five-layer ranking (`R/step8_analysis.R` + `original/Pathway.R`)

## Local-first behavior

Uploaded files and generated results are processed locally by the Shiny process. The application does not need a hosted database. Internet is required for initial package installation and may be required by online pathway resources.

## Important interpretation notes

- The GNN graph is based on expression correlation in the current implementation, not a curated protein interaction graph.
- GEO studies have study-specific designs; do not assume their case/control labels are biologically interchangeable.
- A final score is a prioritization heuristic, not a probability or validated diagnostic classifier.
