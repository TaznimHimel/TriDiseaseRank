# TriDiseaseRank-X

**A local Shiny workspace for multi-dataset gene-expression analysis and exploratory biomarker prioritization.**

TriDiseaseRank-X brings the project's R research scripts into a browser-based interface. It supports multiple count-matrix/metadata CSV pairs and provides a guided workflow for differential expression, WGCNA, cross-dataset comparison, expression-derived graph analysis, machine-learning model comparison, pathway enrichment, and a five-layer evidence ranking.

> **Research-use notice:** Results are exploratory. The bundled GEO studies do not necessarily represent the same biological contrast. A ranked gene is not automatically a clinically validated biomarker. Review study design, metadata, model assumptions, and validation requirements before making scientific claims.

## Features

- Upload multiple counts CSV files and matching metadata CSV files together.
- Validate sample names and dataset pairing.
- Differential expression using edgeR/limma-voom.
- WGCNA co-expression module analysis.
- Cross-dataset significant-gene comparison and evidence ranking.
- TD-CCHI graph analysis based on gene-expression correlations.
- GNN-guided ML analysis using candidate hub genes and the original ML engine.
- GO/KEGG enrichment and optional pathway databases when dependencies are available.
- Five-layer evidence summary and CSV/Excel exports.

## Requirements

- Windows 10/11, 64-bit.
- R for Windows (the project was assembled with R 4.6.x; use the same R minor version where possible).
- Internet connection for initial package installation and online pathway resources.
- Recommended: 8 GB RAM minimum; 16 GB or more for larger WGCNA/GNN/full-ML runs.
- Rtools45 may be needed if a dependency must be built from source. RStudio and VS Code are optional.

## Quick start (Windows)

1. Download the repository using **Code → Download ZIP**, then extract it.
2. Install R for Windows: https://cran.r-project.org/bin/windows/base/
3. Open the extracted repository folder.
4. Double-click `setup.bat` to install required CRAN/Bioconductor packages and attempt Torch runtime setup.
5. Double-click `run_app.bat`.
6. The app opens in your browser. Keep the console window open while using it.

You do **not** need RStudio or VS Code to run the app.

### Manual setup

Run these commands from the repository root:

```powershell
Rscript scripts\\setup_windows.R
Rscript scripts\\check_environment.R
Rscript scripts\\launch_app.R
```

If Torch is installed but LibTorch is not, run this in R:

```r
torch::install_torch()
```

Restart R and run the environment check again.

## Input CSV format

### Counts matrix
- One row per gene.
- Include `gene_id` and `gene_symbol` columns.
- Remaining columns are sample IDs with non-negative raw integer counts.
- Sample IDs must match the `Samples` column in the paired metadata.

### Metadata
- Must contain `Samples` and `Condition`.
- `Condition` must use `Control` and `Case` for the current workflow.
- Each metadata sample must match a counts-matrix column.
- Dataset pairs are matched by filename prefix, e.g. `GSE120852_counts_data.csv` and `GSE120852_metadata.csv`.

The `data/example/` folder is intended for the three GEO example pairs. Upload compatible CSV pairs through the app.

## Suggested workflow

1. Upload and validate dataset pairs.
2. Run DEG.
3. Run WGCNA.
4. Run cross-dataset comparison.
5. Run GNN (quick preview first).
6. Run ML (quick test first; full 113-pipeline mode is computationally expensive).
7. Run pathway enrichment and inspect the five-layer ranking.
8. Download result tables and workbook.

## Project layout

```text
TriDiseaseRank-X/
├── app.R
├── R/                  # Web-integrated analysis modules
├── original/           # Original research R scripts
├── data/example/       # Example GEO CSV pairs
├── www/                # CSS and web assets
├── scripts/            # Setup, environment check, launch
├── setup.bat
├── run_app.bat
└── README.md
```

## Outputs

Generated files are written locally under `TriDiseaseRank_results/`. This directory should not be committed because it can contain large or user-specific analysis outputs.

## Troubleshooting

- **Rscript not found:** Install R and reopen the terminal; check `where Rscript`.
- **Package installation fails:** Read the exact package error. Some source builds require Rtools45. Re-run setup after fixing the dependency.
- **Torch/LibTorch error:** Run `torch::install_torch()`, restart R, then rerun the environment check.
- **App cannot find `R/` or `original/`:** Launch from the repository root using `run_app.bat`.
- **No significant pathway terms:** This can be a valid result. Optional Reactome, Disease Ontology, and WikiPathways resources may be unavailable.
- **Large analyses are slow:** Start with Quick Preview/Quick Test. Full GNN and 113-model ML can require substantial time and memory.

## Scientific and implementation notes

- The current GNN derives a gene-gene network from expression correlations; it does not automatically use a curated STRING, BioGRID, or IntAct protein-protein interaction database.
- The three GEO datasets have different study designs and biological contrasts. Their shared genes should not be interpreted as equivalent evidence for one common disease without study-specific review.
- Internal cross-validation does not replace independent external validation. Upstream feature selection before cross-validation can produce optimistic performance estimates.
- The five-layer score is an exploratory prioritization score, not a calibrated probability or clinical diagnostic score.

## Original research scripts

Original research scripts belong under `original/`; web-integrated modules belong under `R/`. Keep original scripts unchanged where possible.

## License

See [LICENSE](LICENSE). Third-party datasets and package dependencies remain subject to their respective terms and licenses.
