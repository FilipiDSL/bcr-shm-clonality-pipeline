# BCR SHM Clonality Pipeline

📄 **Read the full thesis:** [Computational Modeling Of The Maturation And Clonality Of B Lymphocytes In Colon Adenocarcinoma](docs/Computational_Modeling_Of_The_Maturation_And_Clonality_Of_B_Lymphocytes_In_Colon_Adenocarcinoma.pdf)

## Overview
This repository contains an R pipeline for analyzing B-cell receptor (BCR) repertoires, focusing on somatic hypermutation (SHM), clonality, and isotype switching (CSR) in Colorectal Adenocarcinoma (COAD). The pipeline processes single-cell RNA-seq and V(D)J data to compare tumor and adjacent normal tissues, providing quantitative evidence compatible with mature Tertiary Lymphoid Structures (TLS).

## Pipeline Steps
The analysis is divided into sequential R scripts located in the `scripts/` directory:
- `00_prep.R`: Data preparation and initial filtering of tumor vs. adjacent cells.
- `01_pre_proc.R`: Quality control, normalization, PCA, and UMAP generation.
- `02_heterogeneidade.R`: Analysis of B-cell subtype heterogeneity and distribution.
- `03_clonalidade_do_repertorio_BCR.R`: Clonal expansion and diversity analysis (Shannon/Simpson).
- `04_SHM.R`: Somatic hypermutation analysis and V-gene identity comparisons.
- `05_isotipos.R`: Class switch recombination (CSR) and isotype frequency analysis.
- `06_fig_resumo.R`: Generation of a comprehensive summary mosaic figure.

## Installation
To ensure reproducibility, this project requires specific R packages.

1. Clone the repository:
   ```bash
   git clone https://github.com/FilipiDSL/bcr-shm-clonality-pipeline.git
   cd bcr-shm-clonality-pipeline
   ```

2. Open R or RStudio in the project directory and install the required packages:
   ```R
   install.packages(c("Seurat", "dplyr", "ggplot2", "tidyr", "scales", "patchwork", "vegan", "rstatix", "pheatmap", "yaml"))
   ```

## Quickstart
To run the entire pipeline sequentially, execute the runner script from the project root:

```bash
Rscript scripts/run_all.R
```
Or inside an R session:
```R
source("scripts/run_all.R")
```

## Outputs
All generated outputs are saved in the `results/` directory:
- **Figures**: `results/figures/` (UMAPs, boxplots, heatmaps, summary mosaics in PNG/PDF formats).
- **Tables**: `results/tables/` (Statistical test results, e.g., Wilcoxon tests for diversity).

## Data
Due to privacy and licensing restrictions, the raw `.rds` data files are **not** included in this repository. 

To run the pipeline, you must place the required data files in the `data/raw/` directory:
- `COAD_Tumor_Adjacent.rds` (or the preprocessed `BCR_COAD_preprocessed.rds`)

*See `data/README.md` for detailed instructions on data formatting and acquisition.*

## Configuration
Paths and parameters are centralized in `configs/config.yml`. You can modify this file to adjust input/output directories, file names, or analysis parameters (e.g., random seed, PCA dimensions) without changing the R scripts.

## Reproducibility
- **Environment**: Ensure you have the required packages installed as listed in the Installation section.
- **Seed**: A global seed is set via `configs/config.yml` (default: 123) to ensure reproducible UMAPs and statistical simulations.
- **Paths**: All scripts use relative paths based on the project root.

## Citation
If you use this pipeline in your research, please cite:
> Filipi, L. (2025). Modelagem Computacional de Linfócitos B no COAD. [Undergraduate Thesis].
*(See `CITATION.cff` for more details)*
