# Data Directory

This directory is intended to store the input data for the BCR SHM Clonality Pipeline.

## Important Note on Data Privacy
Due to patient privacy regulations and potential licensing restrictions, the raw and processed `.rds` files containing single-cell RNA-seq and V(D)J data are **not** included in this repository.

## Required Files
To run the pipeline, you must place the following files in the `data/raw/` directory:

1.  **`BCR_COAD_preprocessed.rds`** (Optional): The initial preprocessed Seurat object containing all cell types. If provided, `scripts/00_prep.R` will filter it for "Cancer" and "Adjacent" cells.
2.  **`COAD_Tumor_Adjacent.rds`**: The filtered Seurat object containing only "Cancer" and "Adjacent" cells. If you already have this file, you can skip `00_prep.R` and start from `01_pre_proc.R`.

## Data Format
The input data should be a `Seurat` object saved as an `.rds` file. The object's metadata (`@meta.data`) must contain the following columns for the pipeline to function correctly:

*   `type`: Tissue origin (e.g., "Cancer", "Adjacent").
*   `patient`: Patient identifier.
*   `celltype`: B-cell subtype annotation (required for `02_heterogeneidade.R`).
*   `clone_id`: Clonal lineage identifier (required for `03_clonalidade_do_repertorio_BCR.R`).
*   `locus`: Gene locus (e.g., "IGH", "IGK", "IGL").
*   `v_call`: V gene annotation.
*   `j_call`: J gene annotation.
*   `c_call`: Constant region (isotype) annotation (e.g., "IGHM", "IGHG1").
*   `v_identity`: Percentage of sequence identity to the germline V gene (required for `04_SHM.R`).
*   `mu_freq`: Somatic hypermutation frequency (required for `04_SHM.R`).
*   `productive`: Boolean indicating if the sequence is productive (TRUE/FALSE).
*   `nFeature_RNA`, `nCount_RNA`: Quality control metrics (optional, used in `01_pre_proc.R`).

## Obtaining the Data
*(Provide instructions here on how a user can obtain the data, e.g., "Data is available upon request from the corresponding author" or "Raw sequencing data can be downloaded from GEO under accession number GSEXXXXXX and processed using the standard 10x Genomics Cell Ranger pipeline.")*
