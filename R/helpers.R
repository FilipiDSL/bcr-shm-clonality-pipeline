# ============================================================
# Helper Functions and Setup
# ============================================================

# Load required packages
suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(vegan)
  library(rstatix)
  library(pheatmap)
})

# Load configuration
load_config <- function(config_path = "configs/config.yml") {
  if (!file.exists(config_path)) {
    stop("Configuration file not found at: ", config_path)
  }
  yaml::read_yaml(config_path)
}

# Ensure directories exist
setup_directories <- function(config) {
  dirs <- c(
    config$paths$data_raw,
    config$paths$data_processed,
    config$paths$results_figures,
    config$paths$results_tables,
    config$paths$results_reports
  )
  for (d in dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
}

# Set seed for reproducibility
set_seed <- function(config) {
  set.seed(config$params$seed)
}
