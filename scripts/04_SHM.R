# ============================================================
# SHM
# Projeto: Modelagem Computacional de Linfócitos B no COAD
# Autor: Luiz Filipi
# Data: 20/08/2025
# ============================================================

source("R/helpers.R")
config <- load_config()
set_seed(config)

# ====================================
# Análise BCR: Identidade V e SHM
# ====================================

in_path <- file.path(config$paths$data_raw, config$files$qc_tumor_adj)
if (!file.exists(in_path)) {
  stop("Arquivo não encontrado: ", in_path)
}
bcr_data <- readRDS(in_path)

# ------------------------------------
# 1. Extração dos metadados do objeto
# ------------------------------------
meta <- bcr_data@meta.data

req_cols <- c("patient", "type", "celltype", "clone_id", "c_call", 
              "v_identity", "mu_freq", "umi_count", "productive", "locus")

if (!all(req_cols %in% colnames(meta))) {
  message("Aviso: Colunas necessárias para SHM não encontradas. Pulando análise.")
} else {
  # Selecionar colunas principais
  bcr_df <- meta %>%
    dplyr::select(all_of(req_cols)) %>%
    filter(productive == TRUE, locus == "IGH") %>%
    filter(!is.na(c_call), !is.na(v_identity), !is.na(mu_freq))

  # ------------------------------------
  # 2. Boxplot: Identidade V por isótipo e tipo
  # ------------------------------------
  p1 <- ggplot(bcr_df, aes(x = c_call, y = v_identity, fill = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.3, size = 0.7) +
    scale_fill_manual(values = c("Adjacent" = "#F8766D", "Cancer" = "#00BFC4")) +
    scale_color_manual(values = c("Adjacent" = "#F8766D", "Cancer" = "#00BFC4")) +
    theme_minimal(base_size = 14) +
    labs(title = "Comparação de Identidade V por Isótipo (Tumor vs Adjacent)",
         x = "Isótipo (c_call)",
         y = "Identidade V (% de similaridade com germline)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(config$paths$results_figures, "04_V_Identity_Boxplot.png"), plot = p1, width = 10, height = 6)

  # ------------------------------------
  # 3. Heatmap: SHM (mu_freq) por paciente × isótipo × tipo
  # ------------------------------------
  heatmap_df <- bcr_df %>%
    group_by(type, patient, c_call) %>%
    summarise(mu_mean = mean(mu_freq, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = c_call, values_from = mu_mean, values_fill = 0)

  # --- Tumor ---
  heatmap_tumor <- heatmap_df %>% filter(type == "Cancer") %>% as.data.frame()
  if (nrow(heatmap_tumor) > 0) {
    rownames(heatmap_tumor) <- heatmap_tumor$patient
    heatmap_tumor <- heatmap_tumor[, !(colnames(heatmap_tumor) %in% c("patient","type"))]
    
    png(file.path(config$paths$results_figures, "04_SHM_Heatmap_Tumor.png"), width = 800, height = 600)
    pheatmap(heatmap_tumor,
             scale = "row",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "ward.D2",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             main = "Frequência média de SHM (mu_freq) - TUMOR")
    dev.off()
  }
}

message("✓ 04_SHM.R concluído.")
