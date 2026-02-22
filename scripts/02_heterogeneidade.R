# ============================================================
# Heterogeneidade
# Projeto: Modelagem Computacional de Linfócitos B no COAD
# Autor: Luiz Filipi
# Data: 20/08/2025
# ============================================================

source("R/helpers.R")
config <- load_config()
set_seed(config)

# ------------------------------------------------------------
# 1. Carregar objeto filtrado e processado do Dia 2
# ------------------------------------------------------------
in_path <- file.path(config$paths$data_raw, config$files$qc_tumor_adj)
if (!file.exists(in_path)) {
  stop("Arquivo não encontrado: ", in_path)
}
bcr_obj <- readRDS(in_path)

# ------------------------------------------------------------
# 2. UMAP colorido por subtipos de células B
# ------------------------------------------------------------
if ("celltype" %in% colnames(bcr_obj@meta.data)) {
  p_umap <- DimPlot(bcr_obj, group.by = "celltype", reduction = "umap", 
          label = TRUE, repel = TRUE) +
    ggtitle("UMAP de Subtipos de Células B")
  ggsave(file.path(config$paths$results_figures, "02_UMAP_Subtipos.png"), plot = p_umap, width = 8, height = 6)

  # ------------------------------------------------------------
  # 3. Calcular proporções Tumor vs Adjacent
  # ------------------------------------------------------------
  prop_df <- table(bcr_obj$type, bcr_obj$celltype) %>%
    prop.table(margin = 1) %>%
    as.data.frame()

  colnames(prop_df) <- c("Tissue", "Subtype", "Proportion")

  # ------------------------------------------------------------
  # 4. Gráfico de barras empilhadas
  # ------------------------------------------------------------
  p_bar_fill <- ggplot(prop_df, aes(x = Tissue, y = Proportion, fill = Subtype)) +
    geom_bar(stat = "identity", position = "fill") +
    ylab("Proporção relativa") +
    xlab("Origem") +
    ggtitle("Distribuição de Subtipos de Células B: Tumor vs Adjacent") +
    theme_minimal()
  ggsave(file.path(config$paths$results_figures, "02_Bar_Subtipos_Fill.png"), plot = p_bar_fill, width = 8, height = 6)

  # ------------------------------------------------------------
  # 5. Gráfico de barras lado a lado
  # ------------------------------------------------------------
  p_bar_dodge <- ggplot(prop_df, aes(x = Subtype, y = Proportion, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    ylab("Proporção relativa") +
    xlab("Subtipos de células B") +
    ggtitle("Comparação Tumor vs Adjacent por Subtipo B") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  ggsave(file.path(config$paths$results_figures, "02_Bar_Subtipos_Dodge.png"), plot = p_bar_dodge, width = 10, height = 6)
} else {
  message("Aviso: Coluna 'celltype' não encontrada nos metadados. Pulando gráficos de heterogeneidade.")
}

message("✓ 02_heterogeneidade.R concluído.")
