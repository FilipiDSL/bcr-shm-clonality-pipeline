# ============================================================
# Preparação e Organização
# Projeto: Modelagem Computacional de Linfócitos B no COAD
# Autor: Luiz Filipi
# Data: 19/08/2025
# ============================================================

# Carregar helpers e configurações
source("R/helpers.R")
config <- load_config()
setup_directories(config)
set_seed(config)

# ------------------------------------------------------------
# 1. CARREGAR DADOS .RDS FILTRADOS PARA COAD
# ------------------------------------------------------------
# Verifica se o arquivo bruto original existe. Se não, pula este script
# pois o usuário pode já ter fornecido o arquivo filtrado.
raw_file_path <- file.path(config$paths$data_raw, config$files$raw_bcr_preprocessed)

if (!file.exists(raw_file_path)) {
  message("Arquivo bruto original não encontrado: ", raw_file_path)
  message("Pulando 00_prep.R. Assumindo que ", config$files$filtered_tumor_adj, " já está disponível.")
} else {
  message("Carregando dados brutos...")
  coad <- readRDS(raw_file_path)
  print(coad)
  
  # Ver as colunas de metadados disponíveis
  print(colnames(coad@meta.data))
  
  # ------------------------------------------------------------
  # 2. FILTRAR PARA TUMOR E ADJACENT
  # ------------------------------------------------------------
  # A coluna 'type' contém as categorias: "Cancer", "Adjacent", "Blood", "LN_Met" etc.
  print(table(coad@meta.data$type))
  
  # Filtrar apenas células de tumor e normal adjacente
  coad_tumor_adj <- subset(coad, subset = type %in% c("Cancer", "Adjacent"))
  
  # Salvar objeto filtrado
  out_path <- file.path(config$paths$data_raw, config$files$filtered_tumor_adj)
  saveRDS(coad_tumor_adj, out_path)
  message("Objeto filtrado salvo em: ", out_path)
  
  # ------------------------------------------------------------
  # 3. MINI-RELATÓRIO EXPLORATÓRIO
  # ------------------------------------------------------------
  # Quantas células em cada tipo (Tumor vs Adjacent)?
  print(table(coad_tumor_adj@meta.data$type))
  
  # Quantos pacientes representados em cada tipo?
  print(table(coad_tumor_adj@meta.data$patient, coad_tumor_adj@meta.data$type))
  
  # Quantas células por paciente (resumo rápido)
  resumo_pacientes <- coad_tumor_adj@meta.data %>%
    group_by(patient, type) %>%
    summarise(cells = n(), .groups = "drop") %>%
    arrange(desc(cells)) %>%
    head(10)
  print(resumo_pacientes)
  
  # ------------------------------------------------------------
  # 4. VISUALIZAÇÃO EXPLORATÓRIA
  # ------------------------------------------------------------
  # UMAP já vem calculado no objeto -> plot rápido colorido por tipo
  p <- DimPlot(coad_tumor_adj, group.by = "type", reduction = "umap") +
    ggtitle("Distribuição inicial das células B - Tumor vs Adjacent (COAD)")
  
  ggsave(file.path(config$paths$results_figures, "00_UMAP_inicial.png"), plot = p, width = 8, height = 6)
  message("✓ 00_prep.R concluído.")
}
