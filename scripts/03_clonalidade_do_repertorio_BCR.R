# ============================================================
# Clonalidade do Repertório BCR
# Projeto: Modelagem Computacional de Linfócitos B no COAD
# Autor: Luiz Filipi
# Data: 22/08/2025
# ============================================================

source("R/helpers.R")
config <- load_config()
set_seed(config)

# Carregar objeto já pré-processado
in_path <- file.path(config$paths$data_raw, config$files$qc_tumor_adj)
if (!file.exists(in_path)) {
  stop("Arquivo não encontrado: ", in_path)
}
bcr_data <- readRDS(in_path)

# Frequência de loci (pesada vs leve)
if ("locus" %in% colnames(bcr_data@meta.data)) {
  print(table(bcr_data@meta.data$locus))
}

# Conferir exemplos de genes anotados (V e J)
if (all(c("v_call", "j_call") %in% colnames(bcr_data@meta.data))) {
  print(head(bcr_data@meta.data$v_call))
  print(head(bcr_data@meta.data$j_call))
}

# Resumo dos genes de cadeia constante (isótipos)
if ("c_call" %in% colnames(bcr_data@meta.data)) {
  print(table(bcr_data@meta.data$c_call))
}

# Extrair metadados
if (!all(c("type", "clone_id", "patient") %in% colnames(bcr_data@meta.data))) {
  message("Aviso: Colunas necessárias (type, clone_id, patient) não encontradas. Pulando análise de clonalidade.")
} else {
  meta <- bcr_data@meta.data %>%
    dplyr::select(type, clone_id, patient) %>%
    tibble::rownames_to_column("barcode")

  # ------------------------------------------------------------
  # 1. Top 10 clones expandidos em Tumor vs Adjacent
  # ------------------------------------------------------------
  clone_counts <- meta %>%
    group_by(type, clone_id) %>%
    summarise(Frequency = n(), .groups = "drop") %>%
    arrange(type, desc(Frequency))

  top10_clones <- clone_counts %>%
    group_by(type) %>%
    slice_max(order_by = Frequency, n = 10)

  p_top10 <- ggplot(top10_clones, aes(x = reorder(clone_id, Frequency), y = Frequency, fill = type)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~type, scales = "free_y") +
    labs(title = "Top 10 clones expandidos por tecido",
         x = "Clone ID", y = "Número de células") +
    scale_fill_manual(values = c("Cancer"="#d73027", "Adjacent"="#4575b4")) +
    theme_bw(base_size = 13)
  
  ggsave(file.path(config$paths$results_figures, "03_Top10_Clones.png"), plot = p_top10, width = 10, height = 6)

  # ------------------------------------------------------------
  # 2. Diversidade clonal (Shannon / Simpson) por paciente
  # ------------------------------------------------------------
  calc_diversity <- function(df) {
    freq_table <- table(df$clone_id)
    data.frame(
      Shannon = diversity(freq_table, index = "shannon"),
      Simpson = diversity(freq_table, index = "simpson")
    )
  }

  diversity_stats <- meta %>%
    group_by(type, patient) %>%
    group_modify(~calc_diversity(.x)) %>%
    ungroup()

  # ------------------------------------------------------------
  # 3. Teste estatístico (Wilcoxon)
  # ------------------------------------------------------------
  p_values <- diversity_stats %>%
    tidyr::pivot_longer(cols = c("Shannon","Simpson"),
                        names_to = "Index", values_to = "Value") %>%
    group_by(Index) %>%
    wilcox_test(Value ~ type) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()

  print(p_values)
  write.csv(p_values, file.path(config$paths$results_tables, "03_Diversity_Wilcoxon.csv"), row.names = FALSE)

  # ------------------------------------------------------------
  # 4. Boxplots com jitter (visualização Nature-like)
  # ------------------------------------------------------------
  div_df <- tidyr::pivot_longer(diversity_stats, cols = c("Shannon","Simpson"),
                                names_to = "Index", values_to = "Value")

  p_div <- ggplot(div_df, aes(x = type, y = Value, fill = type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    facet_wrap(~Index, scales = "free_y") +
    labs(title = "Diversidade clonal BCR por tecido",
         x = "", y = "Índice de diversidade") +
    scale_fill_manual(values = c("Cancer"="#d73027", "Adjacent"="#4575b4")) +
    theme_bw(base_size = 13)
  
  ggsave(file.path(config$paths$results_figures, "03_Diversity_Boxplots.png"), plot = p_div, width = 8, height = 6)
}

message("✓ 03_clonalidade_do_repertorio_BCR.R concluído.")
