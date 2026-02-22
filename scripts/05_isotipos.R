# ============================================================
# Mudança de Isótipos (CSR)
# Projeto: Modelagem Computacional de Linfócitos B no COAD
# Autor: Luiz Filipi
# Data: 20/08/2025
# ============================================================

source("R/helpers.R")
config <- load_config()
set_seed(config)

in_path <- file.path(config$paths$data_raw, config$files$qc_tumor_adj)
if (!file.exists(in_path)) {
  stop("Arquivo não encontrado: ", in_path)
}
bcr_data <- readRDS(in_path)

# ------------------------------------
# 1) Extrair metadados relevantes
# ------------------------------------
meta <- bcr_data@meta.data

if (!"c_call" %in% colnames(meta)) {
  message("Aviso: Coluna 'c_call' não encontrada. Pulando análise de isótipos.")
} else {
  iso_df <- meta %>%
    filter(!is.na(c_call)) %>%
    dplyr::select(patient, type, c_call) %>%
    mutate(
      type = factor(type, levels = c("Adjacent", "Cancer")),
      c_call = factor(
        c_call,
        levels = c("IGHD","IGHM","IGHDG","IGHDM",
                   "IGHG1","IGHG2","IGHG3","IGHG4",
                   "IGHA1","IGHA2","IGHE","IGK","IGL")
      )
    )

  # ------------------------------------
  # 2) Frequência absoluta e relativa
  # ------------------------------------
  iso_counts <- iso_df %>%
    group_by(type, c_call) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(type) %>%
    mutate(freq = 100 * n / sum(n)) %>%
    ungroup()

  # ------------------------------------
  # 3) Gráficos lado a lado
  # ------------------------------------
  p_abs <- ggplot(iso_counts, aes(x = c_call, y = n, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Frequência absoluta de isótipos",
         x = "Isótipo", y = "Número de células", fill = "Grupo") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_rel <- ggplot(iso_counts, aes(x = c_call, y = freq, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Frequência relativa de isótipos",
         x = "Isótipo", y = "Proporção (%)", fill = "Grupo") +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(config$paths$results_figures, "05_Isotipos_Absoluto.png"), plot = p_abs, width = 8, height = 6)
  ggsave(file.path(config$paths$results_figures, "05_Isotipos_Relativo.png"), plot = p_rel, width = 8, height = 6)

  # ------------------------------------
  # 4) Teste Qui-quadrado (tabela completa)
  # ------------------------------------
  iso_table <- table(iso_df$type, iso_df$c_call)
  chi_res <- suppressWarnings(chisq.test(iso_table))
  print(chi_res)

  # ------------------------------------
  # 5) Versão com isótipos raros agrupados em "Outros"
  # ------------------------------------
  iso_df_grouped <- iso_df %>%
    mutate(c_call_grouped = ifelse(c_call %in% c("IGHD", "IGHG4"), "Outros", as.character(c_call)),
           c_call_grouped = factor(c_call_grouped))

  iso_counts_grouped <- iso_df_grouped %>%
    group_by(type, c_call_grouped) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(type) %>%
    mutate(freq = 100 * n / sum(n)) %>%
    ungroup()

  p_abs_grouped <- ggplot(iso_counts_grouped, aes(x = c_call_grouped, y = n, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Frequência absoluta de isótipos (agrupado)",
         x = "Isótipo", y = "Número de células", fill = "Grupo") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_rel_grouped <- ggplot(iso_counts_grouped, aes(x = c_call_grouped, y = freq, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Frequência relativa de isótipos (agrupado)",
         x = "Isótipo", y = "Proporção (%)", fill = "Grupo") +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(config$paths$results_figures, "05_Isotipos_Absoluto_Agrupado.png"), plot = p_abs_grouped, width = 8, height = 6)
  ggsave(file.path(config$paths$results_figures, "05_Isotipos_Relativo_Agrupado.png"), plot = p_rel_grouped, width = 8, height = 6)

  # ------------------------------------
  # 6) Testes estatísticos (agrupado)
  # ------------------------------------
  iso_table_grouped <- table(iso_df_grouped$type, iso_df_grouped$c_call_grouped)
  chi_res_grouped <- suppressWarnings(chisq.test(iso_table_grouped))
  print(chi_res_grouped)

  fisher_res <- fisher.test(iso_table_grouped, simulate.p.value = TRUE, B = 10000)
  print(fisher_res)
}

message("✓ 05_isotipos.R concluído.")
