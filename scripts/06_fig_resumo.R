# ============================================
# Dia 7 — Figura-resumo (mosaico TLS)
# Recompõe: Isótipos (%) | V identity (SHM) | Top clones
# Exporta PNG e PDF
# ============================================

source("R/helpers.R")
config <- load_config()
set_seed(config)

in_path <- file.path(config$paths$data_raw, config$files$qc_tumor_adj)
if (!file.exists(in_path)) {
  stop("Arquivo não encontrado: ", in_path)
}
bcr_data <- readRDS(in_path)

# ------------------------------------------------
# 1) Bases de dados a partir do Seurat@meta.data
# ------------------------------------------------
meta <- bcr_data@meta.data

req_cols <- c("type", "c_call", "productive", "locus", "v_identity", "clone_id")
if (!all(req_cols %in% colnames(meta))) {
  message("Aviso: Colunas necessárias não encontradas. Pulando figura resumo.")
} else {
  # Normalizar rótulos de grupo (só por garantia)
  meta$type <- as.character(meta$type)
  meta$type[meta$type %in% c("Tumor", "Cancer", "tumor", "cancer")] <- "Cancer"
  meta$type[meta$type %in% c("Normal", "Adjacent", "adjacent")] <- "Adjacent"
  meta$type <- factor(meta$type, levels = c("Adjacent", "Cancer"))

  # -----------------------------
  # 1A) Dados para isótipos (CSR)
  # -----------------------------
  iso_df <- meta %>%
    filter(!is.na(c_call)) %>%
    select(patient, type, c_call)

  iso_counts <- iso_df %>%
    count(type, c_call, name = "n")

  p_iso_rel <- ggplot(iso_counts, aes(x = type, y = n, fill = c_call)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = percent_format()) +
    labs(title = "CSR: Distribuição de Isótipos (%)",
         x = NULL, y = "Proporção", fill = "Isótipo") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")

  # -------------------------------------------------------
  # 1B) Dados para SHM (V identity) em células B produtivas
  # -------------------------------------------------------
  bcr_shm <- meta %>%
    filter(productive == TRUE, locus == "IGH", !is.na(v_identity),
           !is.na(type))

  p_videntity <- ggplot(bcr_shm, aes(x = type, y = v_identity, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    labs(title = "SHM: V identity por grupo",
         x = NULL, y = "V identity (similaridade à germline)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")

  # ---------------------------------------------
  # 1C) Dados para Top Clones (expansão clonal)
  # ---------------------------------------------
  top10_clones <- meta %>%
    select(type, clone_id) %>%
    filter(!is.na(type), !is.na(clone_id)) %>%
    count(type, clone_id, name = "Frequency") %>%
    arrange(type, desc(Frequency)) %>%
    group_by(type) %>%
    slice_max(order_by = Frequency, n = 10) %>%
    ungroup()

  p_topclones <- ggplot(top10_clones,
                        aes(x = reorder(clone_id, Frequency), y = Frequency, fill = type)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ type, scales = "free_y") +
    labs(title = "Expansão clonal: Top 10 clones por grupo",
         x = "Clone ID", y = "Nº de células") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")

  # ---------------------------------------------
  # 2) Mosaico (Painel A): 3 gráficos lado a lado
  # ---------------------------------------------
  p_mosaic_A <- p_iso_rel | p_videntity | p_topclones
  p_mosaic_A <- p_mosaic_A + plot_annotation(
    title = "Evidências quantitativas compatíveis com TLS maduros no Tumor",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

  p_final <- p_mosaic_A

  # -------------------------
  # 4) Exportar PNG e PDF
  # -------------------------
  ggsave(file.path(config$paths$results_figures, "06_Dia7_TLS_mosaico.png"),
         plot = p_final, width = 16, height = 6.5, dpi = 300)

  ggsave(file.path(config$paths$results_figures, "06_Dia7_TLS_mosaico.pdf"),
         plot = p_final, width = 16, height = 6.5, device = cairo_pdf)

  # (Opcional) salvar os painéis individualmente
  ggsave(file.path(config$paths$results_figures, "06_painel_isotipos.png"),  p_iso_rel,   width = 5.5, height = 5, dpi = 300)
  ggsave(file.path(config$paths$results_figures, "06_painel_videntity.png"), p_videntity, width = 5.5, height = 5, dpi = 300)
  ggsave(file.path(config$paths$results_figures, "06_painel_topclones.png"), p_topclones, width = 6,   height = 5, dpi = 300)
}

message("✓ 06_fig_resumo.R concluído.")
