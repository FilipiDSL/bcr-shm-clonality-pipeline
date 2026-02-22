# ============================================================
# Pré-processamento e QC Tumor vs Adjacent
# Projeto: Modelagem Computacional de Linfócitos B no COAD
# Autor: Luiz Filipi
# Data: 20/08/2025
# ============================================================

source("R/helpers.R")
config <- load_config()
set_seed(config)

# ------------------------------------------------------------
# 1. Carregar objeto filtrado do Dia 1
# ------------------------------------------------------------
in_path <- file.path(config$paths$data_raw, config$files$filtered_tumor_adj)
if (!file.exists(in_path)) {
  stop("Arquivo não encontrado: ", in_path)
}
coad_tumor_adj <- readRDS(in_path)
print(coad_tumor_adj)

# ------------------------------------------------------------
# 2. QC básico (opcional: número de genes, UMI, %MT se disponível)
# ------------------------------------------------------------
if (all(c("nFeature_RNA", "nCount_RNA") %in% colnames(coad_tumor_adj@meta.data))) {
  p_qc <- VlnPlot(coad_tumor_adj, features = c("nFeature_RNA", "nCount_RNA"), group.by = "type")
  ggsave(file.path(config$paths$results_figures, "01_QC_VlnPlot.png"), plot = p_qc, width = 10, height = 5)
}

# ------------------------------------------------------------
# 3. Normalização
# ------------------------------------------------------------
coad_tumor_adj <- NormalizeData(coad_tumor_adj, 
                                normalization.method = config$params$normalization_method, 
                                scale.factor = config$params$scale_factor)

# ------------------------------------------------------------
# 4. Identificação de variáveis + PCA
# ------------------------------------------------------------
coad_tumor_adj <- FindVariableFeatures(coad_tumor_adj, selection.method = "vst", nfeatures = config$params$nfeatures_vst)
coad_tumor_adj <- ScaleData(coad_tumor_adj, verbose = FALSE)
coad_tumor_adj <- RunPCA(coad_tumor_adj, npcs = config$params$npcs, verbose = FALSE)

# ------------------------------------------------------------
# 5. UMAP inicial
# ------------------------------------------------------------
coad_tumor_adj <- RunUMAP(coad_tumor_adj, dims = 1:config$params$umap_dims)
p_umap <- DimPlot(coad_tumor_adj, reduction = "umap", group.by = "type", pt.size = 0.5) +
  ggtitle("COAD: Tumor vs Adjacent")

ggsave(file.path(config$paths$results_figures, "01_UMAP_QC.png"), plot = p_umap, width = 8, height = 6)

# ------------------------------------------------------------
# 6. Salvar objeto processado
# ------------------------------------------------------------
out_path <- file.path(config$paths$data_raw, config$files$qc_tumor_adj)
saveRDS(coad_tumor_adj, out_path)
message("✓ 01_pre_proc.R concluído. Objeto salvo em: ", out_path)
