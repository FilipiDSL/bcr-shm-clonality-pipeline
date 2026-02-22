# ============================================================
# Runner Script
# Executa todo o pipeline em sequência
# ============================================================

message("============================================================")
message("Iniciando Pipeline BCR SHM Clonality")
message("============================================================")

# Carregar helpers e configurações
source("R/helpers.R")
config <- load_config()
setup_directories(config)

# Lista de scripts na ordem de execução
scripts_to_run <- c(
  "scripts/00_prep.R",
  "scripts/01_pre_proc.R",
  "scripts/02_heterogeneidade.R",
  "scripts/03_clonalidade_do_repertorio_BCR.R",
  "scripts/04_SHM.R",
  "scripts/05_isotipos.R",
  "scripts/06_fig_resumo.R"
)

# Executar cada script
for (script in scripts_to_run) {
  if (file.exists(script)) {
    message("\n------------------------------------------------------------")
    message("Executando: ", script)
    message("------------------------------------------------------------")
    
    # Usar tryCatch para capturar erros e parar a execução se necessário
    tryCatch({
      source(script)
    }, error = function(e) {
      stop("Erro ao executar ", script, ":\n", e$message)
    })
    
  } else {
    stop("Script não encontrado: ", script)
  }
}

message("\n============================================================")
message("Pipeline concluído com sucesso!")
message("Resultados salvos em: ", config$paths$results_figures, " e ", config$paths$results_tables)
message("============================================================")
