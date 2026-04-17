# AlphaFold3 pLDDT + Rosetta combined ranking analysis
# Autor: Xzamu · Lab. Fisicoquímica e Ingeniería de Proteínas
#
# 1) Lee todas las subcarpetas de AF3 (t01_*, t02_*, t03_*)
# 2) Extrae atom_plddts de cada full_data_<n>.json
# 3) Promedia por secuencia, ordena, exporta resultados_promedios_plddt.csv
# 4) Toma top 10 → top_10_secuencias_plddt.csv
# 5) Cruza con ranking_top1.tsv de Rosetta y calcula ranking combinado
#    en dos escenarios (50/50 y 60/40) → top50_por_escenario.csv

ruta_base <- "~/Escritorio/Protein_Engineering/ProteinMPNN-outputs_BN/AlphaFold-T01-T03"
carpetas_modelos <- list.dirs(path = ruta_base, recursive = FALSE, full.names = TRUE)

resultados <- data.frame(Modelo = character(),
                         Grupo = character(),
                         Promedio_pLDDT = numeric())

calcular_promedio_plddt <- function(carpeta) {
  archivos <- list.files(path = carpeta,
                         pattern = "full_data_\\d+.json$",
                         full.names = TRUE)
  todos_plddt <- c()
  for (archivo in archivos) {
    json <- jsonlite::fromJSON(archivo)
    plddts <- unlist(json$atom_plddts)
    todos_plddt <- c(todos_plddt, plddts)
  }
  if (length(todos_plddt) == 0) return(NA)
  mean(todos_plddt)
}

for (carpeta in carpetas_modelos) {
  promedio <- calcular_promedio_plddt(carpeta)
  nombre_carpeta <- basename(carpeta)
  grupo <- substr(nombre_carpeta, 1, 3)   # t01, t02, t03
  resultados <- rbind(resultados,
                      data.frame(Modelo = nombre_carpeta,
                                 Grupo = grupo,
                                 Promedio_pLDDT = promedio))
}

resultados <- resultados[order(-resultados$Promedio_pLDDT), ]
write.csv(resultados, file = "resultados_promedios_plddt.csv", row.names = FALSE)

top_10_secuencias <- head(resultados, 10)
write.csv(top_10_secuencias, file = "top_10_secuencias_plddt.csv", row.names = FALSE)

#############################################
# ANÁLISIS DE DATOS FINAL — combinación Rosetta
#############################################
library(readr)
library(dplyr)
library(stringr)
library(tibble)

setwd("~/Escritorio/Protein_Engineering/Fitasa_PyRosetta/input")

top_10_secuencias <- top_10_secuencias %>%
  rename(design = Modelo, plddt_mean = Promedio_pLDDT) %>%
  mutate(design = str_trim(design))

# Lee ranking Rosetta (mejor pose por secuencia)
rosetta <- read.delim("ranking_top1.tsv",
                      header = TRUE,
                      sep = "",
                      comment.char = "#",
                      check.names = FALSE) |>
  dplyr::rename(design_raw = design,
                best_model = top_model,
                rosetta_E  = total_score)

rosetta <- rosetta %>%
  mutate(design = str_remove(design_raw, "_\\d+$"))

merged <- inner_join(top_10_secuencias, rosetta, by = "design")

# Normalización 0-1 de ambas métricas
max_E <- max(merged$rosetta_E)
min_E <- min(merged$rosetta_E)

merged <- merged %>%
  mutate(plddt_scaled   = (plddt_mean - min(plddt_mean)) /
                          (max(plddt_mean) - min(plddt_mean)),
         rosetta_scaled = (max_E - rosetta_E) / (max_E - min_E))

# Dos escenarios de ponderación
w <- tibble(w_rosetta = c(0.50, 0.60),
            w_plddt   = 1 - w_rosetta,
            escenario = c("50-50", "60-40"))

rankings <- lapply(1:nrow(w), function(i) {
  merged %>%
    mutate(score_comb = w$w_rosetta[i] * rosetta_scaled +
                        w$w_plddt[i]   * plddt_scaled) %>%
    arrange(desc(score_comb)) %>%
    mutate(rank = row_number(),
           escenario = w$escenario[i]) %>%
    select(escenario, rank, everything())
}) %>% bind_rows()

top50 <- rankings %>%
  as_tibble() %>%
  filter(rank <= 50) %>%
  select(escenario, rank, design, best_model,
         plddt_mean, rosetta_E, score_comb) %>%
  arrange(escenario, rank)

write_csv(top50, "top50_por_escenario.csv")
cat("Archivo 'top50_por_escenario.csv' guardado en", getwd(), "\n")
