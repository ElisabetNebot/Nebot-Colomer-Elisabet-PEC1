
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install(version = "3.20")

BiocManager::available()

library(readr);library(SummarizedExperiment);library(tidyverse)

human_cachexia <- read_csv("Data/human_cachexia.csv")
#str(human_cachexia)


metadata_cols <- 1:2  ## se seleccionan las dos columnas que contienen la información metadata
expression_data <- as.matrix(human_cachexia[, -(metadata_cols)])  # Los siguientes datos (todos- colmetadata) se convierten a una matriz

rownames(expression_data) <- human_cachexia$`Patient ID` ## Asignar los nombres de las filas como los IDs de los pacientes

# Crear metadata de muestras (colData)
col_metadata <- human_cachexia[, metadata_cols] %>% 
  column_to_rownames(var = "Patient ID")  # Usar 'Patient ID' como rownames

# Metadata descripción del estudio
description = "Cachexia is a complex metabolic syndrome associated with an underlying illness (such as cancer) and characterized by loss of muscle with or without loss of fat mass (Evans et al., 2008). A total of 77 urine samples were collected, 47 from patients with cachexia and 30 from control patients (from the 'specmine.datasets' R package)."

# Crear el objeto SummarizedExperiment
se_object <- SummarizedExperiment(
  assays = list(counts = t(expression_data)),  # Transponer para que filas sean metabolitos
  colData = col_metadata,
  metadata = list(description = description))

se_object

save(se_object, file = "SummarizedExperiment_data.Rda")




