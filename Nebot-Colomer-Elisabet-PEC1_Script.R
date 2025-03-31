
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install(version = "3.20")

BiocManager::available()

library(readr);library(SummarizedExperiment);library(tidyverse); library(xfun); library(ggfortify); library(ggrepel)


#1. Seleccionad y descargad un dataset de metabolómica, que podéis obtener de metabolomicsWorkbench o de este repositorio de GitHub.

human_cachexia <- read_csv("Data/human_cachexia.csv")
str(human_cachexia) # se visualiza la estructura del dataset

#2. Cread un objeto de clase SummarizedExperiment que contenga los datos y los metadatos (información acerca del dataset, sus filas y columnas). La clase SummarizedExperiment es una extensión de ExpressionSet, utilizada por muchas aplicaciones y bases de datos (como es el caso de metabolomicsWorkbench).
#¿Cuáles son sus principales diferencias con la clase ExpressionSet?

# Se seleccionan las dos columnas que contienen la información de los pacientes (metadata)
metadata_cols <- 1:2  

# Se crea la matriz de datos de expresión (excluyendo las columnas de metadata)
expression_data <- as.matrix(human_cachexia[, -(metadata_cols)])  

# Se asignan los nombres de las filas con los IDs de los pacientes
rownames(expression_data) <- human_cachexia$`Patient ID`  

# Se crea los metadatos de las muestras (colData), con los IDs de pacientes como nombres de fila
colData_metadata <- human_cachexia[, metadata_cols] %>%  
  column_to_rownames(var = "Patient ID")  

# Descripción del dataset (metadata general)
descripcion_estudio <- "Cachexia is a complex metabolic syndrome associated with an underlying illness (such as cancer) and characterized by loss of muscle with or without loss of fat mass (Evans et al., 2008). A total of 77 urine samples were collected, 47 from patients with cachexia and 30 from control patients (from the 'specmine.datasets' R package)."

# Creación del objeto SummarizedExperiment
se_object <- SummarizedExperiment(
  assays = list(counts = t(expression_data)),  # Matriz de datos
  colData = colData_metadata, # Metadatos de pacientes
  metadata = list(description = descripcion_estudio))  #  Descripción del estudio

se_object

assay(se_object) # se visualiza la matriz de datos
colData(se_object) #se visualiza la ID de los pacientes y el factor 'Muscle.loss'
metadata(se_object) # se visualiza la descripción del estudio

save(se_object, file = "SummarizedExperiment_data.Rda") # se guarda el 'se_object' como archivo .Rda


#3. Llevad a cabo un análisis exploratorio que os proporcione una visión general del dataset en la línea de lo que hemos visto en las actividades de este reto.

# Análisis univariantes

dim(se_object) # dimensiones de la matriz
names(se_object) # nombres de metabolitos que han sido medidos en cada muestra
summary(assay(se_object, "counts")) # estadística descriptiva de los niveles de cada metabolito en los diferentes pacientes

# Seguidamente para poder graficar con ggplot convierto la matriz en un ‘data.frame’ 
dataframe <- as.data.frame(t(assay(se_object, "counts"))) %>%
  rownames_to_column(var = "Patient_ID") %>%
  pivot_longer(-Patient_ID, names_to = "Metabolito", values_to = "Expresión")

dataframe <- dataframe  %>%
  left_join(as.data.frame(colData(se_object)) %>% 
              rownames_to_column(var = "Patient_ID"), by = "Patient_ID") # se une la columna de 'Muscle.loss'

dim(dataframe) # se visualiza las dimensiones de este dataframe
names(dataframe) # se visualiza los nombres de las columnas del dataframe


# Se crea histogramas por paciente
(histo <- ggplot(dataframe, aes(x = Expresión)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    facet_wrap(~ Patient_ID, scales = "free") + 
    labs(x = "Nivel de expresión", y = "Frecuencia") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)))

ggsave("Outputs/histogramas_expresion.png", plot = histo, width = 30, height = 25, dpi = 300)


# Se crea histogramas por paciente con transformación logaritmica
(histo_trans <- ggplot(dataframe, aes(x = log10(Expresión + 1))) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    facet_wrap(~ Patient_ID, scales = "free") + 
    labs(x = "Nivel de expresión (log transformado)", y = "Frecuencia") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)))

ggsave("Outputs/histogramas_expresion_trans.png", plot = histo_trans, width = 30, height = 25, dpi = 300)


# Se crea un boxplot tranformado para comparar los niveles de expresión por paciente
(boxplot1 <- ggplot(dataframe, aes(x = Patient_ID, y = log10(Expresión + 1), fill = `Muscle.loss`)) +
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = c("cachexic" = "brown4", "control" = "steelblue")) + 
    labs(x = "Paciente", y = "Nivel de expresión", fill = "Muscle.Loss")+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)))

ggsave("Outputs/Boxplot_plot1.png", plot = boxplot1, width = 25, height = 15, dpi = 300)


# Se crea un boxplot tranformado para comparar los niveles de expresión por metabolito
(boxplot2 <- ggplot(dataframe, aes(x = Metabolito, y = log10(Expresión + 1), fill = `Muscle.loss`)) +
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = c("cachexic" = "brown4", "control" = "steelblue")) + 
    labs(x = "Paciente", y = "Nivel de expresión", fill = "Muscle.Loss")+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)))

ggsave("Outputs/Boxplot_plot2.png", plot = boxplot2, width = 20, height = 15, dpi = 300)


# Se crea un boxplot tranformado para comparar los niveles de expresión por afección de caquexia o no
(boxplot3 <- ggplot(dataframe, aes(x = `Muscle.loss`, y = log10(Expresión + 1), fill = `Muscle.loss`)) +
    geom_boxplot(outlier.shape = NA) +  
    scale_fill_manual(values = c("cachexic" = "brown4", "control" = "steelblue")) + 
    labs(x = "Tratamiento", y = "Nivel de expresión (log transformado)", fill = "Muscle.Loss") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    ))

ggsave("Outputs/Boxplot_plot3.png", plot = boxplot3, width = 20, height = 15, dpi = 300)



# Análisis multivariante

#Para visualizar como los pacientes se agrupan espacialmente en función de sus niveles de expresión de metabolitos se realiza una PCA.

# Realizar PCA sobre los datos de expresión de metabolitos
log_counts <- log10(t(assay(se_object, "counts")) + 1) #transformo los datos con el logaritmo en base 10.

pca_res_log <- prcomp(log_counts, scale. = TRUE) # se obtiene las coordenadas de cada PC por metabolito

# Extraer las coordenadas de los primeros dos componentes principales (PC1 y PC2)
pca_df <- as.data.frame(pca_res_log$x)  # Esto extrae la x(coordenadas de PC) y convierte las coordenadas en dataframe
pca_df$Muscle.loss <- colData(se_object)$`Muscle loss`  # Se añade la columna de metadatos (Muscle.loss)

# Extraer las coordenadas de los metabolitos
pca_var <- as.data.frame(pca_res_log$rotation[, 1:2]) #  extrae las columnas del PC1 y PC2 y lo convierte en dataframe
pca_var$Metabolito <- rownames(pca_var)  # Añade una columna con los nombres de metabolitos


## Determinar el porcentaje de variación explicado por cada PC.
std_devs <- pca_res_log$sdev # extrae los desviación estándar del objeto pca_res_log
eigenvalues <- std_devs^2 # cálculo para obtener la varianza.
variance_explained <- eigenvalues / sum(eigenvalues)
cumulative_variance <- cumsum(variance_explained) # porcentaje de variación explicada por los PC.

# PCA del PC1 y Pc2 

(pca_plot_v1 <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Muscle.loss)) + 
    geom_point(size = 3) + 
    scale_colour_manual(values = c("cachexic" = "brown4", "control" = "steelblue")) +
    labs(x = "PC1 (58.53%)",y = "PC2 (4.45%)") +
    theme_bw() +
    theme(text = element_text(size = 16), 
          legend.title = element_text(size = 18), 
          legend.text = element_text(size = 16)))

ggsave("Outputs/pca_plot_v1.png", plot = pca_plot_v1, width = 15, height = 8, dpi = 300)


# PCA añadiendo las coordenadas de los metabolitos
(pca_plot_v2 <- ggplot() +
    geom_point(data = pca_df, aes(x = PC1, y = PC2, colour = Muscle.loss), size = 3) + 
    scale_colour_manual(values = c("cachexic" = "brown4", "control" = "steelblue")) +
    geom_segment(data = pca_var, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),     # Flechas de los metabolitos
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.7) +
    geom_text_repel(data = pca_var, aes(x = PC1 * 5, y = PC2 * 5, label = Metabolito), 
                    size = 5, color = "black", max.overlaps = 20) +     # Etiquetas de los metabolitos
    labs(x = "PC1 (58.53%)",
         y = "PC2 (4.45%)") + 
    theme_bw() +
    theme(text = element_text(size = 16), 
          legend.title = element_text(size = 18), 
          legend.text = element_text(size = 16)))


ggsave("Outputs/pca_plot_v2.png", plot = pca_plot_v2, width = 15, height = 8, dpi = 300)

