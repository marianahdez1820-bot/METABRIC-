library(tidyverse)
library(cluster)
library(ggplot2)
library(factoextra)
library(NbClust)
library(tidyr)
library(pheatmap)

# Leer datos
brca_data <- read_tsv("data_mrna_illumina_microarray.txt", col_names = TRUE)

colnames(brca_data) # Son las muestras con MB
View(brca_data)

# Quitar columnas que no sean numéricas (las dos primeras del nombre y ID)
brca_data.num <- as.matrix(brca_data[, -c(1,2)])  
storage.mode(brca_data.num) <- "double" # forzar a que sea numérico 

# Sacar genes con más variación
var <- apply(brca_data.num, 1, var) # calcular varianza de cada gen en todas las muestras
top  <- order(var, decreasing = TRUE)[1:500] # elegir 500 genes más variables
brca_data.fil <- brca_data.num[top, ] # df solo con esos 500 genes

View(brca_data.fil)

# Voltear para que en filas quede MB y en columnas los genes porque el clustering agrupa filas
data_brca.ord <- t(brca_data.fil)

View(data_brca.ord)

# Normalizar
data_brca.norm <- scale(data_brca.ord)


# ---------------- CLUSTERING JERÁRQUICO ----------------

#Matriz de distancias entre muestras
 dist_brca <- dist(data_brca.norm, method = "euclidean")

# Clustering jerárquico
hc_brca <- hclust(dist_brca, method = "ward.D2")


# Elegir número de clusters
fviz_nbclust(data_brca.norm, FUN = hcut, method = "silhouette")


# Graficar dendrograma
plot(hc_brca, labels = FALSE, hang = -1,
     main = "Clustering jerárquico de muestras METABRIC",
     xlab = "", sub = "")

rect.hclust(hc_brca, k = 2, border = "red")



# Cortar  árbol para obtener clusters
hc.clusters <- cutree(hc_brca, k = 2)



# Ver cuántas muestras quedaron en cada cluster
table(hc.clusters)



# Visualizar clusters
set.seed(123)
fviz_cluster(list(data = data_brca.norm, cluster = hc.clusters),
            geom = "point",
            ellipse.type = "convex",
            main = "Clustering jerárquico",
            xlab = "Dimensión 1",
            ylab = "Dimensión 2",
            palette = "Set1")

?fviz_cluster


# Guardar dataframe 
cluster_df <- data.frame(
  sample_id = rownames(data_brca.norm),
  cluster = hc.clusters)

View(cluster_df)




