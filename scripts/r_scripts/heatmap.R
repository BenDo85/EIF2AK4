setwd(dir = "~/Bureau/projet/data")

install.packages("heatmaply")
library(heatmaply)

luad_data <- readRDS("all_data_luad.rds")
coad_data <- readRDS("all_coad_data.rds")
brca_data <- readRDS("all_data_brca.rds") 
kegg_genes <- read.table("genes.txt",sep=',')



kegg_genes <- as.data.frame(kegg_genes)
kegg_genes <- data.frame(t(kegg_genes[-1]))

#LUAD 

gheatmapdf_luad <- subset(luad_data, luad_data$Hugo_Symbol %in% kegg_genes$X1)

luad_colnames <- unique(gheatmapdf_luad$Hugo_Symbol)
luad_rows <- unique(gheatmapdf_luad$Patient)

heatmap_df_luad <- matrix(nrow = length(luad_rows), ncol = length(luad_colnames))
colnames(heatmap_df_luad) = luad_colnames
rownames(heatmap_df_luad) = luad_rows

for (i in rownames(heatmap_df_luad)){
  for (j in colnames(heatmap_df_luad))
  { 
    if(length(gheatmapdf_luad$MRNA[gheatmapdf_luad$Patient==i & gheatmapdf_luad$Hugo_Symbol == j]) >0) {
      heatmap_df_luad[i,j] <- mean(gheatmapdf_luad$MRNA[gheatmapdf_luad$Hugo_Symbol==j & gheatmapdf_luad$Patient == i])
    }
    
  }
  
}

#COAD

gheatmapdf_coad <- subset(coad_data, coad_data$Hugo_Symbol %in% kegg_genes$X1)

coad_colnames <- unique(gheatmapdf_coad$Hugo_Symbol)
coad_rows <- unique(gheatmapdf_coad$Patient)

heatmap_df_coad <- matrix(nrow = length(coad_rows), ncol = length(coad_colnames))
colnames(heatmap_df_coad) = coad_colnames
rownames(heatmap_df_coad) = coad_rows

for (i in rownames(heatmap_df_coad)){
  for (j in colnames(heatmap_df_coad))
  { 
    if(length(gheatmapdf_coad$MRNA[gheatmapdf_coad$Patient==i & gheatmapdf_coad$Hugo_Symbol == j]) >0) {
      heatmap_df_coad[i,j] <- mean(gheatmapdf_coad$MRNA[gheatmapdf_coad$Hugo_Symbol==j & gheatmapdf_coad$Patient == i])
    }
    
  }
}


#BRCA 

gheatmapdf_brca <- subset(brca_data, brca_data$Hugo_Symbol %in% kegg_genes$X1)

brca_colnames <- unique(gheatmapdf_brca$Hugo_Symbol)
brca_rows <- unique(gheatmapdf_brca$Patient)

heatmap_df_brca <- matrix(nrow = length(brca_rows), ncol = length(brca_colnames))
colnames(heatmap_df_brca) = brca_colnames
rownames(heatmap_df_brca) = brca_rows

for (i in rownames(heatmap_df_brca)){
  for (j in colnames(heatmap_df_brca))
  { 
    if(length(gheatmapdf_brca$MRNA[gheatmapdf_brca$Patient==i & gheatmapdf_brca$Hugo_Symbol == j]) >0) {
      heatmap_df_brca[i,j] <- mean(gheatmapdf_brca$MRNA[gheatmapdf_brca$Hugo_Symbol==j & gheatmapdf_brca$Patient == i])
    }
    
  }
}

par(mfrow=c(2,2))
heatmap(heatmap_df_luad)
heatmap(heatmap_df_coad)
heatmap(heatmap_df_brca)
