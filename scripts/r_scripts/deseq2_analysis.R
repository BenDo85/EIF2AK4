# Ce script fait l'ensemble des analyses avec des graphes et tableau de sortie
# enregistré dans le répertoire result
# Installer et charger les packages nécessaires
required_packages <- c("BiocManager", "dplyr", "data.table", "VennDiagram")

for (package in required_packages) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Il faut au préalable avoir télécharger, regrouper et transformer à partir de GDC cancer (starcount_download.R)
# Ils sont aussi téléchargable à cette adresse https://drive.google.com/drive/folders/1pGrPft4rKYEvKlXKsfDUkKt6ZOexXwv0?usp=sharing
# Les fichiers tsv sont à placer dans le répertoir data
# Il faut aussi téléchargé les études avec download.R

# Charger les données Star Counts de l'étude que l'on souhaite
star_counts_dt <- fread("data/star_count_LUAD.tsv", header = TRUE, sep = '\t', data.table = TRUE, select = c("gene_id", "gene_name", "unstranded", "UUID"))

# On crée un paire entre l'id du gène et son nom usuel
star_counts_dt$gene <- paste(star_counts_dt$gene_id, star_counts_dt$gene_name, sep = "_")

# On réorganise le tableau pour qu'il corresponde au format d'entré de DESeq2
star_counts_dt <- dcast(star_counts_dt, gene ~ UUID, value.var = "unstranded")

# On déplace la colonne des identifiants comme nom des lignes
star_counts_dt <- as.data.frame(star_counts_dt)
rownames(star_counts_dt) <- star_counts_dt$gene
star_counts_dt$gene <- NULL


# Chargement du fichier pour passer des identifiant des échantillon à celui des patients
uuid_to_id <- read.table("data/uuid_to_id.tsv", header = TRUE, sep = '\t')

# On convertit donc le nom des colonnes de notre fichier de comptage en identifiant des patients
colnames(star_counts_dt) <- sapply(colnames(star_counts_dt), function(x) {return (uuid_to_id[uuid_to_id$UUID == x, "ID"])})
rm(uuid_to_id)

# Dans le cas où l'on à plusieurs échantillon pour un seul patient, on choisit de garder le premier
dup_colnames <- duplicated(colnames(star_counts_dt))
star_counts_dt <- star_counts_dt[!dup_colnames]
rm(dup_colnames)


# Chargement des métadonnées (fichier _cna.txt)
# Il faut avoir téléchargé les études avec download.R
meta <- read.table("data/studies_data/luad_tcga/luad_tcga/data_CNA.txt", header = TRUE, sep = '\t')

# On récupère le statut de notre gène d'intéret (EIF2AK4) pour chaque patient
meta <- meta[meta$Hugo_Symbol == "EIF2AK4", c(-1, -2)]

# On transpose le tableau et on ajoute les colonne du sexe des patient et de l'identifiant
meta <- rbind(colnames(meta), meta)
meta <- transpose(meta)
meta <- cbind(meta, Sex=NA)
colnames(meta) = c("ID", "EIF2AK4_status","Sex")

# On transforme les identifiant pour qu'il corresponde a ceux du tableau de comptage
for (the_ID in meta$ID) {
  temp <- substr(the_ID, 1, 12)
  temp <- chartr(".", "-", temp)
  meta[meta$ID == the_ID, "ID"] <- temp
}
rm(temp, the_ID)


# Chargement des données cliniques (fichier patient.txt)
# Il faut avoir téléchargé les études avec download.R
# source("scripts/r_scripts/download.R")
# download_TCGA_Cbio(download_all = FALSE)
# rm(download_TCGA_Cbio, download_TCGA_Cbio_all, download_TCGA_Cbio_interest)

clinical = as.data.frame(fread("data/studies_data/luad_tcga/luad_tcga/data_bcr_clinical_data_patient.txt", header = TRUE, sep = '\t'))

# On récupère le sexe des patients et on le place dans les métadonnées
# ATTENTION AVEC COAD c'est "SEX" et non "Sex"
# ET $PATIENT_ID et non ["Patient Identifier"]
# Version alternative en dessous

for (the_ID in meta$ID) {
  if (!(identical(toupper(clinical[clinical["Patient Identifier"] == the_ID, "Sex"]), character(0)))) {
    meta[meta$ID == the_ID, "Sex"] <- toupper(clinical[clinical["Patient Identifier"] == the_ID, "Sex"]) 
  }
}
rm(clinical, the_ID)

# for (the_ID in meta$ID) {
#   if (!(identical(toupper(clinical[clinical == the_ID, "SEX"]), character(0)))) {
#     meta[meta$ID == the_ID, "Sex"] <- toupper(clinical[clinical$PATIENT_ID == the_ID, "SEX"])
#   }
# }
# rm(clinical, the_ID)

# On supprime les lignes des patient pour lesquels il y a des données manquantes
meta <- na.omit(meta)

# On donne un nom plus explicit au statut d'EIF2AK4
meta$EIF2AK4_status <- factor(meta$EIF2AK4_status,
                              levels = c(-2, -1, 0, 1),
                              labels = c("Deep_deletion", "Shallow_deletion", "Normal_deletion", "Amplification"))
meta$EIF2AK4_status <- as.character(meta$EIF2AK4_status)


# On supprime les patients avec un statut d'amplification
meta = meta[meta$EIF2AK4_status != "Amplification",]

# On fait l'intersection entre les patients des métadonnées et des données de comptage
intersected <- intersect(meta$ID, colnames(star_counts_dt))

# On enlève les individus non-présents dans l'intersection pour les données de comptage
for (the_ID in colnames(star_counts_dt)) {
  if (!(the_ID %in% intersected)) {
    star_counts_dt <- star_counts_dt %>% select(-one_of(the_ID))
  }
}
rm(the_ID)



# On enlève les individus non-présents dans l'intersection pour les métadonnées
meta <- meta[meta$ID %in% intersected,]
rm(the_ID, intersected)

# On supprime la colonne des identifiant car maintenant inutiles
rownames(meta) <- meta$ID
meta$ID <- NULL

# Pour s'assurer que les deux matrices soit utilisable par DESeq2
ncol(star_counts_dt) == nrow(meta) # Même nombre

all(colnames(star_counts_dt) %in% rownames(meta)) # Même valeurs

all(colnames(star_counts_dt) == rownames(meta)) # Même ordre
star_counts_dt <- star_counts_dt[, rownames(meta)] # Pour corriger les désordre


# On sépare les individus de sexe différent
male_idx <- meta$Sex == "MALE"
female_idx <- meta$Sex == "FEMALE"

male_count_data <- star_counts_dt[, male_idx]
male_col_data <- meta[male_idx, ]
rm(male_idx)

female_count_data <- star_counts_dt[, female_idx]
female_col_data <- meta[female_idx, ]
rm(female_idx, meta)

# On prépare la matrice DESeq pour les deux sexes
dds_male <- DESeqDataSetFromMatrix(countData = male_count_data, colData = male_col_data, design = ~ EIF2AK4_status)
dds_female <- DESeqDataSetFromMatrix(countData = female_count_data, colData = female_col_data, design = ~ EIF2AK4_status)
rm(female_col_data, female_count_data, male_col_data, male_count_data)

# On crée les modèle
dds_male <- DESeq(dds_male)
dds_female <- DESeq(dds_female)

# On crée des modèles qui compare deux à deux les modalites pour les deux sexes
# Avec BRCA, il y a trop peu d'homme donc il faut seulement le faire avec les femmes
res_male_normal_vs_deep <- results(dds_male, contrast = c("EIF2AK4_status", "Normal_deletion", "Deep_deletion"))
res_male_normal_vs_shallow <- results(dds_male, contrast = c("EIF2AK4_status", "Normal_deletion", "Shallow_deletion"))
res_male_shallow_vs_deep <- results(dds_male, contrast = c("EIF2AK4_status", "Shallow_deletion", "Deep_deletion"))

res_female_normal_vs_deep <- results(dds_female, contrast = c("EIF2AK4_status", "Normal_deletion", "Deep_deletion"))
res_female_normal_vs_shallow <- results(dds_female, contrast = c("EIF2AK4_status", "Normal_deletion", "Shallow_deletion"))
res_female_shallow_vs_deep <- results(dds_female, contrast = c("EIF2AK4_status", "Shallow_deletion", "Deep_deletion"))
rm(dds_female, dds_male)

# On crée un tableau et une liste de synthèse
# Comportera pour chaque comparaison, le fold change des gènes, 
# 0 si non-significatif et NA si données indisponible
study_table <- data.frame(matrix(ncol = 0, nrow = nrow(star_counts_dt)))
rownames(study_table) <- rownames(star_counts_dt)
# Comportera pour chaque comparaison, le total de gène, nombre de significatif,
# de non-significatif, nombre de gènes up et down régulation et NA
study_list = list()

# On crée le dossier des résulats pour l'étude 
study_name = "LUAD"
result_path = paste0("result/", study_name)
unlink(result_path, recursive = TRUE)
dir.create(result_path, recursive = TRUE)

# Fonction qui traite les données et placer les informations qui nous intéresse
# dans le tableau et la liste précedents
Deseq_analyse <- function(resultat, tableau, liste, col_name, fdr) {
  comp_col = rep(NA, nrow(tableau))
  tableau <- cbind(tableau, comp_col)
  names(tableau)[ncol(tableau)] <- col_name
  
  temp_list = list(total_gene = nrow(resultat))
  
  rows_with_na <- apply(resultat, 1, function(row) any(is.na(row)))
  data_with_na <- resultat[rows_with_na,]
  NA_number = nrow(data_with_na)
  temp_list$NA_number <- NA_number
  
  resultat <- na.omit(resultat)
  
  non_significatif <- resultat[resultat$padj >= fdr, ]
  non_significatif_number = nrow(non_significatif)
  temp_list$non_significatif_number <- non_significatif_number
  
  for (gene in rownames(non_significatif)) {
    tableau[gene, col_name] <- 0
  }
  
  significatif = resultat[resultat$padj <= fdr, ]
  significatif_number = nrow(significatif)
  temp_list$significatif_number <- significatif_number
  
  up = nrow(significatif[significatif$log2FoldChange > 0, ])
  temp_list$up <- up
  
  down = nrow(significatif[significatif$log2FoldChange < 0, ])  
  temp_list$down <- down
  
  for (gene in rownames(significatif)) {
    tableau[gene, col_name] <- significatif[gene, "log2FoldChange"]
  }
  
  liste[[col_name]] <- temp_list
  return(list(tableau, liste))
}

# Le vecteur des nom de comparaison
comparisons <- c("female_normal_vs_deep", "female_normal_vs_shallow", "female_shallow_vs_deep",
                 "male_normal_vs_deep", "male_normal_vs_shallow", "male_shallow_vs_deep")

# La liste des données modèles
input_data <- list(res_female_normal_vs_deep, res_female_normal_vs_shallow, res_female_shallow_vs_deep,
                   res_male_normal_vs_deep, res_male_normal_vs_shallow, res_male_shallow_vs_deep)

# On boucle dessus tout en créant un volcano plot et MAplot pour chaque comparaison
# Tout en remplissant la liste et le tableau de sortie
for (i in 1:length(comparisons)) {
  message(comparisons[i])
  
  sample_path = paste0(result_path, '/', comparisons[i])
  unlink(sample_path, recursive = TRUE)
  dir.create(sample_path, recursive = TRUE)
  
  result <- Deseq_analyse(input_data[[i]], study_table, study_list, comparisons[i], 0.10)
  study_table <- result[[1]]
  study_list <- result[[2]]
  
  png(paste0(sample_path, '/', comparisons[i], '_volcano.png'))
  x = input_data[[i]]$log2FoldChange
  y = -log10(input_data[[i]]$padj)
  plot(x, y, xlab = "log2 Fold Change", ylab = "-log10 Padj", 
       main = paste0(comparisons[i], "_volcano_plot"), 
       sub = "FDR < 0.10 (=1) et log2 Fold Change > 1 ou < -1", 
       xlim = c(-5, 5), ylim = c(0, 5), pch = 3, 
       col = ifelse((x > 1 | x < -1) & y > 1,'red','black')
  )
  abline(v = -1, col="red")
  abline(v = 1, col="red")
  abline(h = 1, col="red")
  dev.off()
  
  png(paste0(sample_path, '/', comparisons[i], '_ma.png'))
  plotMA(input_data[[i]], 
         sub = "En bleu, les gènes avec un FDR < 0.10",
         main = paste0(comparisons[i], "_ma_plot")
  )
  dev.off()
}
rm(i, x, y, sample_path, input_data, star_counts_dt, result)
rm(comparisons)
rm(res_female_normal_vs_deep, res_female_normal_vs_shallow,
   res_female_shallow_vs_deep, res_male_normal_vs_deep,
   res_male_normal_vs_shallow, res_male_shallow_vs_deep)

# Un prétraitement pour transformer la liste en tableau
study_summary <- data.table::copy(study_list)
setDT(study_summary)
study_summary = as.data.frame(study_summary)
rownames(study_summary) <- c("Total", "NA", "Non-significatif", "Significatif", "UpRegulation", "DownRegulation")
list_columns <- sapply(study_summary, is.list)
study_summary[, list_columns] <- lapply(study_summary[, list_columns], as.character)
rm(study_list, list_columns)

# On écrit les deux tableau dans le répertoire résultat de l'étude
write.table(study_summary, file = paste0(result_path, '/', study_name, "_study_summary.tsv"), sep = '\t')
write.table(study_table, file = paste0(result_path, '/', study_name, "_study_table.tsv"), sep = '\t')
rm(Deseq_analyse, study_summary)

# On fait les intersecte pour les différentes comparaisons et pour les deux sexes indépendamment
male_normal_vs_deep = row.names(study_table)[which(study_table[, "male_normal_vs_deep"] != 0 & !is.na(study_table[, "male_normal_vs_deep"]))]
male_normal_vs_shallow = row.names(study_table)[which(study_table[, "male_normal_vs_shallow"] != 0 & !is.na(study_table[, "male_normal_vs_shallow"]))]
male_shallow_vs_deep = row.names(study_table)[which(study_table[, "male_shallow_vs_deep"] != 0 & !is.na(study_table[, "male_shallow_vs_deep"]))]

female_normal_vs_deep = row.names(study_table)[which(study_table[, "female_normal_vs_deep"] != 0 & !is.na(study_table[, "female_normal_vs_deep"]))]
female_normal_vs_shallow = row.names(study_table)[which(study_table[, "female_normal_vs_shallow"] != 0 & !is.na(study_table[, "female_normal_vs_shallow"]))]
female_shallow_vs_deep = row.names(study_table)[which(study_table[, "female_shallow_vs_deep"] != 0 & !is.na(study_table[, "female_shallow_vs_deep"]))]

# Fonction qui nous produit un diagramme de Venn et le place dans le répertoire résultat de l'étude 
venn_tri_diagramm <- function(normal_vs_deep, normal_vs_shallow, shallow_vs_deep, the_study_name, gender, the_path) {
  intersection_1_2 = length(intersect(normal_vs_deep, shallow_vs_deep))
  intersection_1_3 = length(intersect(normal_vs_deep, normal_vs_shallow))
  intersection_2_3 = length(intersect(shallow_vs_deep, normal_vs_shallow))
  intersection_1_2_3 = length(intersect(intersect(normal_vs_deep, normal_vs_shallow), shallow_vs_deep))
  
  png(paste0(the_path, '/', the_study_name, '_', gender , '_venn_diagram.png'))
  
  # Create the Venn diagram
  venn.plot <- draw.triple.venn(area1 = length(normal_vs_deep),
                                area2 = length(shallow_vs_deep),
                                area3 = length(normal_vs_shallow),
                                n12 = intersection_1_2,
                                n23 = intersection_2_3,
                                n13 = intersection_1_3,
                                n123 = intersection_1_2_3,
                                category = c("NORMAL VS DEEP", "SHALLOW VS DEEP", "NORMAL VS SHALLOW"),
                                fill = c("red", "green", "blue"),
                                alpha = c(0.35, 0.35, 0.35),
                                main = paste0(the_study_name, "_", gender, "_venn_diagram")
  )
  
  grid.draw(venn.plot)
  
  dev.off()
}

# On appelle les fonctions avec les bons paramètres
venn_tri_diagramm(male_normal_vs_deep, male_normal_vs_shallow, male_shallow_vs_deep, study_name, "male", result_path)
venn_tri_diagramm(female_normal_vs_deep, female_normal_vs_shallow, female_shallow_vs_deep, study_name, "female", result_path)
rm(female_normal_vs_deep, female_normal_vs_shallow, female_shallow_vs_deep,
   male_normal_vs_deep, male_normal_vs_shallow, male_shallow_vs_deep)
rm(venn_tri_diagramm, study_table, result_path, study_name)
