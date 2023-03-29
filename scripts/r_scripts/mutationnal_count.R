# Comme n'a pas accès à l'ensemble des données pour pouvoir estimer la charge
# mutationnelle, on va à la place faire la somme des mutations détecter par les
# données de séquençage de l'éxome entier (WXS)
# Meilleur que les données Affymétrix car couvrant 1 à 2 % du génome

# Pour le téléchargement des études
source("scripts/r_scripts/download.R")
download_TCGA_Cbio(download_all = FALSE)

if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library("data.table")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("maftools", quietly = TRUE)) {
  BiocManager::install("maftools")
}
library("maftools")

# Les 4 variables suivantes sont à modifier
study_name = "brca_tcga"
mutation_file_name = "data_mutations.txt" # fichier mutation
clinical_file_name = "data_clinical_patient.txt" # fichier des données cliniques
cna_file = "data_cna.txt" # fichier de Copy Number Alteration des gènes

# On récupère le chemin pour l'étude d'intérêt
path = paste0("data/studies_data/", study_name, "/", study_name, "/")

# On récupère le tableau des mutations et on extrait seulement les colonnes nécessaire pas maftools
data_df = fread(paste0(path, mutation_file_name), sep = '\t', header = TRUE, select = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode"))

# On récupère le tableau des données cliniques
clinical_df = fread(paste0(path, clinical_file_name), sep = '\t', header = TRUE)

# On change le nom de la colonne "Patient Identifier" pour quelle corresponde à celle des mutations
names(clinical_df)[names(clinical_df) == "Patient Identifier"] <- "Tumor_Sample_Barcode"

# On crée le crée l'object maf
maf = read.maf(data_df, clinical_df)

# On crée le tableau TMB (Tumor Mutational Burden)
result = tmb(maf, captureSize = 50, logScale = TRUE)

# On récupère le tableau de statut des gène en terme de CNA
meta <- read.table(paste0(path, cna_file), header = TRUE, sep = '\t')

# On extrait le statut de notre gène d'intéret (EIF2AK4) pour chaque patient
meta <- meta[meta$Hugo_Symbol == "EIF2AK4", c(-1, -2)]

# On transpose le tableau et on ajoute les colonne de l'identifiant
meta <- rbind(colnames(meta), meta)
meta <- transpose(meta)
colnames(meta) = c("ID", "EIF2AK4_status")
for (the_ID in meta$ID) {
  temp <- substr(the_ID, 1, 12)
  temp <- chartr(".", "-", temp)
  meta[meta$ID == the_ID, "ID"] <- temp
}
rm(temp, the_ID)

# On donne un nom plus explicit au statut d'EIF2AK4
meta$EIF2AK4_status <- factor(meta$EIF2AK4_status,
                              levels = c(-2, -1, 0, 1),
                              labels = c("Deep_deletion", "Shallow_deletion", "Normal_deletion", "Amplification"))
meta$EIF2AK4_status <- as.character(meta$EIF2AK4_status)

# On transforme les identifiant pour qu'il corresponde a ceux du tableau de mutation
for (the_ID in result$Tumor_Sample_Barcode) {
  temp <- substr(the_ID, 1, 12)
  temp <- chartr(".", "-", temp)
  result[result$Tumor_Sample_Barcode == the_ID, "ID"] <- temp
}
result = result[,-1]
result <- cbind(result, EIF2AK4_status=NA)

# On ajoute le statut d'EIF2AK4 du patient dans les métadonnées
result$EIF2AK4_status <- as.character(result$EIF2AK4_status)
for (patient_ID in result$ID) {
  if (any(meta$ID == patient_ID)) {
    result[result$ID == patient_ID, "EIF2AK4_status"] <- meta[meta$ID == patient_ID, "EIF2AK4_status"]
  } else {
    result[result$ID == patient_ID, "EIF2AK4_status"] <- NA
  }
}

# On supprime les patients avec un statut d'amplification
result = result[result$EIF2AK4_status != "Amplification",]


# On crée le boxplot de l'étude
boxplot(result$total_perMB_log ~ result$EIF2AK4_status, 
        ylab="TMB / 50MB" , xlab="EIF2AK4 status",
        main="Charge mutationnelle de BRCA en fonction du statut EIF2AK4"
)

# On fait une anova et on visualise les p valeur une à une
res.aov <- aov(total_perMB_log ~ EIF2AK4_status, data = result)
TukeyHSD(res.aov)

# On écrit la tableau dans le répertoire résultat
path_table = paste0("result/mutational_tmb/", study_name, "_mutational_count.tsv")
write.table(result, file = path_table , sep = "\t", row.names = FALSE)
