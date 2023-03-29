# Comme n'a pas accès à l'ensemble des données pour pouvoir estimer la charge
# mutationnelle, on va à la place faire la somme des mutations détecter par les
# données de séquençage de l'éxome entier (WXS)
# Meilleur que les données Affymétrix car couvrant 1 à 2 % du génome

# Pour le téléchargemnet des études
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

# Les 3 variables suivantes sont à modifier
study_name = "brca_tcga"
mutation_file_name = "data_mutations.txt"
clinical_file_name = "data_clinical_patient.txt"
cna_file = "data_cna.txt"

path = paste0("data/studies_data/", study_name, "/", study_name, "/")

data_df = fread(paste0(path, mutation_file_name), sep = '\t', header = TRUE, select = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode"))
# A enlever si on souhaite le comptage de tout le génome
# locus_15q15_start = 39800001
# locus_15q15_end = 44500000
#  = "15"
# data_df = data_df[data_df$Start_Position >= locus_15q15_start & data_df$Start_Position <= locus_15q15_end & data_df$Chromosome == chromosome,]

clinical_df = fread(paste0(path, clinical_file_name), sep = '\t', header = TRUE)
names(clinical_df)[names(clinical_df) == "Patient Identifier"] <- "Tumor_Sample_Barcode"


maf = read.maf(data_df, clinical_df)

result = tmb(maf, captureSize = 50, logScale = TRUE)

meta <- read.table(paste0(path, cna_file), header = TRUE, sep = '\t')

# On récupère le statut de notre gène d'intéret (EIF2AK4) pour chaque patient
meta <- meta[meta$Hugo_Symbol == "EIF2AK4", c(-1, -2)]

# On transpose le tableau et on ajoute les colonne du sexe des patient et de l'identifiant
meta <- rbind(colnames(meta), meta)
meta <- transpose(meta)
colnames(meta) = c("ID", "EIF2AK4_status")
for (the_ID in meta$ID) {
  temp <- substr(the_ID, 1, 12)
  temp <- chartr(".", "-", temp)
  meta[meta$ID == the_ID, "ID"] <- temp
}
rm(temp, the_ID)
meta$EIF2AK4_status <- factor(meta$EIF2AK4_status,
                              levels = c(-2, -1, 0, 1),
                              labels = c("Deep_deletion", "Shallow_deletion", "Normal_deletion", "Amplification"))
meta$EIF2AK4_status <- as.character(meta$EIF2AK4_status)

for (the_ID in result$Tumor_Sample_Barcode) {
  temp <- substr(the_ID, 1, 12)
  temp <- chartr(".", "-", temp)
  result[result$Tumor_Sample_Barcode == the_ID, "ID"] <- temp
}
result = result[,-1]
result <- cbind(result, EIF2AK4_status=NA)


result$EIF2AK4_status <- as.character(result$EIF2AK4_status)
for (patient_ID in result$ID) {
  if (any(meta$ID == patient_ID)) {
    result[result$ID == patient_ID, "EIF2AK4_status"] <- meta[meta$ID == patient_ID, "EIF2AK4_status"]
  } else {
    result[result$ID == patient_ID, "EIF2AK4_status"] <- NA
  }
}


result <- merge(result, meta[, c("ID", "EIF2AK4_status")], by = "ID")
result <- result[,-5]
colnames(result)[5] <- "EIF2AK4_status"
result = result[result$EIF2AK4_status != "Amplification",]



boxplot(result$total_perMB_log ~ result$EIF2AK4_status, 
        ylab="TMB / 50MB" , xlab="EIF2AK4 status",
        main="Charge mutationnelle de BRCA en fonction du statut EIF2AK4"
)


res.aov <- aov(total_perMB_log ~ EIF2AK4_status, data = result)
TukeyHSD(res.aov)

path_table = paste0("result/mutational_tmb/", study_name, "_mutational_count.tsv")
write.table(result, file = path_table , sep = "\t", row.names = FALSE)

# Si on souhaite lire 
df = read.table(path_table, header = TRUE, sep = '\t')
