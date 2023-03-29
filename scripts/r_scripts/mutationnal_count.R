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
  install.packages("maftools")
}
library(maftools)

# Les 3 variables suivantes sont à modifier
study_name = "coadread_tcga"
mutation_file_name = "data_mutations.txt"
clinical_file_name = "data_clinical_patient.txt"

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

path_table = paste0("data/mutational_tmb/patient_mutational_count_", study_name, ".tsv")

write.table(result, file = path_table , sep = "\t", row.names = FALSE)

# Si on souhaite lire 
df = read.table(path_table, header = TRUE, sep = '\t')
