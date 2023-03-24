# Comme n'a pas accès à l'ensemble des données pour pouvoir estimer la charge
# mutationnelle, on va à la place faire la somme des mutations détecter par la
# puce Affymetrix SNP 6.0.

if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

library("data.table")


df = fread("data/data_mutations_extended.txt", sep = '\t', header = TRUE, select = c("Tumor_Sample_Barcode", "Mutation_Status", "Start_Position"))

df = df[df$Mutation_Status == "Somatic", c("Tumor_Sample_Barcode", "Start_Position")]

# A enlever si on souhaite le comptage de tout le génome
locus_15q15_start = 39800001
locus_15q15_end = 44500000
df = df[df$Start_Position >= locus_15q15_start & df$Start_Position <= locus_15q15_end, "Tumor_Sample_Barcode"]

df = aggregate(df, by=list(df$Tumor_Sample_Barcode), FUN=length)

colnames(df) = c("Patient_ID", "Mutational_Count")

write.table(df, file = "data/patient_mutational_count.tsv", sep = "\t", row.names = FALSE)
