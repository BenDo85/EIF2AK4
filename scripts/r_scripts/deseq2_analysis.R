if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install.packages("dplyr")

if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

BiocManager::install("BiocParallel")
library("BiocParallel")
options(MulticoreParam=MulticoreParam(workers=12))

library("data.table")
setDTthreads(12)

library("dplyr")

install.packages("VennDiagram")
library("VennDiagram")

star_counts_dt <- fread("data/star_count/star_count_COAD.tsv", header = TRUE, sep = '\t', data.table = TRUE, select = c("gene_id", "gene_name", "unstranded", "UUID"))
star_counts_dt$gene <- paste(star_counts_dt$gene_id, star_counts_dt$gene_name, sep = "_")

star_counts_dt <- dcast(star_counts_dt, gene ~ UUID, value.var = "unstranded")

star_counts_dt <- as.data.frame(star_counts_dt)

# Set gene_id as rownames and remove the gene_id column
rownames(star_counts_dt) <- star_counts_dt$gene
star_counts_dt$gene <- NULL


uuid_to_id <- read.table("data/uuid_to_id.tsv", header = TRUE, sep = '\t')


colnames(star_counts_dt) <- sapply(colnames(star_counts_dt), function(x) {return (uuid_to_id[uuid_to_id$UUID == x, "ID"])})

# Remove duplicate, choose the first one
dup_colnames <- duplicated(colnames(star_counts_dt))
star_counts_dt <- star_counts_dt[!dup_colnames]
rm(dup_colnames)

meta <- read.table("data/studies_data/coadread_tcga/coadread_tcga/data_cna.txt", header = TRUE, sep = '\t')

meta <- meta[meta$Hugo_Symbol == "EIF2AK4", c(-1, -2)]

meta <- rbind(colnames(meta), meta)

meta <- transpose(meta)
meta <- cbind(meta, Sex=NA)
colnames(meta) = c("ID", "EIF2AK4_status","Sex")

for (the_ID in meta$ID) {
  temp <- substr(the_ID, 1, 12)
  temp <- chartr(".", "-", temp)
  meta[meta$ID == the_ID, "ID"] <- temp
}
rm(temp, the_ID)

clinical = read.table("data/studies_data/coadread_tcga/coadread_tcga/data_clinical_patient.txt", header = TRUE, sep = '\t')

for (the_ID in meta$ID) {
  if (!(identical(toupper(clinical[clinical$PATIENT_ID == the_ID, "SEX"]), character(0)))) {
    meta[meta$ID == the_ID, "Sex"] <- toupper(clinical[clinical$PATIENT_ID == the_ID, "SEX"])
  }
}
rm(clinical, the_ID)
meta <- na.omit(meta)


meta$EIF2AK4_status <- factor(meta$EIF2AK4_status,
                             levels = c(-2, -1, 0, 1),
                             labels = c("Deep_deletion", "Shallow_deletion", "Normal_deletion", "Amplification"))

# Convert back to character if needed
meta$EIF2AK4_status <- as.character(meta$EIF2AK4_status)


# Choose to keep or not the Amplification
meta = meta[meta$EIF2AK4_status != "Amplification" & meta$Sex == "FEMALE",]

intersected <- intersect(meta$ID, colnames(star_counts_dt))

# filter star_counts_dt sur les individus présents dans l'intersect
for (the_ID in colnames(star_counts_dt)) {
  if (!(the_ID %in% intersected)) {
    star_counts_dt <- star_counts_dt %>% select(-one_of(the_ID))
  }
}
rm(the_ID)

# filters meta sur les individus présents dans l'intersect
for (the_ID in meta$ID) {
  if (!(the_ID %in% intersected)) {
    meta <- meta[!meta$ID == the_ID,]
  }
}
rm(the_ID)

rownames(meta) <- meta$ID
meta$ID <- NULL

# Pour s'assurer que les matrices soit utilisable par DESeq2
ncol(star_counts_dt) == nrow(meta) # Même nombre

all(colnames(star_counts_dt) %in% rownames(meta)) # Même valeurs

all(colnames(star_counts_dt) == rownames(meta)) # Même ordre
star_counts_dt <- star_counts_dt[, rownames(meta)]

Sex = meta$Sex
EIF2AK4_status = meta$EIF2AK4_status

dds <- DESeqDataSetFromMatrix(countData = star_counts_dt, colData = meta, design = ~ EIF2AK4_status)

test = DESeq(dds, test="Wald", minReplicatesForReplace = Inf)
results = results(test, contrast=c("EIF2AK4_status", "Normal_deletion", "Deep_deletion"), parallel = TRUE, cooksCutoff = FALSE, independentFiltering = FALSE)

fdr = 0.01

COAD_FEMALE_SHALLOW_DEEP  <- rep(NA, length(row.names(full_df)))
full_df <- cbind(full_df, COAD_FEMALE_SHALLOW_DEEP )

nrow(results) # total number of genes


rows_with_na <- apply(results, 1, function(row) any(is.na(row)))
data_with_na <- results[rows_with_na,]
nrow(data_with_na) # total number of missing value
for (gene in rownames(data_with_na)) {
  full_df[gene, "COAD_FEMALE_SHALLOW_DEEP "] <- NA
}

results <- na.omit(results)

non_significatif <- results[results$padj >= fdr, ]
nrow(results[results$padj >= fdr, ]) # Nombre de non-significatif
for (gene in rownames(non_significatif)) {
  full_df[gene, "COAD_FEMALE_SHALLOW_DEEP "] <- 0
}

significatif = results[results$padj <= fdr, ]
nrow(significatif) # Nombre de significatif total

nrow(significatif[significatif$log2FoldChange > 0, ]) # Nombre de significatif up

nrow(significatif[significatif$log2FoldChange < 0, ]) # Nombre de significatif down


for (gene in rownames(significatif)) {
  full_df[gene, "COAD_FEMALE_SHALLOW_DEEP "] <- significatif[gene, "log2FoldChange"]
}

plot(results$log2FoldChange, -log10(results$padj))
plotMA(results)


COAD_FEMALE_SHALLOW_DEEP = rownames(significatif)
COAD_FEMALE_NORMAL_SHALLOW = rownames(significatif)
COAD_FEMALE_NORMAL_DEEP = rownames(significatif)

# Calculate intersections
intersection_1_2 = length(intersect(COAD_FEMALE_NORMAL_DEEP, COAD_FEMALE_SHALLOW_DEEP))
intersection_1_3 = length(intersect(COAD_FEMALE_NORMAL_DEEP, COAD_FEMALE_NORMAL_SHALLOW))
intersection_2_3 = length(intersect(COAD_FEMALE_SHALLOW_DEEP, COAD_FEMALE_NORMAL_SHALLOW))
intersection_1_2_3 = length(intersect(intersect(COAD_FEMALE_NORMAL_DEEP, COAD_FEMALE_NORMAL_SHALLOW), COAD_FEMALE_SHALLOW_DEEP))

# Create the Venn diagram
venn.plot <- draw.triple.venn(area1 = length(COAD_FEMALE_NORMAL_DEEP),
                              area2 = length(COAD_FEMALE_SHALLOW_DEEP),
                              area3 = length(COAD_FEMALE_NORMAL_SHALLOW),
                              n12 = intersection_1_2,
                              n23 = intersection_2_3,
                              n13 = intersection_1_3,
                              n123 = intersection_1_2_3,
                              category = c("NORMAL VS DEEP", "SHALLOW VS DEEP", "NORMAL VS SHALLOW"),
                              fill = c("red", "green", "blue"),
                              alpha = c(0.35, 0.35, 0.35))

# Display the plot
grid.draw(venn.plot)




write.table(full_df, file = "full_df.tsv", sep = "\t", row.names = TRUE)

full_df = data.frame(matrix(ncol = 0, nrow = nrow(results)))
row.names(full_df) = rownames(results)

full_df = read.table("full_df.tsv", sep = '\t', header = TRUE)
