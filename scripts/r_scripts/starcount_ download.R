if (!require("BiocManager")) {
  install.packages("BiocManager")
}
library("BiocManager")

if (!require("GenomicDataCommons")) {
  BiocManager::install('GenomicDataCommons')
}
library("GenomicDataCommons")

status()

ahi = read.table("data/uuid_to_id.tsv", header = TRUE, sep = '\t')

star_count_df = data.frame(matrix(ncol = 11, nrow = 0))
colnames(star_count_df) <- c("gene_id", "gene_name", "gene_type", "unstranded", "stranded_first", "stranded_second", "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded", "UUID", "ID")

project = 'TCGA-BRCA'

for (project in c('TCGA-BRCA')) {
  ge_manifest <- files() |>
    filter( cases.project.project_id == project ) |>
    filter( type == 'gene_expression' ) |>
    filter( access == 'open') |>
    manifest()
  
  destdir <- paste0(getwd(), "/data/starCount/", project, "/DL")
  
  dir.create(destdir, recursive = TRUE)
  
  for (uuid in ge_manifest$id) {
    path = gdcdata(uuid)
    file.copy(paste0(path), destdir)
    
    file_name = paste0(destdir, '/', ge_manifest[ge_manifest$id == uuid, "file_name"])
    
    headers = read.table(file_name, skip = 0, header = F, nrows = 1, as.is = T)
    headers = unlist(c(headers, "UUID", "ID"))
    
    df = read.table(file_name, skip = 6, header = F, sep = "\t", quote = "")
    df <- cbind(df, rep(uuid, nrow(df)))
    df <- cbind(df, rep(ahi[ahi$UUID == uuid, "ID"], nrow(df)))
    
    colnames(df) = headers
    
    star_count_df <- rbind(df, star_count_df)
  }
}

write.table(star_count_df, file = "star_count_BRCA.tsv", sep = '\t', row.names = FALSE)
