# Ce fichier télécharge les Star Counts de chaque patient pour un ou plusieurs projet
# Ces fichiers sont ensuite fusionné dans un grand final qui facilitera nos analyses

if (!require("BiocManager")) {
  install.packages("BiocManager")
}
library("BiocManager")

if (!require("GenomicDataCommons")) {
  BiocManager::install('GenomicDataCommons')
}
library("GenomicDataCommons")

# Les fichiers sont téléchargable à cette adresse https://drive.google.com/drive/folders/1pGrPft4rKYEvKlXKsfDUkKt6ZOexXwv0?usp=sharing

# Pour vérifier le statut d'accès à l'api de GDC
status()

# Tableau pour permettre la conversion d'un indentifiant d'échantillon vers
# un patient
uuid_to_id = read.table("data/uuid_to_id.tsv", header = TRUE, sep = '\t')

# Attention, consommme beaucoup de RAM
projects = c('LUAD', 'COAD', 'BRCA')

# Pour chaque projet
for (project in projects) {
  
  # On préconstruit le tableau regroupant les Star Counts
  star_count_df = data.frame(matrix(ncol = 11, nrow = 0))
  colnames(star_count_df) <- c("gene_id", "gene_name", "gene_type", "unstranded", "stranded_first", "stranded_second", "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded", "UUID", "ID")
  
  # On récupère les liens vers les fichiers que l'on souhaite télécharger
  ge_manifest <- files() |>
    filter( cases.project.project_id == paste0("TCGA-", project) ) |>
    filter( type == 'gene_expression' ) |>
    filter( access == 'open') |>
    manifest()
  
  # On crée le répertoire de destination
  destdir <- paste0(getwd(), "/data/starCount/TCGA-", project, "/DL")
  dir.create(destdir, recursive = TRUE)
  
  # Pour chaque fichier
  for (uuid in ge_manifest$id) {
    # On les téléchargent
    path = gdcdata(uuid)
    
    # On le copie dans le bon répertoire
    file.copy(paste0(path), destdir)
    file_name = paste0(destdir, '/', ge_manifest[ge_manifest$id == uuid, "filename"])
    
    # Il est traité pour en extraitre les bonnes lignes
    headers = read.table(file_name, skip = 0, header = F, nrows = 1, as.is = T)
    headers = unlist(c(headers, "UUID", "ID"))
    
    # Il est ensuite ajouté au tableau complet
    df = read.table(file_name, skip = 6, header = F, sep = "\t", quote = "")
    df <- cbind(df, rep(uuid, nrow(df)))
    df <- cbind(df, rep(uuid_to_id[uuid_to_id$UUID == uuid, "ID"], nrow(df)))
    
    colnames(df) = headers
    
    star_count_df <- rbind(df, star_count_df)
  }
  # On écrit le tableau final
  write.table(star_count_df, file = paste0("data/star_count_", project, ".tsv"), sep = '\t', row.names = FALSE)
  
}
