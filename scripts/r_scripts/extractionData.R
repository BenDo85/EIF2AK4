if (!require("BiocManager")) {
  install.packages("BiocManager")
}
library("BiocManager")

if (!require("GenomicDataCommons")) {
  BiocManager::install('GenomicDataCommons')
}
library("GenomicDataCommons")

extractDatas <- function(listeEtude, label, gene) {
  dataMutations <- NULL
  for (el in listeEtude){
    sample <- allSamples(cbio, el) #sélection de tous les échantillons
    study <- paste(el, label, sep="")
    #mutations <- mutationData(cbio, molecularProfileIds = study, entrezGeneIds = gene, sampleIds = unlist(sample$sampleId))[[1]]
    mutations <- getDataByGenes(cbio, studyId =el, genes = gene, by = "hugoGeneSymbol", molecularProfileIds = study)[[1]]
    if (is.null(dataMutations)) {
      # si dataMutations est NULL, créer une nouvelle frame de données avec les colonnes de mutations
      dataMutations <- mutations
    } else {
      # si dataMutations n'est pas NULL, lier les mutations à la frame de données existante
      common_cols <- intersect(names(mutations), names(dataMutations))
      dataMutations <- dataMutations[, common_cols]
      mutations <- mutations[, common_cols]
      dataMutations <- rbind(dataMutations, mutations)
    }
    
  }
  return(dataMutations)
}

listeEtude <- c("luad_tcga", "brca_tcga", "coadread_tcga")
genes <- listeGene$Nom #la liste de gènes mais marche aussi pour un seul gène
dataMutations <- extractDatas(listeEtude, "_mutations", genes)
dataGistic <- extractDatas(listeEtude, "_gistic", genes)
dataMRNA <- extractDatas(listeEtude, "_mrna", genes)
mols