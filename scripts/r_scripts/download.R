# Il faut que le répertoire de travail soit celui du projet (EIF2AK4)

# Pour le téléchargement des données TGCA par cBioPortal

download_TCGA_Cbio_all <- function(cbio) {
  studies <- getStudies(cbio, buildReport = TRUE)
  
  studies <- data.frame(studies$name, studies$studyId)
  
  firehose_studies <- studies[grepl("TCGA, Firehose Legacy", studies$studies.name), ]
  
  for (studyId in firehose_studies$studies.studyId) {
    if (!(file.exists(paste0("data/studies_data/", studyId)))) {
      path <- downloadStudy(studyId)
      path <- untarStudy(path)
      file.rename(path, paste0("data/studies_data/", studyId))
    }
  }
}

download_TCGA_Cbio_interest <- function(cbio) {
  for (studyId in c('brca_tcga', 'luad_tcga', 'coadread_tcga')) {
    if (!(file.exists(paste0("data/studies_data/", studyId)))) {
      path <- downloadStudy(studyId)
      path <- untarStudy(path)
      file.rename(path, paste0("data/studies_data/", studyId))
    }
  }
}

download_TCGA_Cbio <- function(download_all) {
  stopifnot(is.logical(download_all))
  
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  if (!require("cBioPortalData", quietly = TRUE)) {
    BiocManager::install("cBioPortalData")
  }
  
  library(cBioPortalData)
  
  cbio <- cBioPortal()
  
  if (download_all == TRUE) {
    download_TCGA_Cbio_all(cbio)
  } else {
    download_TCGA_Cbio_interest(cbio)
  }
  
  print("Les fichiers ont été extraits dans data/studies_data")
}


