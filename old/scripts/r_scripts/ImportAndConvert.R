# Après certaines étapes, les objets qui ne sont plus nécessaires sont supprimés
# Cela permet de garder de la place en mémoire et de mieux se répérer dans l'environnement de travail

# Installations si nécessaire
if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

if("GenomicDataCommons" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("GenomicDataCommons")
}

if("BiocParallel" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("GenomicDataCommons")
} # Je sais pas si ce package fonctionne vraiment



# Activations des bibliothèques
library("GenomicDataCommons")
library("BiocParallel")
library("dplyr")


# Pour vérifier l'accès à GDC
status()


# Activer le téléchargement en multithreading
register(MulticoreParam())


# Filtrage des recherches
ge_manifest <- files() |>
  filter( cases.project.project_id == 'TCGA-LUAD') |>
  filter( type == 'clinical_supplement' ) |>
  filter( data_format == 'bcr biotab' ) |>
  manifest()


# Préciser le répertoire de destination
destdir <- "/home/benjamin/Bureau/aaaaa"


# Pour récupérer les ID
ge_manifest = ge_manifest$id[grep(".txt$", ge_manifest$filename)]


# Lancer le téléchargement
fnames <- lapply(head(ge_manifest), gdcdata)
rm(ge_manifest)


# Déplace les fichiers téléchargés du répertoire temporaire vers celui destinataire
file.copy(paste0(fnames), destdir)
file.remove(paste0(fnames))
rm(fnames)


# Récupère l'ensemble des fichiers
# Les convertis en dataframe
# Saute les lignes bizarres (2 et 3 en particulier)
# Remplace les valeurs inconnus/invalides par NA
# Aggrège les valeurs répété en un vecteur (un peu comme GROUP BY en SQL)
# Place l'ensemble dans une liste de dataframe
li = list()
for (file in list.files(destdir)) {
  filename = paste(destdir, '/', file, sep = '')
  headers = read.table(filename, skip = 0, header = F, nrows = 1, as.is = T)
  df = read.table(filename, skip = 3, header = F, sep = "\t", quote = "")
  colnames(df) = headers
  df[df=="[Not Available]" | df=="[Not Evaluated]" | df=="[Unknown]" | df=="[Not Applicable]"] <- NA
  df = aggregate(
    by = df[colnames(df)[1]],
    x = df[colnames(df)[-c(1)]],
    FUN = "c"
  )
  li = append(li, list(df))
}
rm(df, headers, file, filename)


# Récupère le dataframe avec le plus de lignes qui correspond donc à celui avec tout les patients
# le place dans la valeur main_frame
# Le supprime de la liste
i = 0
best_index = 0
max_rows = 0
for (df in li) {
  i = i + 1
  if (nrow(df) > max_rows) {
    best_index = i
    max_rows = nrow(df)
  }
}
main_frame = data.frame(li[best_index])
li = li[-best_index];
rm(df, i, max_rows, best_index)


# Exporte les colonnes exclusives aux dataframes de la liste vers le dataframe principal
# On initialise toutes les valeurs d'une nouvelle colonne par non-disponible (NA)
# On convertit les valeurs groupés en string (expression compliqué avec paste sur la ligne 109)
# Elles sont regroupées selon l'UUID du patient
for (df in li) {
  for (col_name in colnames(data.frame(df))) {
    if (!(col_name %in% colnames(main_frame))) {
      main_frame[col_name] <- NA
      for (patient in df[[colnames(df)[1]]]) {
        main_frame[main_frame[colnames(main_frame)[1]] == patient, col_name] = paste(df[df[colnames(df)[1]] == patient, col_name], collapse = ',')
      }
    }
  }
}
rm(df, col_name, patient, li)


# La fin de l'importation est en cours de finalisation


# Utiliser ça pour convertir les valeurs groupé (par exemple : "(\"Clinial Progressive Disease\", \"Clinial Progressive Disease\")") en vecteur
aya = eval(str2lang(string))


# Exporter au format TSV dans le répertoire de travail
write.table(main_frame, file='test.tsv', sep='\t', row.names = F)
