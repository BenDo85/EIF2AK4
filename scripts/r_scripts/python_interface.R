# Ce script fait l'interface avec les fonctions pythons ayant pour but d'extraire 
# des gènes spécifiques (surtout EIF2AK4)
# Il faut que le répertoire de travail soit celui du projet (EIF2AK4)
# Python3 doit aussi être installé
# Télécharger le .gff3 du chromosome souhaité depuis 
# https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/
# puis le placer dans le répertoire data

if(.Platform$OS.type == "windows") {
  system('python scripts\\python_scripts\\closest_genes.py -g EIF2AK4 -c 15')
} else {
  system('python3 scripts/python_scripts/closest_genes.py -g EIF2AK4 -c 15')
}

# La liste triée des gènes les plus proches d'EIF2AK4
df <- read.table("data/result_closest_genes.tsv", header = TRUE, sep = '\t')

# Pour récuper les 10 gènes les plus proches (en excluant EIF2AK4 lui-même)
df[2:11,]

# Commentaires :
# On récupère ce qui semble être un doublon : PAK6 et BUB1B-PAK6.
# En réalite BUB1B-PAK6 est une gènes "readthrough".
# Il est lié à un processus de traduction particulier qui fait que la protéine
# synthétiser ne s'arrête pas au codon stop habituel (contrairement à PAK6).
# A noter que ces gènes sont étudier pour le cancer.
# source : https://pubmed.ncbi.nlm.nih.gov/19951906/


if(.Platform$OS.type == "windows") {
  system('python.exe scripts\\python_scripts\\all_gene_in_interval.py -c 15 --start 39800001 --stop 44500000 -s')
} else {
  system('python3 scripts/python_scripts/all_gene_in_interval.py -c 15 --start 39800001 --stop 44500000 -s')
}

# La liste des gènes complétement inclus dans le locus 15q15
df2 <- read.table("data/result_gene_in_interval.tsv", header = TRUE, sep = '\t')
df2
