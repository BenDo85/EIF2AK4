# Le .gff3 est téléchargé depuis https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/
# Choisir le bon chromosome
# Exécuter depuis le répertoire du projet (EIF2AK4)
# Sur linux : python3 python_scripts/get_neighbors_genes.py -c 15 -g EIF2AK4 -i 2000000
# Sur windows : python.exe python_scripts\get_neighbors_genes.py -c 15 -g EIF2AK4 -i 2000000

from argparse import ArgumentParser
from os import path
from urllib import request
import gzip


parser = ArgumentParser()

parser.add_argument(
    "-i",
    "--intervalle",
    type=int,
    help="l'intervalle considérée autour du gène",
    required=True,
)
parser.add_argument(
    "-g",
    "--gene",
    type=str,
    help="le nom du gène pour lequel on va chercher ses voisins (par défaut=EIF2AK4)",
    required=True,
)
parser.add_argument(
    "-c",
    "--chromosome",
    type=str,
    help="Nom du chromosome (par défaut=15)",
    required=True,
)
parser.add_argument(
    "-s",
    "--strict",
    action="store_true",
    help="Exclure les gènes qui ne sont pas complètement inclus dans l'intervalle (par défaut=False)",
    default=False,
)

args = parser.parse_args()


local_file = "data/" + "Homo_sapiens.GRCh38.109.chromosome." + "15" + ".gff3.gz"

if path.exists(local_file) == False:
    remote_url = (
        "https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/"
        + "Homo_sapiens.GRCh38.109.chromosome."
        + args.chromosome
        + ".gff3.gz"
    )
    request.urlretrieve(remote_url, local_file)


def extract_gff3_genes(path: str) -> list:
    genes = []
    file = gzip.open(path, "rt")

    for line in file.readlines():
        temp_line = line.split("\t")
        if len(temp_line) > 1 and temp_line[2] == "gene":
            temp = temp_line[8].split(";")
            temp_line.pop(8)
            temp_line.extend(temp)
            if "Name=" in temp_line[9]:
                genes.append(temp_line)

    file.close()
    return genes


def get_interval_from_gene(genes_list: list, gene_name: str, interval: int) -> tuple:
    start = None
    end = None
    for gene in genes_list:
        if gene[9] == "Name=" + gene_name:
            start = int(gene[3]) - interval
            end = int(gene[4]) + interval
            break

    if start == None or end == None:
        raise Exception("Le gène de référence n'est pas trouvé")

    return (start, end)


def extract_gene_from_interval(
    genes_list: list, lower: int, upper: int, strict: bool
) -> list:
    matching_genes = []
    if strict:
        for gene in genes_list:
            start = int(gene[3])
            end = int(gene[4])
            if start >= lower and end <= upper:
                matching_genes.append(gene)
    else:
        for gene in genes_list:
            start = int(gene[3])
            end = int(gene[4])
            if (
                (start <= lower and end >= lower)
                or (start <= upper and end >= upper)
                or (start >= lower and end <= upper)
            ):
                matching_genes.append(gene)
    return matching_genes


genes_list = extract_gff3_genes(local_file)

interval = get_interval_from_gene(genes_list, args.gene, args.intervalle)

matching_genes = extract_gene_from_interval(
    genes_list, interval[0], interval[1], args.strict
)

print(matching_genes[0][8][matching_genes[0][8].find(":") + 1 :])

with open("data/result_neighbors_genes.tsv", "w+") as file:
    file.write("Name" + "\t" + "ID\n")
    for gene in matching_genes:
        file.write(
            gene[9][gene[9].find("=") + 1 :]
            + "\t"
            + gene[8][gene[8].find(":") + 1 :]
            + "\n"
        )

with open("data/result_neighbors_genes.tsv", "r") as file:
    print(file.read())
    print(f"Au total : {len(matching_genes)} gènes récupérés")
