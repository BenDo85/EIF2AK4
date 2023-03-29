# Le .gff3 est téléchargé depuis https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/
# Choisir le bon chromosome
# Exécuter depuis le répertoire du projet (EIF2AK4)
# Sur linux : python3 python_scripts/closest_genes.py -g EIF2AK4 -c 15
# Sur windows : python.exe python_scripts\closest_genes.py -g EIF2AK4 -c 15

from argparse import ArgumentParser
from os import path
from urllib import request
import gzip


parser = ArgumentParser()

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


def get_interval_from_gene(genes_list: list, gene_name: str) -> tuple:
    start = None
    end = None
    for gene in genes_list:
        if gene[9] == "Name=" + gene_name:
            start = int(gene[3])
            end = int(gene[4])
            break

    if start == None or end == None:
        raise Exception("Le gène de référence n'est pas trouvé")

    return (start, end)


def extract_gene_from_interval(genes_list: list, lower: int, upper: int) -> list:
    matching_genes = []
    for gene in genes_list:
        start = int(gene[3])
        end = int(gene[4])
        distance = 0
        if end <= lower:
            distance = lower - end
        else:
            distance = start - upper
        if distance < 0:
            distance = 0
        if distance < 222911:
            pass
        temp = gene
        temp.append(distance)
        matching_genes.append(temp)
    return matching_genes


genes_list = extract_gff3_genes(local_file)

interval = get_interval_from_gene(genes_list, args.gene)

matching_genes = extract_gene_from_interval(genes_list, interval[0], interval[1])

matching_genes.sort(key=lambda x: x[-1])

with open("data/result_closest_genes.tsv", "w+") as file:
    file.write("Nom\tID\tDistance\n")
    for gene in matching_genes:
        gene_name = gene[9][gene[9].find("=") + 1 :]
        file.write(
            gene_name
            + "\t"
            + gene[8][gene[8].find(":") + 1 :]
            + "\t"
            + str(gene[-1])
            + "\n"
        )

with open("data/result_closest_genes.tsv", "r") as file:
    print(file.read())
