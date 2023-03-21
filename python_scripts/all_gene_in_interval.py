# Télécharger le .gff3 depuis https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/
# Choisir le bon chromosome
# Exécuter depuis le répertoire parent
# Les données du .gff3 doivent se trouver dans le répertoire data
# Le début et la fin correspond à la position du locus 15q15
# Sur linux : python3 python_scripts/closest_genes.py -c 15 --start 39800001 --stop 44500000
# Sur windows : python.exe python_scripts\all_gene_in_interval.py -c 15 --start 39800001 --stop 44500000

from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    "--start",
    type=int,
    help="Début de l'intervalle",
    required=True,
)
parser.add_argument(
    "--stop",
    type=int,
    help="Fin de l'intervalle",
    required=True,
)
parser.add_argument(
    "-c",
    "--chromosome",
    type=str,
    help="Nom du chromosome",
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


def extract_gff3_genes(path: str) -> list:
    genes = []
    file = open(path, "r")

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


genes_list = extract_gff3_genes(
    "data/Homo_sapiens.GRCh38.109.chromosome." + args.chromosome + ".gff3"
)

matching_genes = extract_gene_from_interval(
    genes_list, args.start, args.stop, args.strict
)

with open("data/result_gene_in_interval.tsv", "w+") as file:
    file.write("Nom\tID\tStart\tEnd\n")
    for gene in matching_genes:
        gene_name = gene[9][gene[9].find("=") + 1 :]
        file.write(
            gene_name
            + "\t"
            + gene[8][gene[8].find(":") + 1 :]
            + "\t"
            + gene[3]
            + "\t"
            + gene[4]
            + "\n"
        )

with open("data/result_gene_in_interval.tsv", "r") as file:
    print(file.read())