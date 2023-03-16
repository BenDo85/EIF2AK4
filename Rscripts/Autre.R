if(.Platform$OS.type == "windows") {
  system('python python_scripts\\closest_genes.py -g EIF2AK4 -c 15')
} else {
  system('python3 python_scripts/closest_genes.py -g EIF2AK4 -c 15')
}

# La liste triée des gènes les plus proches d'EIF2AK4
df = read.table("python_scripts/result.tsv", header = TRUE, sep = '\t')

