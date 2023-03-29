setwd(dir = "~/Bureau/projet/data")
library(ggvenn)

#Fonction d'analyse différentielle
luad_data <- readRDS("all_data_luad.rds")
coad_data <- readRDS("all_coad_data.rds")
brca_data <- readRDS("all_data_brca.rds") 

differential_analysis <- function(df,list_genes){
  df_gene_pval=data.frame() #Initialisation du dataframe des p-val associée à chaque gène
  x=0
  for (genes in list_genes){
    test_gene <- df[df$Hugo_Symbol==genes,] #df de tous les patients du gène "genes" 
    if (length(test_gene$CNA_EIF2AK4=="1") > 100 && length(test_gene$CNA_EIF2AK4=="0") >100 ) #Condition pour que le nombre de patients par gène soit supérieur à 100 dans les deux groupes
    {
      test_gene_01 <- test_gene[test_gene$CNA_EIF2AK4=="1",] #Groupe EIF2AK4 délétion
      test_gene_00 <- test_gene[test_gene$CNA_EIF2AK4=="0",] #Groupe diploïde 
      tp_wcox <- wilcox.test(test_gene_00$MRNA,test_gene_01$MRNA) #Test de wilcoxon
      x=tp_wcox$p.value #On récupère la p-val
      if (is.na(x)){x=0.7} #Si le test n'a pas marché ou n'a pas de p-val pour une quelconque raison on prend une p-val non significative
      if (mean(test_gene_01$MRNA) <= mean(test_gene_00$MRNA)){level = "down"} #On a parfois plusieurs valeurs pour un même combo de patient et gène, on prend donc la moyenne dans ces cas là.
      else {level = "up"}
      df_gene_pval <- rbind(df_gene_pval,c(genes,x,level)) #On incrémente le dataframe avec la p-val 
      remove(test_gene) #On supprime les df de l'environnement 
      remove(test_gene_01)
      remove(test_gene_00)
    }
  }
  names(df_gene_pval) = c("gene", "pval" , "type")
  return (df_gene_pval)
}

liste_genes_coad <- unique(coad_data$Hugo_Symbol)
liste_genes_luad <- unique(luad_data$Hugo_Symbol)
liste_genes_brca <- unique(brca_data$Hugo_Symbol)

coad_pval <- differential_analysis(coad_data,liste_genes_coad)
luad_pval <- differential_analysis(luad_data,liste_genes_luad)
brca_pval <- differential_analysis(brca_data,liste_genes_brca)

#Il est important de sauvegarder les données de p-val car l'analyse différentielle prend plusieurs dizaines de minutes. 
#Ici, on charge des données déjà analysées. 
coad_pval <- readRDS("coad_pval.rds")
luad_pval <- readRDS("luad_pval.rds")
brca_pval <- readRDS("brca_pval.rds")

adjusted_coad_pval <- coad_pval
adjusted_coad_pval$pval <- p.adjust(coad_pval$pval,method = "fdr")
adjusted_brca_pval <- brca_pval
adjusted_brca_pval$pval <- p.adjust(brca_pval$pval,method = "fdr")
adjusted_luad_pval <- luad_pval
adjusted_luad_pval$pval <- p.adjust(luad_pval$pval,method = "fdr")

#Liste des gènes différentiellement exprimée avec une p-val<0.05

coad_0.05 <- adjusted_coad_pval[adjusted_coad_pval$pval<0.05,]
luad_0.05 <- adjusted_luad_pval[adjusted_luad_pval$pval<0.05,]
brca_0.05 <- adjusted_brca_pval[adjusted_brca_pval$pval<0.05,]

#Liste des gènes différentiellement exprimée avec une p-val<0.01

coad_0.01 <- adjusted_coad_pval[adjusted_coad_pval$pval<0.01,]
luad_0.01 <- adjusted_luad_pval[adjusted_luad_pval$pval<0.01,]
brca_0.01 <- adjusted_brca_pval[adjusted_brca_pval$pval<0.01,]



#Data pour le graph de venn

data_venn_0.05 <- list(coad=coad_0.05$gene,luad=luad_0.05$gene,brca=brca_0.05$gene)

#Liste des gènes significatifs communs aux 3 cancers

genes_3_cancers <- intersect(intersect(data_venn_0.05$coad,data_venn_0.05$luad),data_venn_0.05$brca)

#Caractéristiques patients pour chaque gène.
intersected_luad <- subset(luad_0.05, luad_0.05$gene %in% genes_3_cancers)
intersected_coad <- subset(coad_0.05, coad_0.05$gene %in% genes_3_cancers)
intersected_brca <- subset(brca_0.05, brca_0.05$gene %in% genes_3_cancers)

ggvenn(data_venn_0.05)

