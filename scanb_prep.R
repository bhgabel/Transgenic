library(tidyverse)

#data downloaded from https://data.mendeley.com/datasets/w8tsvd5zrf/1

clinical <- read.table("SCAN-B/Clinical Metadata/SCANB_clinical_metadata.txt",
                       sep='\t', header=T)

#get gene ids for subset of genes of interest: MLH1, MSH2, BRCA1, and PAM50
genes <- read.table("SCAN-B/Gene Annotations/Gene.ID.ann.txt",
                    sep='\t', header=T)

selection <- read.csv("pam50_gene_list.csv")

gene.ids <- genes$Gene.ID[genes$EntrezGene %in% selection$EntrezGene.ID |
                          genes$Gene.Name %in% c("MLH1", "MSH2", "BRCA1")]



#gene expression table too large to load all at once :)
rna <- read.table("SCAN-B/Gene-expression matrix (protocol-adjusted)/bc4genematrix_noNeg.txt",
                  sep='\t', header=T, nrows=0)

skip <- 0
for(i in 1:14){
  temp <- read.table("SCAN-B/Gene-expression matrix (protocol-adjusted)/bc4genematrix_noNeg.txt",
                      sep='\t', header=T, skip=i*skip, nrows=1000)
  
  temp <- filter(temp, Gene.ID %in% gene.ids)
  rna <- rbind(rna, temp)
}
rna <- read.table("SCAN-B/Gene-expression matrix (protocol-adjusted)/bc4genematrix_noNeg.txt",
                  sep='\t', header=T)

scanb <- 