library(tidyverse)

mrna <- read.table("cptac_mRNA_FPKM.tsv",
                   sep = "\t", header=T)

colnames(mrna)[2:4] <- c("SAMPLE_ID", "MLH1_RNA", "MSH2_RNA")

protein <- read.table("cptac_protein.tsv", header=T, sep='\t')
colnames(protein)[2:4] <- c("SAMPLE_ID", "MLH1_prot", "MSH2_prot")

cna <- read.table("cptac_cna.tsv", sep='\t', header=T)
colnames(cna)[2:4] <- c("SAMPLE_ID", "MLH1_cna", "MSH2_cna")


mutations <- read.table("cptac_mutations.tsv", header=T, sep='\t')
colnames(mutations)[2:4] <- c("SAMPLE_ID", "MLH1_mut", "MSH2_mut")

cptac <- full_join(mrna, protein)
cptac <- full_join(cptac, cna)
cptac <- full_join(cptac, mutations)


rm(mrna, protein, cna, mutations)
cptac$STUDY_ID <- NULL

cptac$MLH1_cna <- as.numeric(cptac$MLH1_cna)
cptac$MSH2_cna <- as.numeric(cptac$MSH2_cna)

cptac <- drop_na(cptac)

summary(cptac)

cptac$MLH1_RNA_low <- (cptac$MLH1_RNA <= quantile(cptac$MLH1_RNA, probs=0.12))
cptac$MLH1_prot_low <- (cptac$MLH1_prot <= quantile(cptac$MLH1_prot, probs=0.12))
cptac$MLH1_cna_low <- (cptac$MLH1_cna < 0)
cptac$MLH1_mut_bool <- !(cptac$MLH1_mut == "WT" | cptac$MLH1_mut == "NP")
cptac$MLH1_dna_low <- (cptac$MLH1_cna_low | cptac$MLH1_mut_bool)

cptac$MSH2_RNA_low <- (cptac$MSH2_RNA <= quantile(cptac$MSH2_RNA, probs=0.08))
cptac$MSH2_prot_low <- (cptac$MSH2_prot <= quantile(cptac$MSH2_prot, probs=0.08))
cptac$MSH2_cna_low <- (cptac$MSH2_cna < 0)
cptac$MSH2_mut_bool <- !(cptac$MSH2_mut == "WT" | cptac$MSH2_mut == "NP")
cptac$MSH2_dna_low <- (cptac$MSH2_cna_low | cptac$MSH2_mut_bool)

summary(cptac)

write.csv(cptac, "cptac_MMR_perturbance.csv", row.names=F)

cptac$SAMPLE_ID[cptac$MLH1_dna_low | cptac$MLH1_RNA_low | cptac$MLH1_prot_low]
