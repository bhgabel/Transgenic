library(tidyverse)

mrna <- read_delim("TCGA mRNA expression (RNA Seq V2 RSEM).txt",
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select("SAMPLE_ID", "MLH1", "MSH2")

colnames(mrna) <- c("SAMPLE_ID", "MLH1_RNA", "MSH2_RNA")

protein <- read.csv("tcga_tcpa_rppa.csv", header=T)
colnames(protein)[2:4] <- c("SAMPLE_ID", "MLH1_prot", "MSH2_prot")

cna <- read.table("tcga_cna.tsv", sep='\t', header=T)
colnames(cna)[2:4] <- c("SAMPLE_ID", "MLH1_cna", "MSH2_cna")


mutations <- read.table("tcga_mutations.tsv", header=T)
colnames(mutations)[2:4] <- c("SAMPLE_ID", "MLH1_mut", "MSH2_mut")

tcga <- full_join(mrna, protein)
tcga <- full_join(tcga, cna)
tcga <- full_join(tcga, mutations)


rm(mrna, protein, cna, mutations)
tcga$STUDY_ID <- NULL
tcga$sample_id <- NULL

tcga$MLH1_cna <- as.numeric(tcga$MLH1_cna)
tcga$MSH2_cna <- as.numeric(tcga$MSH2_cna)

tcga <- drop_na(tcga)

summary(tcga)

tcga$MLH1_RNA_low <- (tcga$MLH1_RNA <= quantile(tcga$MLH1_RNA, probs=0.12))
tcga$MLH1_prot_low <- (tcga$MLH1_prot <= quantile(tcga$MLH1_prot, probs=0.12))
tcga$MLH1_cna_low <- (tcga$MLH1_cna < 0)
tcga$MLH1_mut_bool <- !(tcga$MLH1_mut == "WT" | tcga$MLH1_mut == "NP")
tcga$MLH1_dna_low <- (tcga$MLH1_cna_low | tcga$MLH1_mut_bool)

tcga$MSH2_RNA_low <- (tcga$MSH2_RNA <= quantile(tcga$MSH2_RNA, probs=0.08))
tcga$MSH2_prot_low <- (tcga$MSH2_prot <= quantile(tcga$MSH2_prot, probs=0.08))
tcga$MSH2_cna_low <- (tcga$MSH2_cna < 0)
tcga$MSH2_mut_bool <- !(tcga$MSH2_mut == "WT" | tcga$MSH2_mut == "NP")
tcga$MSH2_dna_low <- (tcga$MSH2_cna_low | tcga$MSH2_mut_bool)

summary(tcga)

write.csv(tcga, "TCGA_MMR_perturbance.csv", row.names=F)

tcga$SAMPLE_ID[tcga$MLH1_dna_low | tcga$MLH1_RNA_low | tcga$MLH1_prot_low]
