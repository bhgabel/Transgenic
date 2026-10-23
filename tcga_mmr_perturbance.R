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


#Mutatations in BRCA1 / BRCA2
brca <- read.table("tcga_brca_mutations.tsv", header=T, sep='\t')

brca[1] <- NULL
for(i in 1:length(brca[[1]])){
  if(brca[i, 2] == "NP"){
    brca[i, 2] <- NA
  }
  else if(brca[i, 3] == "NP"){
    brca[i, 3] <- NA
  }
}
brca <- drop_na(brca)

summary(brca)
gdc <- left_join(gdc, brca, by=c("Sample.ID"="SAMPLE_ID"))

cna <- read.table("tcga_brca_cna.tsv", sep='\t', header=T)
cna[1] <- NULL

cna$BRCA1_cna_loss <- (cna$BRCA1 < 0)
cna$BRCA2_cna_loss <- (cna$BRCA2 < 0)

colnames(cna)[2:3] <- c("BRCA1_cna", "BRCA2_cna")

brca <- left_join(brca, cna)

brca$combo <- factor(NA, levels=c("None", "Mutation", "CNA", "Both"))
for(i in 1:length(brca[[1]])){
  #both false
  if(!(brca[i, 4] | brca[i, 5]) & !(brca[i, 8] | brca[i, 9])){
    brca$combo[[i]] <- "None"
  }
  #both true
  else if((brca[i, 4] | brca[i, 5]) & (brca[i, 8] | brca[i, 9])){
    brca$combo[[i]] <- "Both"
  }
  #only mutation
  else if((brca[i, 4] | brca[i, 5]) & !(brca[i, 8] | brca[i, 9])){
    brca$combo[[i]] <- "Mutation"
  }
  #only cna
  else if(!(brca[i, 4] | brca[i, 5]) & (brca[i, 8] | brca[i, 9])){
    brca$combo[[i]] <- "CNA"
  }
  else{
    brca$combo <- NA
  }
}

table(gdc$PAM50, gdc$BRCA1_mut | gdc$BRCA2_mut)
table(gdc$MLH1_low, gdc$BRCA1_mut | gdc$BRCA2_mut)

ggplot(data=as.data.frame(table(gdc$PAM50, gdc$BRCA1_mut | gdc$BRCA2_mut)),
       aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity", position="fill", color="black") +
  geom_text(stat="identity", aes(label=Freq),
            position=position_fill(vjust=0.1)) +
  theme_classic() +
  scale_y_continuous(labels=c("0%", "10%", "20%", "30%")) +
  coord_cartesian(ylim=c(0, 0.3)) +
  labs(title="BRCA Mutations Across PAM50 Subtypes (TCGA)", x=NULL, y="Proportion",
       fill="BRCA Mutation") +
  scale_fill_discrete(labels=c("FALSE"="Absent", "TRUE"="Present"))
  
