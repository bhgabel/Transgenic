library(ggbiplot)
library(readxl)
library(tidyverse)
library(genefu)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

#import datasets
mouse <- read_xlsx("~/Documents/Haricharan/Transgenic/TgMiceRNA.xlsx", sheet="zscore")
mouse$Gene.name <- str_to_upper(mouse$Gene.name)
pam50 <- read.csv("~/Documents/Haricharan/PAM50 gene list.csv", header=F)
pam50$V1[16] <- "ORC6"
colnames(pam50) <- "Genes"

#tidy data, observations are each row instead of column
mt <- t(mouse[,3:length(mouse)])
rownames(mt) <- colnames(mouse)[3:9]
colnames(mt) <- mouse$Gene.name

#remove columns without any values
mt <- mt[, colSums(mt) != 0]
mt <- data.frame("sample" = factor("MLT", levels=c("MLT", "MST")), mt)
mt$sample[5:7] <- "MST"

#add pam50 subtypes to see if pca matches
data("pam50.robust")
pam50.pred <- molecular.subtyping(sbt.model = "pam50", data = mt[,-1],
                    annot = pam50,
                    do.mapping = FALSE)
as.data.frame(pam50.pred$subtype)
mt <- add_column(mt, "pam50" = pam50.pred$subtype, .after=1)
mouse.pam50 <- c("Basal", "Normal", "Basal", "LumB", "Her2", "LumA", "LumB")

#subset only coding genes
mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")
coding.genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name",
                                              "transcript_biotype"),
                        filters = c("transcript_biotype"),
                        values = "protein_coding", mart = mart)
mt <- mt[colnames(mt) %in% coding.genes$external_gene_name]

#subset without pam50 genes
mt.pam50 <- dplyr::select(mt, !any_of(pam50$Genes))

#import pt data
tcga <- read.table("~/Documents/Haricharan/Transgenic/brca_tcga_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt", header=T)
tcga <- drop_na(tcga)

#change from Entrez id to gene symbol
genes <- clusterProfiler::bitr(tcga$Entrez_Gene_Id, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)

tcga$Entrez_Gene_Id <- as.character(tcga$Entrez_Gene_Id)
tcga <- left_join(genes, tcga, by=c("ENTREZID" = "Entrez_Gene_Id"))

#assign row names and transpose df
rownames(tcga) <- tcga$SYMBOL
tcga <- as.data.frame(t(tcga[3:length(tcga)]))

#subset for coding genes and intersection of mouse/tcga datasets
tcga <- tcga[colnames(tcga) %in% coding.genes$external_gene_name]
tcga <- tcga[colnames(tcga) %in% colnames(mt)]
mt <- mt[colnames(mt) %in% colnames(tcga)]

#assign pam50 subtypes
pam50.pred <- molecular.subtyping(sbt.model = "pam50", data = tcga,
                                  annot = pam50,
                                  do.mapping = FALSE)

tcga <- add_column(tcga, "pam50" = pam50.pred$subtype, .before=1)

mlh1 <- quantile(tcga$MLH1, 0.12)
msh2 <- quantile(tcga$MSH2, 0.08)

#subset for patients with either low MLH1 or MSH2 (but not both)
tcga.mmr <- subset(tcga, xor(MLH1 <= mlh1,  MSH2 <= msh2))

#not doing this
#sort genes by variance, and subset for 1000 highest
# vars <- apply(mt[,-c(1,2)], 2, var)
# top.genes <- head(order(vars, decreasing=T), 500)
# top.mouse <- colnames(mt)[top.genes + 2]
# top.genes <- top.mouse[top.mouse %in% colnames(tcga)]

#combine tcga and mouse data
mt <- add_column(mt, pam50=mouse.pam50, .before=1)
total <- rbind(mt, tcga)
total <- data.frame(type = factor("Human", levels=c("Mouse", "Human")),
                    pam50 = factor("Normal", levels=c("LumA", "LumB", "Her2",
                                                      "Normal", "Basal")),
                    total)
total$type[1:7] <- "Mouse"
total$pam50[1:7] <- mouse.pam50
total$pam50[8:length(total$pam50)] <- tcga$pam50

#combine MMR loss and mouse data
mmr <- rbind(mt, tcga.mmr)
mmr <- data.frame(type = factor("Human", levels=c("Mouse", "Human")),
                  MMR = factor("MLH1", levels=c("MLH1", "MSH2")), mmr)
mmr$type[1:7] <- "Mouse"


mmr$MMR <- ifelse(mmr$MLH1 <= mlh1, "MLH1", "MSH2")
mmr$MMR[1:4] <- "MLH1"
mmr$MMR[5:7] <- "MSH2"

mmr <- mutate(mmr, MMR = as.factor(MMR), pam50 = as.factor(pam50))
mmr <- add_column(mmr, combo = as.character(mmr$MMR:mmr$pam50), .after=3)
for(i in 1:length(mmr$combo)){
  if(mmr$MMR[[i]] == "MLH1"){
    if(mmr$pam50[[i]] == "LumA" | mmr$pam50[[i]] == "LumB"){
      mmr$combo[[i]] <- "MLH1:Lum"
    }
  }
  else{
    if(mmr$pam50[[i]] == "LumA" | mmr$pam50[[i]] == "LumB"){
      mmr$combo[[i]] <- "MSH2:Lum"
  }
}}

write.csv(mmr, "mouse_pca_prepped.csv")


pc2 <- prcomp(mt[,3:length(mt)], center=T)
pc3 <- prcomp(mt.pam50[,3:length(mt.pam50)], center=T, scale.=T)

#validation clustering
pc.tcga <- prcomp(dplyr::select(tcga, any_of(pam50$Genes)), center=T, scale.=T)

#all patients ± pam50
combo <- prcomp(total[,-c(1,2)], center=T, scale.=T)
notpam <- prcomp(dplyr::select(total[,-c(1,2)], !any_of(pam50$Genes)), center=T, scale.=T)

#MMR low
pc.mmr <- prcomp(mmr[,-c(1:4)], center=T, scale.=T)
pc.mmr.pam50 <- prcomp(dplyr::select(mmr[,-c(1:4)], !any_of(pam50$Genes)), center=T, scale.=T)

pc.mmr <- prcomp(dplyr::select(mmr[,-c(1:4)], any_of(pam50$Genes)), center=T, scale.=T)

#plot PCA's
ggbiplot::ggbiplot(pc3, obs.scale=1, var.scale=1, groups=mt$pam50,
              ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="not pam50")


ggbiplot::ggbiplot(pc.tcga, obs.scale=1, var.scale=1, groups=tcga$pam50,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="TCGA - PAM50 validation clustering")

# all pt ± pam50 ----
ggbiplot::ggbiplot(combo, var.scale=1, groups=total$type,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

ggbiplot::ggbiplot(combo, var.scale=1, groups=total$pam50,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

ggbiplot::ggbiplot(notpam, var.scale=1, groups=total$type,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga - not pam50")

ggbiplot::ggbiplot(notpam, var.scale=1, groups=total$pam50,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

# MMR low ----
ggbiplot::ggbiplot(pc.mmr, var.scale=1, groups=mmr$type,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

ggbiplot::ggbiplot(pc.mmr, var.scale=1, groups=mmr$pam50,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

ggbiplot::ggbiplot(pc.mmr, var.scale=1, groups=mmr$combo,
                   ellipse=F, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="Only PAM50 Genes") +
  geom_point(aes(shape=mmr$type, color=mmr$combo, size=mmr$type)) +
  scale_shape_manual(values=c(3,1)) + scale_size_manual(values=c("Mouse"=5,"Human"=1))

#exclude pam50 genes
ggbiplot::ggbiplot(pc.mmr.pam50, var.scale=1, groups=mmr$type,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga-pam50")

ggbiplot::ggbiplot(pc.mmr.pam50, var.scale=1, groups=mmr$pam50,
                   ellipse=T, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga-pam50")

ggbiplot::ggbiplot(pc.mmr.pam50, var.scale=1, groups=mmr$combo,
                   ellipse=F, circle=F, ellipse.prob=0.68, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="Not PAM50 Genes") +
  geom_point(aes(shape=mmr$type, color=mmr$combo, size=mmr$type)) +
  scale_shape_manual(values=c(3,1)) + scale_size_manual(values=c("Mouse"=5,"Human"=1))
