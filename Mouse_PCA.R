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

#subtypes provided by Megha
mouse.pam50 <- c("Basal", "Normal", "Basal", "LumB", "Her2", "LumA", "LumB")

#subset only coding genes
mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")
coding.genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name",
                                              "transcript_biotype"),
                        filters = c("transcript_biotype"),
                        values = "protein_coding", mart = mart)
mt <- mt[colnames(mt) %in% coding.genes$external_gene_name]

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
data("pam50.robust")
pam50.pred <- molecular.subtyping(sbt.model = "pam50", data = tcga,
                                  annot = pam50,
                                  do.mapping = FALSE)

tcga <- add_column(tcga, "pam50" = pam50.pred$subtype, .before=1)

mlh1 <- quantile(tcga$MLH1, 0.12)
msh2 <- quantile(tcga$MSH2, 0.08)

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
                    MMR = factor("Rest", levels=c("MLH1", "MSH2", "Rest")),
                    total)

total$type[1:7] <- "Mouse"

#remove double loss
total <- subset(total, (type == "Mouse") | (MLH1 > mlh1 | MSH2 > msh2))
for(i in 8:length(total$MMR)){
  if(total$MLH1[[i]] <= mlh1){
    total$MMR[[i]] <- "MLH1"
  } else if(total$MSH2[[i]] <= msh2){
    total$MMR[[i]] <- "MSH2"
  } else{
    total$MMR[[i]] <- "Rest"
  }
}
#set mouse factors
total$MMR[1:4] <- "MLH1"
total$MMR[5:7] <- "MSH2"

total <- mutate(total, MMR = as.factor(MMR), pam50 = as.factor(pam50))
total <- tibble::add_column(total, combo = as.character(total$MMR:total$pam50), .after=3)
for(i in 1:length(total$combo)){
  if(total$MMR[[i]] == "MLH1"){
    if(total$pam50[[i]] == "LumA" | total$pam50[[i]] == "LumB"){
      total$combo[[i]] <- "MLH1:Lum"
    }
  }
  else if(total$MMR[[i]] == "MSH2"){
    if(total$pam50[[i]] == "LumA" | total$pam50[[i]] == "LumB"){
      total$combo[[i]] <- "MSH2:Lum"
    }
  }
  else{
    if(total$pam50[[i]] == "LumA" | total$pam50[[i]] == "LumB"){
      total$combo[[i]] <- "Rest:Lum"
    }
  }
}

#combine MMR loss and mouse data
mmr <- subset(total, MMR != "None")


write.csv(total, "mouse_pca_prepped.csv")

#total <- read.csv("mouse_pca_prepped.csv")
#total[1] <- NULL



#all patients ± pam50
pc.total <- prcomp(total[,-c(1:4)], center=T, scale.=T)
pc.total.notpam <- prcomp(dplyr::select(total[,-c(1:4)], !any_of(pam50$Genes)), center=T, scale.=T)
pc.total.pam <- prcomp(dplyr::select(total[,-c(1:4)], any_of(pam50$Genes)), center=T, scale.=T)

#MMR low
pc.mmr <- prcomp(dplyr::select(mmr[,-c(1:4)], !any_of(pam50$Genes)), center=T, scale.=T)

pc.mmr.pam50 <- prcomp(dplyr::select(mmr[,-c(1:4)], any_of(pam50$Genes)), center=T, scale.=T)


#basal pt - mlh1 vs not
basal <- subset(total, pam50 == "Basal" & type != "Mouse")
basal <- add_column(basal, id = row_number(basal$type), .before=1)
pc.basal <- prcomp(dplyr::select(basal[,-c(1:4)], any_of(pam50$Genes)), center=T, scale.=T)
#all pt - mlh1 + basal vs not

#luminal pt - msh2 vs not
luminal <- subset(total, (pam50 == "LumA" | pam50 == "LumB") & type != "Mouse")
pc.luminal <- prcomp(dplyr::select(luminal[,-c(1:4)], any_of(pam50$Genes)), center=T, scale.=T)
#all pt - msh2 + lum vs not

#plot PCA's

#validation clustering
pc.tcga <- prcomp(dplyr::select(tcga, any_of(pam50$Genes)), center=T, scale.=T)
ggbiplot::ggbiplot(pc.tcga, obs.scale=1, var.scale=1, groups=tcga$pam50,
                   ellipse=T, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="TCGA - PAM50 validation clustering")

# all pt ± pam50 ----
ggbiplot::ggbiplot(pc.total, var.scale=1, groups=total$type,
                   ellipse=T, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

ggbiplot::ggbiplot(pc.total, var.scale=1, groups=total$pam50,
                   ellipse=T, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

ggbiplot::ggbiplot(pc.total.notpam, var.scale=1, groups=total$type,
                   ellipse=T, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga - not pam50")

ggbiplot::ggbiplot(pc.total.notpam, var.scale=1, groups=total$pam50,
                   ellipse=T, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="mouse+tcga")

# MMR low ----
ggbiplot::ggbiplot(pc.mmr.pam50, var.scale=1, groups=mmr$combo,
                   ellipse=F, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="Only PAM50 Genes") +
  geom_point(aes(shape=mmr$type, color=mmr$combo, size=mmr$type)) +
  scale_shape_manual(values=c(3,1)) + scale_size_manual(values=c("Mouse"=5,"Human"=1))

#exclude pam50 genes
ggbiplot::ggbiplot(pc.mmr, var.scale=1, groups=mmr$combo,
                   ellipse=F, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="Not PAM50 Genes") +
  geom_point(aes(shape=mmr$type, color=mmr$combo, size=mmr$type, label=row_number(mmr$combo))) +
  scale_shape_manual(values=c(3,1)) + scale_size_manual(values=c("Mouse"=5,"Human"=1))

#basal - mlh1 vs not
ggbiplot::ggbiplot(pc.basal, var.scale=1, groups=basal$MMR,
         ellipse=F, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="Basal Patients - Only PAM50") +
  scale_color_manual(values=c("None"="grey", "MLH1"="red", "MSH2"="blue"), labels=c("None"="Rest"))

#basal + mlh1 vs rest
ggbiplot::ggbiplot(pc.total, var.scale=1, groups=(total$combo == "MLH1:Basal"),
                   ellipse=F, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="MLH1:Basal vs Rest") +
  geom_point(aes(shape=total$type, color=(total$combo == "MLH1:Basal"), size=total$type)) +
  scale_shape_manual(values=c("Mouse"=3,"Human"=1)) + scale_size_manual(values=c("Mouse"=5,"Human"=1))

#luminal - msh2 vs rest
ggbiplot::ggbiplot(pc.luminal, var.scale=1, groups=luminal$MMR,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0.5) +
  theme_bw() + labs(title="Luminal Patients - Only PAM50") +
  scale_color_manual(values=c("None"="grey", "MLH1"="red", "MSH2"="blue"), labels=c("None"="Rest"))

#luminal + msh2 vs rest
ggbiplot::ggbiplot(pc.total, var.scale=1, groups=(total$combo == "MSH2:Lum"),
                   ellipse=F, circle=F, var.axes=F, varname.size=0) +
  theme_bw() + labs(title="MSH2:Lum vs Rest") +
  geom_point(aes(shape=total$type, color=(total$combo == "MSH2:Lum"), size=total$type)) +
  scale_shape_manual(values=c("Mouse"=3,"Human"=1)) + scale_size_manual(values=c("Mouse"=5,"Human"=1))


ggbiplot::ggbiplot(pc.total, var.scale=1, groups=total$pam50,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0.4) +
  theme_bw() + labs(title="All Patients") +
  geom_point(aes(shape=(total$combo == "MLH1:Basal" | total$combo == "MSH2:Lum"), color=total$pam50,
                 size=(total$combo == "MLH1:Basal" | total$combo == "MSH2:Lum"))) +
  scale_shape_manual(values=c("TRUE"=3,"FALSE"=1)) + scale_size_manual(values=c("TRUE"=5,"FALSE"=1)) +
  theme(legend.position="none")

ggbiplot::ggbiplot(pc.total, var.scale=1, groups=k.3,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0.4) +
  theme_bw() + labs(title="All Patients") +
  geom_point(aes(shape=(total$combo == "MLH1:Basal" | total$combo == "MSH2:Lum"), color=k.3,
                 size=(total$combo == "MLH1:Basal" | total$combo == "MSH2:Lum"))) +
  scale_shape_manual(values=c("TRUE"=3,"FALSE"=1)) + scale_size_manual(values=c("TRUE"=5,"FALSE"=1)) +
  theme(legend.position="none")



#kmeans -----
wss <- numeric(12)
for(i in 1:12){
  k.out <- kmeans(total[,-c(1:4)], centers=i, nstart=20)
  wss[i] <- k.out$tot.withinss
}
rm(k.out)
wss.df <- tibble(clusters=1:12, wss=wss)

ggplot(data=wss.df, aes(x=clusters, y=wss)) + geom_point() + geom_line()

k.3 <- kmeans(total[,-c(1:6)], centers=3, nstart=20)$cluster
k.3 <- as.factor(k.3)
k.4 <- kmeans(total[,-c(1:6)], centers=4, nstart=20)$cluster
k.4 <- as.factor(k.4)
