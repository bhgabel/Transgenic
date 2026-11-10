library(tidyverse)
library(ggbiplot)
library("DESeq2")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

tcga <- read.csv(file="tcga_pam50.csv", header=T)

rna <- read.table(file="brca_tcga_gdc/data_mrna_seq_read_counts.txt", header=T,
                  sep='\t')

genes <- clusterProfiler::bitr(rna$Entrez_Gene_Id, fromType="ENTREZID",
                               toType="SYMBOL", OrgDb=org.Hs.eg.db)

rna$Entrez_Gene_Id <- as.character(rna$Entrez_Gene_Id)
rna <- left_join(genes, rna, by=c("ENTREZID" = "Entrez_Gene_Id"))
rownames(rna) <- rna$SYMBOL
rna <- data.frame(t(rna[-(1:2)]))

rna <- rownames_to_column(rna, var="Patient ID")

rna$Patient.ID <- str_sub(rna$`Patient ID`, end=-5) %>% 
  str_replace_all(pattern="\\.", replacement="-")

total <- left_join(tcga, rna)
total <- distinct(total, Patient.ID, .keep_all=T)
total <- drop_na(total, "Patient ID")

c.cycle <- dplyr::select(total, 1:3, any_of(ccycle$Gene))
c.cycle <- as.matrix(t(c.cycle[-(1:3)]))
vst <- vst(c.cycle, nsub=100)

cellcycle.vst <- as.data.frame(t(c.cycle))
cellcycle.vst <- data.frame(total$Patient.ID, total$PAM50, cellcycle.vst)

pc.ccycle <- prcomp(cellcycle.vst[,-(1:3)], center=T, scale.=T)
ggbiplot(pc.ccycle, groups=NULL,
         ellipse=T, ellipse.fill=F,
         var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=cellcycle.vst$PAM50, shape=cellcycle.vst$MMR), alpha=0.8) +
  scale_color_manual(values=c('#CC79A7', '#c29f1f', '#27AD81', '#31688E', '#471164')) +
  scale_shape_manual(values=c("MLH1"=3, "MSH2"=24, "None"=1),
                     labels=c("MLH1 Low", "MSH2 Low", "Rest")) +
  theme_bw() + labs(title="TCGA by Cell Cycle Genes", color="PAM50", shape="MMR")

ggsave("Images/TCGA_PCA_cell_cycle.tiff", dpi=600)

temp <- subset(cellcycle.vst, MMR != "None")
pc.ccycle <- prcomp(temp[,-(1:3)], center=T, scale.=T)
ggbiplot(pc.ccycle, groups=NULL,
         ellipse=T, ellipse.fill=F,
         var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=temp$PAM50, shape=temp$MMR), alpha=1, size=2) +
  scale_color_manual(values=c('#CC79A7', '#c29f1f', '#27AD81', '#31688E', '#471164')) +
  scale_shape_manual(values=c("MLH1"=3, "MSH2"=24),
                     labels=c("MLH1 Low", "MSH2 Low")) +
  theme_bw() + labs(title="TCGA by Cell Cycle Genes", color="PAM50", shape="MMR")

ggsave("Images/TCGA_PCA_cell_cycle.tiff", dpi=600)


pc.pam <- prcomp(gdc[, (22:71)], center=T, scale.=T)
ggbiplot(pc.pam, groups=gdc$MMR, ellipse=F, var.axes=F, varname.size=0)


basal <- subset(gdc, PAM50 == "Basal")
pc.basal <- prcomp(basal[, (22:71)], center=T, scale.=T)
ggbiplot::ggbiplot(pc.basal, var.scale=1,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=basal$MMR)) +
  theme_bw() +
  labs(title="Basal Patients by PAM50 Genes", color="MMR") +
  scale_color_manual(values=c("None"="grey", "MLH1"='#3953A4', "MSH2"='#FECB67'),
                     labels=c("None"="Rest", "MLH1"="MLH1 Low", "MSH2"="MSH2 Low"))


lum <- subset(gdc, PAM50 == "LumA")
pc.lum <- prcomp(lum[, (22:71)], center=T, scale.=T)
ggbiplot::ggbiplot(pc.lum, var.scale=1,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=lum$MMR, fill=lum$MMR), shape=21, alpha=0.8) +
  theme_bw() +
  labs(title="Luminal A Patients by PAM50 Genes", fill="MMR") +
  scale_color_manual(values=c("None"="grey", "MLH1"='#3953A4', "MSH2"='#FECB67'),
                     guide="none") +
  scale_fill_manual(values=c("None"="white",  "MLH1"='#3953A4', "MSH2"='#FECB67'),
                    labels=c("None"="Rest", "MLH1"="MLH1 Low", "MSH2"="MSH2 Low"))

ggsave(file="Images/TCGA_PCA_lumA.tiff", dpi=600)


basal <- subset(cellcycle.vst, PAM50 == "Basal" & !is.na(MMR))
pc.basal <- prcomp(basal[, -(1:3)], center=T, scale.=T)
ggbiplot::ggbiplot(pc.basal, var.scale=1,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=basal$MMR)) +
  theme_bw() +
  labs(title="Basal Patients by Cell Cycle Genes", color="MMR") +
  scale_color_manual(values=c("None"="grey", "MLH1"='#3953A4', "MSH2"='#FECB67'),
                     labels=c("None"="Rest", "MLH1"="MLH1 Low", "MSH2"="MSH2 Low"))


lum <- subset(cellcycle.vst, PAM50 == "LumA" & !is.na(MMR))
pc.lum <- prcomp(lum[, -(1:3)], center=T, scale.=T)
ggbiplot::ggbiplot(pc.lum, var.scale=1,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=lum$MMR, fill=lum$MMR), shape=21, alpha=0.8) +
  theme_bw() +
  labs(title="Luminal A Patients by Cell Cycle Genes", fill="MMR") +
  scale_color_manual(values=c("None"="grey", "MLH1"='#3953A4', "MSH2"='#FECB67'),
                     guide="none") +
  scale_fill_manual(values=c("None"="white",  "MLH1"='#3953A4', "MSH2"='#FECB67'),
                    labels=c("None"="Rest", "MLH1"="MLH1 Low", "MSH2"="MSH2 Low"))




cellcycle <- read.table(file="tcga cell cycle mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).tsv",
                        sep='\t', header=T)

cellcycle <- left_join(cellcycle, select(gdc, "Sample.ID", "MMR", "PAM50"), by=c("SAMPLE_ID"="Sample.ID"))

cellcycle <- drop_na(cellcycle, MMR)

cellcycle$MMR <- factor(cellcycle$MMR, levels=c("None", "MLH1", "MSH2"))

sink("TCGA_heatmap_stats.txt")
print("Average Gene Expression by MMR")
ccmeans <- cellcycle %>% group_by(MMR) %>% summarise_at(vars(CCNC:CCNB1), mean, na.rm=T)
ccmeans <- as.data.frame(ccmeans) %>% column_to_rownames(var="MMR")

basal <- cellcycle %>% filter(PAM50 == "Basal") %>% group_by(MMR) %>%
  summarise_at(vars(CCNC:CCNB1), mean, na.rm=T)
basal <- as.data.frame(basal) %>% column_to_rownames(var="MMR")

luma <- cellcycle %>% filter(PAM50 == "LumA") %>% group_by(MMR) %>%
  summarise_at(vars(CCNC:CCNB1), mean, na.rm=T)
luma <- as.data.frame(luma) %>% column_to_rownames(var="MMR")

lumb <- cellcycle %>% filter(PAM50 == "LumB") %>% group_by(MMR) %>%
  summarise_at(vars(CCNC:CCNB1), mean, na.rm=T)
lumb <- as.data.frame(lumb) %>% column_to_rownames(var="MMR")

her2 <- cellcycle %>% filter(PAM50 == "Her2") %>% group_by(MMR) %>%
  summarise_at(vars(CCNC:CCNB1), mean, na.rm=T)
her2 <- as.data.frame(her2) %>% column_to_rownames(var="MMR")

norm <- cellcycle %>% filter(PAM50 == "Normal") %>% group_by(MMR) %>%
  summarise_at(vars(CCNC:CCNB1), mean, na.rm=T)
norm <- as.data.frame(norm) %>% column_to_rownames(var="MMR")


print("Wilcox test for each Gene")
lapply(cellcycle[1:10], function(col){
       pairwise.wilcox.test(col, cellcycle$MMR)})
sink()

cellcycle <- column_to_rownames(cellcycle, var="SAMPLE_ID") %>% arrange("MMR")


annot_col <- select(cellcycle, "MMR", "PAM50")
annot_row <- data.frame(Phase=factor(c("G0/G1", "G0/G1", "G1/S", "G1/S", "G1/S",
                                "G1/S", "G1/S/G2", "S/G2", "M", "M"),
                                levels=c("G0/G1", "G1/S", "G1/S/G2", "S/G2", "M"),
                                ordered=T),
                        row.names=colnames(cellcycle)[1:10])

color <- list("MMR"=c("None"="grey", "MLH1"="blue", "MSH2"="yellow"),
              "PAM50"=c("Basal"='#E3E418', "Normal"='#27AD81',
                        "Her2"='#31688E', "LumA"='#443A83', "LumB"='#471164'))

pheatmap(t(cellcycle[c(1:10)]), scale="row",
         color=colorRampPalette(c("blue", "white", "red"))(100),
         show_colnames=F, cluster_rows=F, #cluster_cols=F,
         annotation_row=annot_row, legend=F,
         annotation_col=annot_col, annotation_colors=color,
         main="Differential Cell Cycle Markers in TCGA",
         filename="Images/TCGA_heatmap_cellcycle_clustered.png")

pheatmap(mat=t(ccmeans),
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows=F, cluster_cols=F,
         annotation_row=annot_row, annotation_names_row=F,
         display_numbers=T, number_color="black",
         labels_col=c("MLH1", "MSH2", "Rest"), angle_col=0,
         main="Mean Expression by MMR Deficiency",
        filename="Images/TCGA_heatmap_cellcycle_mean.png")

pheatmap(mat=t(basal),
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows=F, cluster_cols=F,
         annotation_row=annot_row, annotation_names_row=F,
         display_numbers=T, number_color="black",
         labels_col=c("MLH1", "MSH2", "Rest"), angle_col=0,
         main="Basal Mean Expression by MMR Deficiency",
         filename="Images/TCGA_heatmap_cellcycle_basal.png")


pc.cc <- prcomp(cellcycle[, c(3:12)], center=T, scale.=T)

ggbiplot::ggbiplot(pc.cc, var.scale=1,
                   ellipse=F, circle=F, var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=cellcycle$MMR, fill=cellcycle$MMR), shape=21, alpha=0.8) +
  theme_bw() +
  labs(title="TCGA by Cell Cycle Genes", fill="MMR") +
  scale_color_manual(values=c("None"="grey", "MLH1"='#3953A4', "MSH2"='#FECB67'),
                     guide="none") +
  scale_fill_manual(values=c("None"="white",  "MLH1"='#3953A4', "MSH2"='#FECB67'),
                    labels=c("None"="Rest", "MLH1"="MLH1 Low", "MSH2"="MSH2 Low"))

ggsave("Images/TCGA_PCA_cellcycle_mmr.png", dpi=600)

ggbiplot::ggbiplot(pc.cc, var.scale=1,
                   groups=cellcycle$PAM50,
                   ellipse=F, ellipse.linewidth=0.8, ellipse.fill=F,
                   circle=F, var.axes=F, varname.size=0, alpha=0) +
  geom_point(aes(color=cellcycle$PAM50, shape=cellcycle$MMR, alpha=cellcycle$MMR)) +
  theme_bw() +
  labs(title="TCGA by Cell Cycle Genes", color="PAM50", shape="MMR") +
  scale_color_manual(values=c('#CC79A7', '#c29f1f', '#27AD81', '#31688E', '#471164')) +
  scale_shape_manual(values=c("None"=21, "MLH1"=3, "MSH2"=24),
                     labels=c("Rest", "MLH1 Low", "MSH2 Low")) +
  scale_alpha_manual(values=c("None"=0.4, "MLH1"=0.9, "MSH2"=0.9), guide="none")

ggsave("Images/TCGA_PCA_cellcycle_pam50.png", dpi=600)
