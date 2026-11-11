library(tidyverse)
library(genefu)

#data downloaded from https://data.mendeley.com/datasets/w8tsvd5zrf/1

clinical <- read.table("SCAN-B/Clinical Metadata/SCANB_clinical_metadata.txt",
                       sep='\t', header=T)

#get gene ids for subset of genes of interest: MLH1, MSH2, BRCA1, and PAM50
genes <- read.table("SCAN-B/Gene Annotations/Gene.ID.ann.txt",
                    sep='\t', header=T)

selection <- read.csv("pam50_gene_list.csv")

gene.ids <- genes$Gene.ID[genes$EntrezGene %in% selection$EntrezGene.ID |
                          genes$Gene.Name %in% c("MLH1", "MSH2", "BRCA1")]



#gene expression table is REALLY big
rna <- read.table("SCAN-B/Gene-expression matrix (protocol-adjusted)/bc4genematrix_noNeg.txt",
                  sep='\t', header=T)

#subset for only genes of interest
rna <- filter(rna, Gene.ID %in% gene.ids)

#get gene names from gene id
rna <- left_join(rna, genes[, 1:2])

#transpose and remove gene id
rna <- column_to_rownames(rna, var='Gene.Name')

rna <- as.data.frame(t(rna))

rna <- rna[-1, ]

rna <- rownames_to_column(rna, var='StringTieName')

#combine rna and clinical data
scanb <- left_join(rna, clinical)

#cast all values to numeric
scanb[, 2:54] <- lapply(scanb[, 2:54], as.numeric)

rm(clinical, genes, rna, selection, gene.ids)


#classify patients with PAM50
data("pam50.robust")
gene_info <- read.csv("pam50_gene_list.csv")

#rename genes to match subtyping call
temp <- scanb %>%
  rename(all_of(c("CDCA1" = "NUF2", "KNTC2" = "NDC80", "ORC6L" = "ORC6"))) %>% 
  dplyr::select(any_of(gene_info$probe))

pam50.pred <- molecular.subtyping(sbt.model = "pam50", data = temp,
                                  annot = gene_info, do.mapping = TRUE)

scanb$PAM50 <- factor(pam50.pred$subtype,
                      levels=c("Basal", "Normal", "Her2", "LumA", "LumB"))

write.csv(pam50.pred$subtype.proba, file="scanb_pam50_prob.csv", row.names=F)
rm(gene_info, pam50.robust, pam50.pred, temp)

#set factors for low RNA expression


scanb$MLH1_low <- (scanb$MLH1 <= quantile(scanb$MLH1, probs=0.12, na.rm=T))
scanb$MSH2_low <- (scanb$MSH2 <= quantile(scanb$MSH2, probs=0.08, na.rm=T))
scanb$BRCA1_low <- (scanb$BRCA1 <= quantile(scanb$BRCA1, probs=0.10, na.rm=T))

#set MMR factor, NA for both low
scanb$MMR <- factor(NA, levels=c("MLH1", "MSH2", "None"))
for(i in 1:length(scanb$MMR)){
  if(is.na(scanb$MLH1_low[[i]])){
    scanb$MMR[[i]] <- NA
  }
  else if(scanb$MLH1_low[[i]] & scanb$MSH2_low[[i]]){
    scanb$MMR[[i]] <- NA
  }
  else if(scanb$MLH1_low[[i]]){
    scanb$MMR[[i]] <- "MLH1"
  }
  else if(scanb$MSH2_low[[i]]){
    scanb$MMR[[i]] <- "MSH2"
  }
  else{
    scanb$MMR[[i]] <- "None"
  }
}


scanb$MLH1_low.subtype <- FALSE
basal <- quantile(subset(scanb, PAM50=="Basal")$MLH1, probs=0.12)
normal <- quantile(subset(scanb, PAM50=="Normal")$MLH1, probs=0.12)
her2 <- quantile(subset(scanb, PAM50=="Her2")$MLH1, probs=0.12)
luma <- quantile(subset(scanb, PAM50=="LumA")$MLH1, probs=0.12)
lumb <- quantile(subset(scanb, PAM50=="LumB")$MLH1, probs=0.12)

for(i in 1:length(scanb$MLH1)){
  if(is.na(scanb$PAM50[[i]])){
    #pass
  }
  else if(scanb$PAM50[[i]] == "Basal"){
    scanb$MLH1_low.subtype[[i]] <- (scanb$MLH1[[i]] <= basal)
  }
  else if(scanb$PAM50[[i]] == "Normal"){
    scanb$MLH1_low.subtype[[i]] <- (scanb$MLH1[[i]] <= normal)
  }
  else if(scanb$PAM50[[i]] == "Her2"){
    scanb$MLH1_low.subtype[[i]] <- (scanb$MLH1[[i]] <= her2)
  }
  else if(scanb$PAM50[[i]] == "LumA"){
    scanb$MLH1_low.subtype[[i]] <- (scanb$MLH1[[i]] <= luma)
  }
  else if(scanb$PAM50[[i]] == "LumB"){
    scanb$MLH1_low.subtype[[i]] <- (scanb$MLH1[[i]] <= lumb)
  }
}

scanb$MSH2_low.subtype <- FALSE
basal <- quantile(subset(scanb, PAM50=="Basal")$MSH2, probs=0.08)
normal <- quantile(subset(scanb, PAM50=="Normal")$MSH2, probs=0.08)
her2 <- quantile(subset(scanb, PAM50=="Her2")$MSH2, probs=0.08)
luma <- quantile(subset(scanb, PAM50=="LumA")$MSH2, probs=0.08)
lumb <- quantile(subset(scanb, PAM50=="LumB")$MSH2, probs=0.08)

for(i in 1:length(scanb$MSH2)){
  if(is.na(scanb$PAM50[[i]])){
    #pass
  }
  else if(scanb$PAM50[[i]] == "Basal"){
    scanb$MSH2_low.subtype[[i]] <- (scanb$MSH2[[i]] <= basal)
  }
  else if(scanb$PAM50[[i]] == "Normal"){
    scanb$MSH2_low.subtype[[i]] <- (scanb$MSH2[[i]] <= normal)
  }
  else if(scanb$PAM50[[i]] == "Her2"){
    scanb$MSH2_low.subtype[[i]] <- (scanb$MSH2[[i]] <= her2)
  }
  else if(scanb$PAM50[[i]] == "LumA"){
    scanb$MSH2_low.subtype[[i]] <- (scanb$MSH2[[i]] <= luma)
  }
  else if(scanb$PAM50[[i]] == "LumB"){
    scanb$MSH2_low.subtype[[i]] <- (scanb$MSH2[[i]] <= lumb)
  }
}

scanb$MMR.subtype <- factor("None", levels=c("MLH1", "MSH2", "None"))
for(i in 1:length(scanb$MMR.subtype)){
  if(scanb$MLH1_low[[i]] & scanb$MSH2_low[[i]]){
    scanb$MMR.subtype[[i]] <- NA
  } else if(scanb$MLH1_low[[i]]){
    scanb$MMR.subtype[[i]] <- "MLH1"
  } else if(scanb$MSH2_low[[i]]){
    scanb$MMR.subtype[[i]] <- "MSH2"
  }
}


write.csv(scanb, file="scanb_prepped.csv", row.names=F)
