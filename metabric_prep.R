library(readr)
library(tidyverse)
library(genefu)

# mrna <- read_delim("Metabric mRNA expression z-scores relative to all samples (log microarray).txt",
#                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

mrna <- read_delim("Metabric mRNA expression (Illumina HT-12 v3 microarray).txt",
                   delim="\t", escape_double = FALSE, trim_ws = TRUE)

clinical <- read_delim("metabric_clinical_data.tsv", delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

meta <- left_join(mrna, clinical, by=c("SAMPLE_ID" = "Sample ID"))

remove(mrna, clinical)

meta <- dplyr::select(meta, c(2:54, "Mutation Count", "Pam50 + Claudin-low subtype",
                       "ER Status", "HER2 Status", "PR Status", "TMB (nonsynonymous)",
                       "Relapse Free Status", "Relapse Free Status (Months)"))

colnames(meta)[c(1, 54:61)] <- c("Sample ID", "Mutations", "Pam50", "ER", "HER2",
                                 "PR", "TMB", "Relapse.Status", "Relapse.Months")

#remove NA
meta <- drop_na(meta, MLH1)
# Establish factors
#MMR low
MLH1_prob <- quantile(meta$MLH1, prob=0.12, na.rm=T)
MSH2_prob <- quantile(meta$MSH2, prob=0.08, na.rm=T)
meta$MLH1_low <- (meta$MLH1 <= MLH1_prob)
meta$MSH2_low <- (meta$MSH2 <= MSH2_prob)
meta$MMR <- factor("None", levels=c("MLH1", "MSH2", "None"))
for(i in 1:length(meta$MMR)){
  if(meta$MLH1_low[[i]] && meta$MSH2_low[[i]]){
    meta$MMR[[i]] <- NA
  } else if(meta$MLH1_low[[i]]){
    meta$MMR[[i]] <- "MLH1"
  } else if(meta$MSH2_low[[i]]){
    meta$MMR[[i]] <- "MSH2"
  }
}
meta <- drop_na(meta, "MMR")

meta$MLH1_low.subtype <- FALSE
basal <- quantile(subset(meta, PAM50=="Basal")$MLH1, probs=0.12)
normal <- quantile(subset(meta, PAM50=="Normal")$MLH1, probs=0.12)
her2 <- quantile(subset(meta, PAM50=="Her2")$MLH1, probs=0.12)
luma <- quantile(subset(meta, PAM50=="LumA")$MLH1, probs=0.12)
lumb <- quantile(subset(meta, PAM50=="LumB")$MLH1, probs=0.12)

for(i in 1:length(meta$MLH1)){
  if(is.na(meta$PAM50[[i]])){
    #pass
  }
  else if(meta$PAM50[[i]] == "Basal"){
    meta$MLH1_low.subtype[[i]] <- (meta$MLH1[[i]] <= basal)
  }
  else if(meta$PAM50[[i]] == "Normal"){
    meta$MLH1_low.subtype[[i]] <- (meta$MLH1[[i]] <= normal)
  }
  else if(meta$PAM50[[i]] == "Her2"){
    meta$MLH1_low.subtype[[i]] <- (meta$MLH1[[i]] <= her2)
  }
  else if(meta$PAM50[[i]] == "LumA"){
    meta$MLH1_low.subtype[[i]] <- (meta$MLH1[[i]] <= luma)
  }
  else if(meta$PAM50[[i]] == "LumB"){
    meta$MLH1_low.subtype[[i]] <- (meta$MLH1[[i]] <= lumb)
  }
}

meta$MSH2_low.subtype <- FALSE
basal <- quantile(subset(meta, PAM50=="Basal")$MSH2, probs=0.08)
normal <- quantile(subset(meta, PAM50=="Normal")$MSH2, probs=0.08)
her2 <- quantile(subset(meta, PAM50=="Her2")$MSH2, probs=0.08)
luma <- quantile(subset(meta, PAM50=="LumA")$MSH2, probs=0.08)
lumb <- quantile(subset(meta, PAM50=="LumB")$MSH2, probs=0.08)

for(i in 1:length(meta$MSH2)){
  if(is.na(meta$PAM50[[i]])){
    #pass
  }
  else if(meta$PAM50[[i]] == "Basal"){
    meta$MSH2_low.subtype[[i]] <- (meta$MSH2[[i]] <= basal)
  }
  else if(meta$PAM50[[i]] == "Normal"){
    meta$MSH2_low.subtype[[i]] <- (meta$MSH2[[i]] <= normal)
  }
  else if(meta$PAM50[[i]] == "Her2"){
    meta$MSH2_low.subtype[[i]] <- (meta$MSH2[[i]] <= her2)
  }
  else if(meta$PAM50[[i]] == "LumA"){
    meta$MSH2_low.subtype[[i]] <- (meta$MSH2[[i]] <= luma)
  }
  else if(meta$PAM50[[i]] == "LumB"){
    meta$MSH2_low.subtype[[i]] <- (meta$MSH2[[i]] <= lumb)
  }
}

meta$MMR.subtype <- factor("None", levels=c("MLH1", "MSH2", "None"))
for(i in 1:length(meta$MMR.subtype)){
  if(meta$MLH1_low[[i]] && meta$MSH2_low[[i]]){
    meta$MMR.subtype[[i]] <- NA
  } else if(meta$MLH1_low[[i]]){
    meta$MMR.subtype[[i]] <- "MLH1"
  } else if(meta$MSH2_low[[i]]){
    meta$MMR.subtype[[i]] <- "MSH2"
  }
}

#Hormone status
meta$HR <- factor(NA, levels=c('HR+/HER2-', 'HR+/HER2+', 'HR-/HER2+', 'TNBC'))
for(i in 1:length(meta$HR)){
  if(is.na(meta$ER[[i]]) | is.na(meta$PR[[i]]) | is.na(meta$HER2[[i]])){
    meta$HR[[i]] <- NA
  } else if(meta$ER[[i]] == 'Negative' & meta$PR[[i]] == 'Negative' &
            meta$HER2[[i]] == 'Negative'){
    meta$HR[[i]] <- 'TNBC'
  } else if((meta$ER[[i]] == 'Positive' | meta$PR[[i]] == 'Positive') &
            meta$HER2[[i]] == 'Negative'){
    meta$HR[[i]] <- 'HR+/HER2-'
  } else if((meta$ER[[i]] == 'Positive' | meta$PR[[i]] == 'Positive') &
            meta$HER2[[i]] == 'Positive'){
    meta$HR[[i]] <- 'HR+/HER2+'
  } else if(meta$ER[[i]] == 'Negative' & meta$PR[[i]] == 'Negative' &
            meta$HER2[[i]] == 'Positive'){
    meta$HR[[i]] <- 'HR-/HER2+'
  } else {
    meta$HR[[i]] <- NA
  }
}

#pam50
meta$Pam50 <- as.factor(meta$Pam50)

#reclassify pam50 without claudin-low type
data("pam50.robust")
gene_info <- read.csv("pam50_gene_list.csv")
meta <- rename(meta, all_of(c("CDCA1" = "NUF2", "KNTC2" = "NDC80", "ORC6L" = "ORC6")))
pam50.pred <- molecular.subtyping(sbt.model = "pam50", data = meta[,c(4:53)],
                                  annot = gene_info, do.mapping = TRUE)

meta$PAM50 <- as.factor(pam50.pred$subtype)

meta$PAM50 <- factor(meta$Pam50, levels=c("Basal", "Normal", "Her2", "LumA", "LumB"))

for(i in 1:length(meta$PAM50)){
  if(meta$Pam50[[i]] == "claudin-low"){
    meta$PAM50[[i]] <- pam50.pred$subtype[[i]]
  }
  else if(meta$Pam50[[i]] == "NC"){
    meta$PAM50[[i]] <- pam50.pred$subtype[[i]]
  }
  else{
    meta$PAM50[[i]] <- meta$Pam50[[i]]
  }
}

#remove intermediate data frames
rm(pam50.pred, pam50.robust, mrna, clinical, gene_info)

write.csv(meta, "metabric_prepped.csv", row.names=F)
