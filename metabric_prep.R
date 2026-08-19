library(readr)
library(tidyverse)
library(genefu)

mrna <- read_delim("Metabric mRNA expression z-scores relative to all samples (log microarray).txt",
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE)

clinical <- read_delim("metabric_clinical_data.tsv", delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

meta <- left_join(mrna, clinical, by=c("SAMPLE_ID" = "Sample ID"))

remove(mrna, clinical)

meta <- dplyr::select(meta, c(2:54, "Mutation Count", "Pam50 + Claudin-low subtype",
                       "ER Status", "HER2 Status", "PR Status", "TMB (nonsynonymous)",
                       "Relapse Free Status", "Relapse Free Status (Months)"))

colnames(meta)[c(1, 54:61)] <- c("Sample ID", "Mutations", "Pam50", "ER", "HER2",
                                 "PR", "TMB", "Relapse.Status", "Relapse.Months")

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

#Hormone status
meta$ER <- as.factor(meta$ER)
meta$HER2 <- as.factor(meta$HER2)
meta$PR <- as.factor(meta$PR)
meta$HR <- factor('HR+/HER2+', levels=c('HR+/HER2+', 'HR+/HER2-', 'HR-/HER2+', 'TNBC'))
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
meta$pam50_f <- factor("Other", levels=c("Basal", "Luminal", "Other"))
for(i in 1:length(meta$Pam50)){
  if(is.na(meta$Pam50[[i]])){
    meta$pam50_f[[i]] <- "Other"
  } else if(meta$Pam50[[i]] == "Basal"){
    meta$pam50_f[[i]] <- "Basal"
  } else if(meta$Pam50[[i]] == "LumA" | meta$Pam50[[i]] == "LumB"){
    meta$pam50_f[[i]] <- "Luminal"
  } else{
    meta$pam50_f[[i]] <- "Other"
  }
}

#reclassify pam50 without claudin-low type
data("pam50.robust")
gene_info <- read.csv("pam50_gene_list.csv")
meta <- rename(meta, all_of(c("CDCA1" = "NUF2", "KNTC2" = "NDC80", "ORC6L" = "ORC6")))
pam50.pred <- molecular.subtyping(sbt.model = "pam50", data = meta[,c(4:53)],
                                  annot = gene_info, do.mapping = TRUE)

meta$PAM50 <- as.factor(pam50.pred$subtype)

meta <- drop_na(meta)

write.csv(meta, "metabric_prepped.csv", row.names=F)
