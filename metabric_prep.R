library(readr)
library(tidyverse)

mrna <- read_delim("Metabric mRNA expression z-scores relative to all samples (log microarray).txt", 
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE)

clinical <- read_delim("metabric_clinical_data.tsv", delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

meta <- left_join(mrna, clinical, by=c("SAMPLE_ID" = "Sample ID"))

remove(mrna, clinical)

meta <- select(meta, c("SAMPLE_ID", "MLH1", "MSH2", "Mutation Count", "Pam50 + Claudin-low subtype",
                       "ER Status", "HER2 Status", "PR Status", "TMB (nonsynonymous)"))

colnames(meta) <- c("Sample ID", "MLH1", "MSH2", "Mutations", "PAM50", "ER", "HER2", "PR", "TMB")

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
  if(meta$ER[[i]] == 'Negative' & meta$PR[[i]] == 'Negative' & meta$HER2[[i]] == 'Negative'){
    meta$HR[[i]] <- 'TNBC'
  } else if((meta$ER[[i]] == 'Positive' | meta$PR[[i]] == 'Positive') & meta$HER2[[i]] == 'Negative'){
    meta$HR[[i]] <- 'HR+/HER2-'
  } else if((meta$ER[[i]] == 'Positive' | meta$PR[[i]] == 'Positive') & meta$HER2[[i]] == 'Positive'){
    meta$HR[[i]] <- 'HR+/HER2+'
  } else{
    meta$HR[[i]] <- 'HR-/HER2+'
  }
}

#pam50
meta$PAM50 <- as.factor(meta$PAM50)
meta$pam50_f <- factor("Other", levels=c("Basal", "Luminal", "Other"))
for(i in 1:length(meta$PAM50)){
  if(is.na(meta$PAM50[[i]])){
    meta$pam50_f[[i]] <- "Other"
  } else if(meta$PAM50[[i]] == "Basal"){
    meta$pam50_f[[i]] <- "Basal"
  } else if(meta$PAM50[[i]] == "LumA" | meta$PAM50[[i]] == "LumB"){
    meta$pam50_f[[i]] <- "Luminal"
  } else{
    meta$pam50_f[[i]] <- "Other"
  }
}

meta <- drop_na(meta)

write.csv(meta, "metabric_prepped.csv")
