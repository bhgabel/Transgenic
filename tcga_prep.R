library(readr)
library(readxl)
library(tidyverse)
library(TCGAbiolinks)

#Creating and subsetting data frames ----
#FPKM z-scores gene expression
mrna <- read_delim("pam50 mRNA expression fpkm Zscores.tsv", 
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(c('SAMPLE_ID', 'MLH1', 'MSH2'))

#Clinical info
clinical <- read_delim("brca_tcga_gdc_clinical_data.tsv", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(c("Patient ID", "Sample ID", "Mutation Count", "Fraction Genome Altered"))

#HR status                          
firehose <- read_excel("~/Documents/Haricharan/PanCancer/TCGA_Firehose.xlsx", sheet = "Clinical") %>% 
  dplyr::select(c("Sample ID", "ER Status By IHC", "IHC-HER2", "PR status by ihc"))

#Combine into one data frame, rename for pam50 subtyping
gdc <- right_join(clinical, mrna, by=c("Sample ID" = "SAMPLE_ID"))
#rename(all_of(c("CDCA1" = "NUF2", "KNTC2" = "NDC80", "ORC6L" = "ORC6")))
# gdc <- left_join(gdc, mrna, by=c("Sample ID" = "SAMPLE_ID"))

#Add MSI MANTIS scores
mantis <- read_excel("MSI_Mantis.xlsx") %>%
  dplyr::select(c("Case ID", "MANTIS Score"))

gdc <- left_join(gdc, mantis, by=c("Patient ID" = "Case ID"))

gdc <- left_join(gdc, firehose)

#Add CIN data
CIN_data <- read_excel("CIN_data.xlsx", sheet = "ST_20_TCGA_Activities_scaled")
colnames(CIN_data)[1] <- 'Patient ID'

gdc <- left_join(gdc, CIN_data[c('Patient ID', 'Total', 'CX1', 'CX6')])


# gene_info <- data.frame(colnames(firehose)[-c(1:5)], pam50$centroids.map$EntrezGene.ID)
# colnames(gene_info) <- c("probe", "EntrezGene.ID")
# 
# pam50_pred <- molecular.subtyping(sbt.model="pam50", data=firehose, annot=gene_info, do.mapping=T, verbose=T)

#Pull pam50 subtypes from TCGA package
subtypes <- TCGAquery_subtype("brca") %>% dplyr::select(c("patient", "BRCA_Subtype_PAM50"))

#Combine to main data frame
gdc <- left_join(gdc, subtypes, by=c("Patient ID" = "patient"))

#Rename columns
colnames(gdc) <- c('Patient ID', 'Sample ID', 'Mutations', 'FGA', 'MLH1', 'MSH2', 'MANTIS', 'ER', 'HER2', 'PR', 'CIN', 'CX1', 'CX6', 'PAM50')

#remove primary data frames, since everything joined into gdc
remove(list = c('CIN_data', 'clinical', 'firehose', 'mantis', 'mrna', 'mrna_pam50', 'subtypes'))

#Assign factors ----
#Assign factors + basal/lum
gdc$PAM50 <- as.factor(gdc$PAM50)
gdc$pam50_f <- factor("Other", levels=c("Basal", "Luminal", "Other"))
for(i in 1:length(gdc$PAM50)){
  if(is.na(gdc$PAM50[[i]])){
    gdc$pam50_f[[i]] <- "Other"
  } else if(gdc$PAM50[[i]] == "Basal"){
    gdc$pam50_f[[i]] <- "Basal"
  } else if(gdc$PAM50[[i]] == "LumA" | gdc$PAM50[[i]] == "LumB"){
    gdc$pam50_f[[i]] <- "Luminal"
  } else{
    gdc$pam50_f[[i]] <- "Other"
  }
}

gdc <- drop_na(gdc, "MLH1", "MSH2", "Mutations")

#Assign factors for low gene expression and remove combined loss
gdc$MLH1_low <- (gdc$MLH1 <= quantile(gdc$MLH1, probs=0.12, na.rm=T))
gdc$MSH2_low <- (gdc$MSH2 <= quantile(gdc$MSH2, probs=0.08, na.rm=T))
gdc$p53_low <- (gdc$p53 <= quantile(gdc$p53, probs=0.5, na.rm=T))
gdc$ERBB2_low <- (gdc$ERBB2 <= quantile(gdc$ERBB2, probs=0.5, na.rm=T))
gdc$ESR1_low <- (gdc$ESR1 <= quantile(gdc$ESR1, probs=0.5, na.rm=T))
gdc$MMR <- factor("None", levels=c("MLH1", "MSH2", "None"))
for(i in 1:length(gdc$MMR)){
  if(gdc$MLH1_low[[i]] && gdc$MSH2_low[[i]]){
    gdc$MMR[[i]] <- NA
  } else if(gdc$MLH1_low[[i]]){
    gdc$MMR[[i]] <- "MLH1"
  } else if(gdc$MSH2_low[[i]]){
    gdc$MMR[[i]] <- "MSH2"
  }
}
gdc <- drop_na(gdc, "MMR")


#Assign factors for quintiles [0-20, 20-40,...]
MLH1_probs <- quantile(gdc$MLH1, probs=c(0.2, 0.4, 0.6, 0.8), na.rm=T)
MSH2_probs <- quantile(gdc$MSH2, probs=c(0.2, 0.4, 0.6, 0.8), na.rm=T)
gdc$MLH1_quin <- factor("1st", levels=c("1st", "2nd", "3rd", "4th", "5th"))
gdc$MSH2_quin <- factor("1st", levels=c("1st", "2nd", "3rd", "4th", "5th"))
for(i in 1:length(gdc$MLH1_quin)){
  if(gdc$MLH1[[i]] < MLH1_probs[[1]]){
    gdc$MLH1_quin[[i]] <- "1st"
  } else if(gdc$MLH1[[i]] <= MLH1_probs[[2]]){
    gdc$MLH1_quin[[i]] <- "2nd"
  } else if(gdc$MLH1[[i]] <= MLH1_probs[[3]]){
    gdc$MLH1_quin[[i]] <- "3rd"
  } else if(gdc$MLH1[[i]] <= MLH1_probs[[4]]){
    gdc$MLH1_quin[[i]] <- "4th"
  } else{
    gdc$MLH1_quin[[i]] <- "5th"}
}

for(i in 1:length(gdc$MSH2_quin)){
  if(gdc$MSH2[[i]] < MSH2_probs[[1]]){
    gdc$MSH2_quin[[i]] <- "1st"
  } else if(gdc$MSH2[[i]] <= MSH2_probs[[2]]){
    gdc$MSH2_quin[[i]] <- "2nd"
  } else if(gdc$MSH2[[i]] <= MSH2_probs[[3]]){
    gdc$MSH2_quin[[i]] <- "3rd"
  } else if(gdc$MSH2[[i]] <= MSH2_probs[[4]]){
    gdc$MSH2_quin[[i]] <- "4th"
  } else{
    gdc$MSH2_quin[[i]] <- "5th"}
}

#Assign MSI factors
#per paper, high cutoff=0.4
#low cutoff = ???, trying top 5% 0.3559
MSI_prob <- quantile(gdc$MANTIS, probs=0.95, na.rm=T)
gdc$MSI <- factor("MSS", levels=c("MSI-H", 'MSI-L', "MSS"))
for (i in 1:length(gdc$MSI)){
  if(is.na(gdc$MANTIS[[i]])){
    #pass; do nothing
  }else if(gdc$MANTIS[[i]] >= 0.4){
    gdc$MSI[[i]] <- "MSI-H"
  }else if(gdc$MANTIS[[i]] >= MSI_prob){
    gdc$MSI[[i]] <- 'MSI-L'
  }
}
table(gdc$MSI, gdc$MMR)

#Assign hormone factors - HR+, HER2+
gdc$HR <- factor('HR+/HER2+', levels=c('HR+/HER2+', 'HR+/HER2-', 'HR-/HER2+', 'TNBC'))
for(i in 1:length(gdc$HR)){
  if(gdc$ER[[i]] == 'Negative' & gdc$PR[[i]] == 'Negative' & gdc$HER2[[i]] == 'Negative'){
    gdc$HR[[i]] <- 'TNBC'
  } else if((gdc$ER[[i]] == 'Positive' | gdc$PR[[i]] == 'Positive') & gdc$HER2[[i]] == 'Negative'){
    gdc$HR[[i]] <- 'HR+/HER2-'
  } else if((gdc$ER[[i]] == 'Positive' | gdc$PR[[i]] == 'Positive') & gdc$HER2[[i]] == 'Positive'){
    gdc$HR[[i]] <- 'HR+/HER2+'
  } else{
    gdc$HR[[i]] <- 'HR-/HER2+'
  }
}

#CIN factor using CX1 metric, high cutoff = ? 0.6584 is total 75 percentile
gdc$CIN_fga <- factor('NA', levels=c('CIN Low', 'CIN High', 'NA'))
gdc$CIN_cx1 <- factor('NA', levels=c('CIN Low', 'CIN High', 'NA'))
FGA_prob <- quantile(gdc$FGA, probs=0.75, na.rm=T)
CX1_prob <- quantile(gdc$CX1, probs=0.75, na.rm=T)
#FGA
for(i in 1:length(gdc$`Patient ID`)){
  if(is.na(gdc$FGA[[i]])){
    gdc$CIN_fga[[i]] <- 'NA'
  } else if(gdc$FGA[[i]] > FGA_prob){
    gdc$CIN_fga[[i]] <- 'CIN High'
  } else{
    gdc$CIN_fga[[i]] <- 'CIN Low'
  }
}
#CX1
for(i in 1:length(gdc$`Patient ID`)){
  if(is.na(gdc$CX1[[i]])){
    gdc$CIN_cx1[[i]] <- 'NA'
  } else if(gdc$CX1[[i]] > CX1_prob){
    gdc$CIN_cx1[[i]] <- 'CIN High'
  } else{
    gdc$CIN_cx1[[i]] <- 'CIN Low'
  }
}


write.csv(gdc, file="tcga_prepped.csv")