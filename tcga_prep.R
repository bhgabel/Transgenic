library(readr)
library(readxl)
library(tidyverse)
library(TCGAbiolinks)

#Creating and subsetting data frames ----
#mrna z-scores gene expression
mrna <- read_delim("TCGA mRNA expression (RNA Seq V2 RSEM).txt",
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(-c("STUDY_ID"))

#Clinical info
clinical <- read_delim("brca_tcga_clinical_data.tsv", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(c("Patient ID", "Sample ID", "Mutation Count", "Fraction Genome Altered",
                  "Disease Free (Months)","Disease Free Status",
                  "ER Status By IHC", "IHC-HER2", "PR status by ihc"))

colnames(clinical)[3:9] <- c('Mutations', 'FGA',
                             'Disease.Free.Months', 'Disease.Free.Status',
                             'ER', 'HER2', 'PR')

#Combine into one data frame, rename for pam50 subtyping
gdc <- right_join(clinical, mrna, by=c("Sample ID" = "SAMPLE_ID"))
#rename(all_of(c("CDCA1" = "NUF2", "KNTC2" = "NDC80", "ORC6L" = "ORC6")))

#Add MSI MANTIS scores
mantis <- read_excel("MSI_Mantis.xlsx") %>%
  dplyr::select(c("Case ID", "MANTIS Score"))
colnames(mantis)[2] <- "MANTIS"

gdc <- left_join(gdc, mantis, by=c("Patient ID" = "Case ID"))


#Pull pam50 subtypes from TCGA package
subtypes <- TCGAquery_subtype("brca") %>% dplyr::select(c("patient", "BRCA_Subtype_PAM50"))
colnames(subtypes)[2] <- "PAM50"

#Combine to main data frame
gdc <- left_join(gdc, subtypes, by=c("Patient ID" = "patient"))

#remove primary data frames, since everything joined into gdc
remove(list = c('CIN_data', 'clinical', 'mantis', 'mrna', 'mrna_pam50', 'subtypes'))

#Assign factors ----
gdc <- drop_na(gdc, MLH1)
#Assign factors for low gene expression and remove combined loss
#trying new cutoff
gdc$MLH1_low <- (gdc$MLH1 <= quantile(gdc$MLH1, probs=0.12, na.rm=T))
gdc$MSH2_low <- (gdc$MSH2 <= quantile(gdc$MSH2, probs=0.08, na.rm=T))
gdc$MMR <- factor("None", levels=c("MLH1", "MSH2", "None"))
for(i in 1:length(gdc$MMR)){
  if(is.na(gdc$MLH1_low[[i]])){
    gdc$MMR[[i]] <- NA
  } else if(gdc$MLH1_low[[i]] && gdc$MSH2_low[[i]]){
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
gdc$HR <- factor('HR+/HER2-', levels=c('HR+/HER2-', 'HR+/HER2+', 'HR-/HER2+', 'TNBC'))
for(i in 1:length(gdc$HR)){
  if(is.na(gdc$ER[[i]]) | is.na(gdc$PR[[i]]) | is.na(gdc$HER2[[i]])){
      gdc$HR[[i]] <- NA
  } else if(gdc$ER[[i]] == 'Negative' & gdc$PR[[i]] == 'Negative' & gdc$HER2[[i]] == 'Negative'){
      gdc$HR[[i]] <- 'TNBC'
  } else if((gdc$ER[[i]] == 'Positive' | gdc$PR[[i]] == 'Positive') & gdc$HER2[[i]] == 'Negative'){
      gdc$HR[[i]] <- 'HR+/HER2-'
  } else if((gdc$ER[[i]] == 'Positive' | gdc$PR[[i]] == 'Positive') & gdc$HER2[[i]] == 'Positive'){
      gdc$HR[[i]] <- 'HR+/HER2+'
  } else if(gdc$ER[[i]] == 'Negative' & gdc$PR[[i]] == 'Negative' & gdc$HER2[[i]] == 'Positive'){
      gdc$HR[[i]] <- 'HR-/HER2+'
  } else {
      gdc$HR[[i]] <- NA
  }
}


write.csv(gdc, file="tcga_prepped.csv")
