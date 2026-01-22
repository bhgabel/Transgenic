# pam50 probes
"""ACTR3B, ANLN, BAG1, BCL2, BIRC5, BLVRA, CCNB1, CCNE1, CDC20, CDC6, CDCA1, 
CDH3, CENPF, CEP55, CXXC5, EGFR, ERBB2, ESR1, EXO1, FGFR4, FOXA1, FOXC1, 
GPR160, GRB7, KIF2C, KNTC2, KRT14, KRT17, KRT5, MAPT, MDM2, MELK, MIA, 
MKI67, MLPH, MMP11, MYBL2, MYC, NAT1, ORC6L, PGR, PHGDH, PTTG1, RRM2, 
SFRP1, SLC39A6, TMEM45B, TYMS, UBE2C, UBE2T
CDCA1 -> NUF2, KNTC2 -> NDC80, ORC6L -> ORC6"""

library(readr)
library(readxl)
library(tidyverse)
library(TCGAbiolinks)

# data("pam50")
# data("pam50.robust")

#Creating and subsetting data frames ----
#FPKM z-scores gene expression
mrna <- read_delim("mRNA expression Zscores.txt", 
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
       dplyr::select(-c('STUDY_ID', 'MSH2', 'MLH1'))

#Clinical info
clinical <- read_delim("brca_tcga_gdc_clinical_data.tsv", 
                          delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
            dplyr::select(c("Patient ID", "Sample ID", "Mutation Count", "Fraction Genome Altered"))

#HR status                          
firehose <- read_excel("~/Documents/Haricharan/PanCancer/TCGA_Firehose.xlsx", sheet = "Clinical") %>% 
            dplyr::select(c("Sample ID", "ER Status By IHC", "IHC-HER2", "PR status by ihc"))

#Combine into one data frame, rename for pam50 subtyping
gdc <- right_join(clinical, mrna_pam50, by=c("Sample ID" = "SAMPLE_ID"))
       #rename(all_of(c("CDCA1" = "NUF2", "KNTC2" = "NDC80", "ORC6L" = "ORC6")))
gdc <- left_join(gdc, mrna, by=c("Sample ID" = "SAMPLE_ID"))

#Add MSI MANTIS scores
mantis <- read_excel("MSI_Mantis.xlsx") %>%
          dplyr::select(c("Case ID", "MANTIS Score"))

gdc <- left_join(gdc, mantis, by=c("Patient ID" = "Case ID"))

gdc <- left_join(gdc, firehose)

#Add CIN data
CIN_data <- read_excel("CIN_data.xlsx", sheet = "ST_20_TCGA_Activities_scaled")
colnames(CIN_data)[1] <- 'Patient ID'

gdc <- left_join(gdc, CIN_data[c('Patient ID', 'Total', 'CX1', 'CX6')])

#remove primary data frames, since everything joined into gdc
remove(list = c('CIN_data', 'clinical', 'firehose', 'mantis', 'mrna', 'mrna_pam50', 'subtypes'))


# gene_info <- data.frame(colnames(firehose)[-c(1:5)], pam50$centroids.map$EntrezGene.ID)
# colnames(gene_info) <- c("probe", "EntrezGene.ID")
# 
# pam50_pred <- molecular.subtyping(sbt.model="pam50", data=firehose, annot=gene_info, do.mapping=T, verbose=T)

#Pull pam50 subtypes from TCGA package
subtypes <- TCGAquery_subtype("brca") %>% dplyr::select(c("patient", "BRCA_Subtype_PAM50"))

#Combine to main data frame
gdc <- left_join(gdc, subtypes, by=c("Patient ID" = "patient"))

#Rename columns
colnames(gdc) <- c('Patient ID', 'Sample ID', 'Mutations', 'FGA', 'MLH1', 'MSH2', 'p53', 'ERBB2', 'ESR1', 'MANTIS', 'ER', 'HER2', 'PR', 'CIN', 'CX1', 'CX6', 'PAM50')

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
gdc$MSI <- factor("MSS", levels=c("MSI-H", 'MSI-L', "MSS"))
for (i in 1:length(gdc$MSI)){
  if(is.na(gdc$MANTIS[[i]])){
    #pass; do nothing
  }else if(gdc$MANTIS[[i]] >= 0.4){
    gdc$MSI[[i]] <- "MSI-H"
  }else if(gdc$MANTIS[[i]] >= 0.3559){
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

#CIN factor, high cutoff = ? 1.97 is total 75 percentile
gdc$CIN_f <- factor('NA', levels=c('CIN High', 'CIN Low', 'NA'))
for(i in 1:length(gdc$CIN)){
  if(is.na(gdc$CIN[[i]])){
    gdc$CIN_f[[i]] <- 'NA'
  } else if(gdc$CIN[[i]] > 1.97){
    gdc$CIN_f[[i]] <- 'CIN High'
  } else{
    gdc$CIN_f[[i]] <- 'CIN Low'
  }
}



#T-tests for non-normal data ----
pairwise.wilcox.test(gdc$Mutations, gdc$pam50_f)
pairwise.wilcox.test(gdc$Mutations, gdc$MLH1_quin)
pairwise.wilcox.test(gdc$Mutations, gdc$MSH2_quin)
pairwise.wilcox.test(gdc$Mutations, gdc$MSI)

wilcox.test(gdc$Mutations[gdc$MLH1_low], gdc$Mutations[!gdc$MLH1_low])
wilcox.test(gdc$Mutations[gdc$MSH2_low], gdc$Mutations[!gdc$MSH2_low])

#mlh.df <- data.frame(group1="FALSE", group2="TRUE", p=p1$p.value, y.position=1500)

kruskal.test(data=gdc, Mutations ~ MLH1_quin)
kruskal.test(data=gdc, Mutations ~ MSH2_quin)
kruskal.test(data=gdc, Mutations ~ MSI)

kruskal.test(data=subset(gdc, pam50_f=="Basal"), Mutations ~ MMR)
kruskal.test(data=subset(gdc, pam50_f=="Luminal"), Mutations ~ MMR)
kruskal.test(data=subset(gdc, PAM50=="LumA"), Mutations ~ MMR)
kruskal.test(data=subset(gdc, PAM50=="LumB"), Mutations ~ MMR)

wilcox.test(gdc$Mutations[gdc$MLH1_low & gdc$pam50_f=="Luminal"],
            gdc$Mutations[!gdc$MLH1_low & gdc$pam50_f=="Luminal"])
wilcox.test(gdc$Mutations[gdc$MSH2_low & gdc$pam50_f=="Luminal"],
            gdc$Mutations[!gdc$MSH2_low & gdc$pam50_f=="Luminal"])

wilcox.test(gdc$Mutations[gdc$MLH1_low & gdc$pam50_f=="Basal"],
            gdc$Mutations[!gdc$MLH1_low & gdc$pam50_f=="Basal"])
wilcox.test(gdc$Mutations[gdc$MSH2_low & gdc$pam50_f=="Basal"],
            gdc$Mutations[!gdc$MSH2_low & gdc$pam50_f=="Basal"])

wilcox.test(gdc$Mutations[gdc$MLH1_low & gdc$pam50_f=="Other"],
            gdc$Mutations[!gdc$MLH1_low & gdc$pam50_f=="Other"])
wilcox.test(gdc$Mutations[gdc$MSH2_low & gdc$pam50_f=="Other"],
            gdc$Mutations[!gdc$MSH2_low & gdc$pam50_f=="Other"])

wilcox.test(gdc$Mutations[gdc$MSI=="MSS"], gdc$Mutations[gdc$MSI=="MSI-H"])
wilcox.test(gdc$Mutations[gdc$MSI=="MSS" & gdc$MLH1_low],
            gdc$Mutations[gdc$MSI=="MSI-H" & gdc$MLH1_low])
wilcox.test(gdc$Mutations[gdc$MSI=="MSS" & !gdc$MLH1_low],
            gdc$Mutations[gdc$MSI=="MSI-H" & !gdc$MLH1_low])
wilcox.test(gdc$Mutations[gdc$MSI=="MSI-H" & gdc$MLH1_low],
            gdc$Mutations[gdc$MSI=="MSI-H" & !gdc$MLH1_low])
wilcox.test(gdc$Mutations[gdc$MSI=="MSS" & gdc$MLH1_low],
            gdc$Mutations[gdc$MSI=="MSS" & !gdc$MLH1_low])

wilcox.test(gdc$Mutations[gdc$MSI=="MSS" & gdc$MSH2_low],
            gdc$Mutations[gdc$MSI=="MSI-H" & gdc$MSH2_low])
wilcox.test(gdc$Mutations[gdc$MSI=="MSS" & !gdc$MSH2_low],
            gdc$Mutations[gdc$MSI=="MSI-H" & !gdc$MSH2_low])
wilcox.test(gdc$Mutations[gdc$MSI=="MSI-H" & gdc$MSH2_low],
            gdc$Mutations[gdc$MSI=="MSI-H" & !gdc$MSH2_low])
wilcox.test(gdc$Mutations[gdc$MSI=="MSS" & gdc$MSH2_low],
            gdc$Mutations[gdc$MSI=="MSS" & !gdc$MSH2_low])

wilcox.test(gdc$Mutations[gdc$CIN_f=='CIN High'], gdc$Mutations[gdc$CIN_f=='CIN Low'])

#chi square tests for stacked columns
table(gdc$MSI, gdc$MMR)
chisq.test(table(gdc$MSI, gdc$MMR))
chisq.test(c(28,20,172), p=c(69/701, 44/701, 588/701)) #p=0.06435, MSI-L vs MSS


#Plots ----
##Boxplots ----
##MLH1 low vs rest
ggplot(data=gdc, aes(x=MLH1_low, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="MLH1 Low (<12%)")
#p + add_pvalue(mlh.df)

#MSH2 low vs rest
ggplot(data=gdc, aes(x=MSH2_low, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="MSH2 Low (<8%)")

##PAM50 subtypes, MLH1 vs MSH2 vs rest
ggplot(data=gdc, aes(x=MMR, y=Mutations, fill=PAM50)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Mutation Count per Subtype with Low Gene Expression")

ggplot(data=gdc, aes(x=MMR, y=Mutations, fill=pam50_f)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Mutation Count per Subtype with Low Gene Expression")

ggplot(data=gdc, aes(x=HR, y=Mutations, fill=MMR)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Mutation Count per Subtype with Low Gene Expression")

###Basal
ggplot(data=subset(gdc, gdc$PAM50=="Basal"),
       aes(x=MMR, y=Mutations, fill=PAM50)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

###Luminal A + Luminal B
ggplot(data=subset(gdc, gdc$PAM50=="LumA" | gdc$PAM50=="LumB"),
       aes(x=MMR, y=Mutations, fill=PAM50)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

##Quintiles
ggplot(data=gdc, aes(x=MLH1_quin, y=Mutations, fill=MLH1_quin)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

ggplot(data=gdc, aes(x=MSH2_quin, y=Mutations, fill=MSH2_quin)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

ggplot(data=subset(gdc, pam50_f=="Luminal"), aes(x=MLH1_quin, y=Mutations, fill=MLH1_quin)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Luminal - Quintiles")

ggplot(data=subset(gdc, pam50_f=="Luminal"), aes(x=MSH2_quin, y=Mutations, fill=MSH2_quin)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Luminal - Quintiles")

ggplot(data=subset(gdc, pam50_f=="Basal"), aes(x=MLH1_quin, y=Mutations, fill=MLH1_quin)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Basal - Quintiles")

ggplot(data=subset(gdc, pam50_f=="Basal"), aes(x=MSH2_quin, y=Mutations, fill=MSH2_quin)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Basal - Quintiles")

##MLH1 / MSH2 low per subtype
ggplot(data=gdc, aes(x=pam50_f, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()
ggplot(data=gdc, aes(x=pam50_f, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

##MSI
ggplot(data=gdc, aes(x=MSI, y=Mutations, fill=MSI)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="MSI + TMB")

ggplot(data=gdc, aes(x=MSI, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="MSI + TMB")

ggplot(data=gdc, aes(x=MSI, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="MSI + TMB")

ggplot(data=gdc, aes(x=MSI, y=Mutations, fill=pam50_f)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="MSI + TMB")

#CIN
ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=CIN_f, y=Mutations, fill=CIN_f)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="CIN + TMB")

##Stacked Column ----
ggplot(data=subset(gdc, MMR != "None"), aes(x=HR, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=gdc, aes(x=HR, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=gdc, aes(x=MSI, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=gdc, aes(x=HR, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="MSI Classification")

ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=CIN_f, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=HR, fill=CIN_f)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=HR, fill=CIN_f)) +
  geom_bar(position="fill", stat="count") + facet_grid(MMR ~ .) +
  theme_classic() + labs(title="CIN by HR and MMR low")

#MSI by HR subtype
ggplot(data=subset(gdc, HR == 'HR+/HER2+'), aes(x=MMR, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="HR+/HER2+")

ggplot(data=subset(gdc, HR == 'HR+/HER2-'), aes(x=MMR, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="HR+/HER2-")

ggplot(data=subset(gdc, HR == 'HR-/HER2+'), aes(x=MMR, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="HR-/HER2+")

ggplot(data=subset(gdc, HR == 'TNBC'), aes(x=MMR, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="TNBC")

ggplot(data=gdc, aes(x=HR, fill=MSI)) +
  geom_bar(position="fill", stat="count") + facet_grid(MMR ~ .) +
  theme_classic() + labs(title="MSI by HR and MMR low")

#CIN scatter plots
# mutations < 100 just to visualize data since range too wide
ggplot(data=subset(gdc, Mutations <= 100), aes(x=CX1, y=Mutations)) +
  geom_point()

ggplot(data=subset(gdc, Mutations <= 100), aes(x=FGA, y=Mutations)) +
  geom_point()