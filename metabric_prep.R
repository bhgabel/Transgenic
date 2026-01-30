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

write.csv("metabric_prepped.csv")