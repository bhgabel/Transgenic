library(tidyverse)
library(scales)
library(ggbiplot)

scanb <- read.csv("scanb_prepped.csv", header=T)

scanb$PAM50 <- factor(scanb$PAM50, levels=c("Basal", "Normal", "Her2", "LumA", "LumB"))

#enrichment by MMR

ggplot(data=subset(scanb, !is.na(MMR)), aes(x=MMR, fill=PAM50)) +
  geom_bar(position="fill", stat="count", color="black") +
  # geom_text(stat="count", position=position_fill(vjust=0.5),
  #           aes(label=after_stat(count)), color="white") +
  theme_classic() +
  labs(title="SCAN-B PAM50 Enrichment by MMR Loss", y="Percentage", x=NULL) +
  scale_x_discrete(labels=c("Low MLH1", "Low MSH2", "Rest")) +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  scale_y_continuous(expand=c(0,0), labels=percent) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))
ggsave("Images/SCANB_enrichment.tiff", dpi=600)

sink(file="SCANB_enrichment_stats.txt")
print("Counts")
table(scanb$PAM50, scanb$MMR)
print("Proportions")
round(proportions(table(scanb$PAM50, scanb$MMR), margin=2)*100, digits=1)
print("Overall")
chisq.test(table(scanb$PAM50, scanb$MMR))
print("MLH1 low vs rest - Basal")
fisher.test(table(scanb$PAM50 == "Basal", scanb$MMR == "MLH1"))
print("MSH2 low vs rest - Luminal A")
fisher.test(table(scanb$PAM50 == "LumA", scanb$MMR == "MSH2"))
print("MSH2 low vs rest - Luminal:")
fisher.test(table(scanb$PAM50 == "LumA" | scanb$PAM50 == "LumB", scanb$MMR == "MSH2"))
sink()
