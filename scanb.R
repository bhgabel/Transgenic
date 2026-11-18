library(tidyverse)
library(scales)
library(ggbiplot)

scanb <- read.csv("scanb_prepped.csv", header=T)

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


#Gene expression by PAM50
write.csv(dplyr::select(scanb, 1, MLH1, MSH2, PAM50), file="scanb_pam50.csv", row.names=F)

#MLH1
cutoff <- quantile(scanb$MLH1, probs=0.12)

ggplot(data=scanb, aes(x=PAM50, y=MLH1, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  coord_cartesian(ylim=c(0, 75)) +
  labs(title="SCAN-B - MLH1 Expression by Subtype", x=NULL, y="FPKM") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

scanb %>% group_by(PAM50) %>% summarise(mean = mean(MLH1))

#MSH2
cutoff <- quantile(scanb$MSH2, probs=0.08)
ggplot(data=scanb, aes(x=PAM50, y=MSH2, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  coord_cartesian(ylim=c(0, 75)) +
  labs(title="TCGA - MSH2 Expression by Subtype", x=NULL, y="FPKM") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

scanb %>% group_by(PAM50) %>% summarise(mean = mean(MSH2))
