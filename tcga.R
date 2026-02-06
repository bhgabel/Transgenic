# pam50 probes
"""ACTR3B, ANLN, BAG1, BCL2, BIRC5, BLVRA, CCNB1, CCNE1, CDC20, CDC6, CDCA1, 
CDH3, CENPF, CEP55, CXXC5, EGFR, ERBB2, ESR1, EXO1, FGFR4, FOXA1, FOXC1, 
GPR160, GRB7, KIF2C, KNTC2, KRT14, KRT17, KRT5, MAPT, MDM2, MELK, MIA, 
MKI67, MLPH, MMP11, MYBL2, MYC, NAT1, ORC6L, PGR, PHGDH, PTTG1, RRM2, 
SFRP1, SLC39A6, TMEM45B, TYMS, UBE2C, UBE2T
CDCA1 -> NUF2, KNTC2 -> NDC80, ORC6L -> ORC6"""

library(tidyverse)
library(ggsignif)

gdc <- read.csv("tcga_prepped.csv")

#Statistical tests ----
pairwise.wilcox.test(gdc$Mutations, gdc$pam50_f)
pairwise.wilcox.test(gdc$Mutations, gdc$MLH1_quin)
pairwise.wilcox.test(gdc$Mutations, gdc$MSH2_quin)
pairwise.wilcox.test(gdc$Mutations, gdc$MSI)

wilcox.test(gdc$Mutations[gdc$CIN_cx1 == 'CIN High'], gdc$Mutations[gdc$CIN_cx1 == 'CIN Low'])
wilcox.test(gdc$Mutations[gdc$CIN_fga == 'CIN High'], gdc$Mutations[gdc$CIN_fga == 'CIN Low'])

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

#assumptions not met
cor.test(gdc$Mutations, gdc$CX1) #sig
cor.test(gdc$Mutations, gdc$FGA) #not sig


#Plots ----
##Boxplots ----
##MLH1 low vs rest
ggplot(data=gdc, aes(x=MLH1_low, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Normal", "TRUE"="Low (<12%)")) +
  labs(x="MLH1 Gene Expression", y="Mutation Count (log10)") +
  scale_fill_discrete(guide="none", palette=c("grey", "steelblue")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2))) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="***", textsize=5)

#fold change
gdc %>% group_by(MLH1_low) %>% summarise(mean = mean(Mutations))

#ggsave(filename="Images/MLH1_low.tiff", height=5, width=5)

#MSH2 low vs rest
ggplot(data=gdc, aes(x=MSH2_low, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Normal", "TRUE"="Low (<8%)")) +
  labs(x="MSH2 Gene Expression", y="Mutation Count (log10)") +
  scale_fill_discrete(guide="none", palette=c("grey", "yellow")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2))) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="n.s.", textsize=4)

ggsave(filename="Images/test.tiff", height=5, width=5)

##PAM50 subtypes, MLH1 vs MSH2 vs rest
ggplot(data=gdc, aes(x=MMR, y=Mutations, fill=PAM50)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Mutation Count per Subtype with Low Gene Expression")

ggplot(data=gdc, aes(x=MMR, y=Mutations, fill=pam50_f)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Mutation Count per Subtype with Low Gene Expression")

ggplot(data=gdc, aes(x=pam50_f, y=Mutations, fill=MMR)) +
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
ggplot(data=subset(gdc, CIN_fga != 'NA'), aes(x=CIN_fga, y=Mutations, fill=CIN_fga)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="CIN + TMB")

ggplot(data=subset(gdc, CIN_cx1 != 'NA'), aes(x=CIN_cx1, y=Mutations, fill=CIN_cx1)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + labs(title="CIN + TMB")

##Stacked Column ----
ggplot(data=subset(gdc, MMR != "None"), aes(x=HR, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=subset(gdc, MMR != "None" & PAM50 != "NA"), aes(x=MMR, fill=PAM50)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")
table(gdc$PAM50, gdc$MMR)

ggplot(data=gdc, aes(x=HR, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=gdc, aes(x=MSI, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=gdc, aes(x=pam50_f, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="MSI Status by Subtype")

ggplot(data=gdc, aes(x=HR, fill=MSI)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="MSI Classification")


#MSI by subtype
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

ggplot(data=gdc, aes(x=pam50_f, fill=MSI)) +
  geom_bar(position="fill", stat="count") + facet_grid(MMR ~ .) +
  theme_classic() + labs(title="MSI by HR and MMR low")

#CIN
ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=CIN_f, fill=MMR)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=HR, fill=CIN_f)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Low mRNA Expression")

ggplot(data=subset(gdc, CIN_f != 'NA'), aes(x=HR, fill=CIN_f)) +
  geom_bar(position="fill", stat="count") + facet_grid(MMR ~ .) +
  theme_classic() + labs(title="CIN by HR and MMR low")

# mutations < 100 just to visualize data since range too wide
ggplot(data=subset(gdc, Mutations <= 100 & !is.na(CX1)), aes(x=CX1, y=Mutations)) +
  geom_point() + geom_smooth(method='lm', se=F)

ggplot(data=subset(gdc, Mutations <= 100), aes(x=FGA, y=Mutations)) +
  geom_point() + geom_smooth(method='lm', se=F)

ggplot(data=subset(gdc, !is.na(CX1)), aes(x=CX1, y=FGA)) + geom_point() + geom_smooth(method='lm', se=F)


aneuploidy_score <- read_excel("aneuploidy_score.xlsx")
cin <- left_join(gdc, aneuploidy_score, by=c("Sample ID"="Sample"))
ggplot(data=cin, aes(x=FGA, y=`AneuploidyScore(AS)`)) + geom_point() + geom_smooth(method='lm', se=F)
ggplot(data=cin, aes(x=`AneuploidyScore(AS)`, y=FGA)) + geom_point() + geom_smooth(method='lm', se=F)
ggplot(data=subset(cin, Mutations <= 100), aes(x=`AneuploidyScore(AS)`, y=Mutations)) + geom_point() + geom_smooth(method='lm', se=F)

cor.test(cin$Mutations, cin$`AneuploidyScore(AS)`)
