library(tidyverse)
library(ggsignif)
library(car)

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

#TNBC pt, MLH1 low vs not
wilcox.test(gdc$Mutations[gdc$MLH1_low & gdc$HR=="TNBC"],
            gdc$Mutations[!gdc$MLH1_low & gdc$HR=="TNBC"])
wilcox.test(gdc$MANTIS[gdc$MLH1_low & gdc$HR=="TNBC"],
            gdc$MANTIS[!gdc$MLH1_low & gdc$HR=="TNBC"])

wilcox.test(gdc$Mutations[gdc$CIN_f=='CIN High'], gdc$Mutations[gdc$CIN_f=='CIN Low'])

#chi square tests for stacked columns
table(gdc$MSI, gdc$MMR)
chisq.test(table(gdc$MSI, gdc$MMR))
chisq.test(c(28,20,172), p=c(69/701, 44/701, 588/701)) #p=0.06435, MSI-L vs MSS

#assumptions not met
cor.test(gdc$Mutations, gdc$CX1) #sig
cor.test(gdc$Mutations, gdc$FGA) #not sig

mlh1.lm <- lm(data=subset(gdc, MMR == "MLH1"), log(Mutations) ~ HR)
mlh1.glm <- glm(data=subset(gdc, MMR == "MLH1"), Mutations ~ HR, family=quasipoisson)
plot(mlh1.lm); anova(mlh1.lm)
kruskal.test(data=subset(gdc, MMR == "MLH1"), Mutations ~ HR)

msi.lm <- lm(data=subset(gdc, MMR == "MLH1"), MANTIS ~ HR)
msi.glm <- glm(data=subset(gdc, MMR == "MLH1"), Mutations ~ HR, family=quasipoisson)
plot(msi.lm); anova(msi.lm)
kruskal.test(data=subset(gdc, MMR == "MLH1"), Mutations ~ HR)


#Plots ----
##Boxplots ----
##MLH1 low vs rest
ggplot(data=gdc, aes(x=MLH1_low, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Normal", "TRUE"="Low (<12%)")) +
  labs(x="MLH1 Gene Expression", y="Mutation Count (log10)", title="Tumor Mutation Burden - TCGA") +
  scale_fill_discrete(guide="none", palette=c("grey", "steelblue")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2)),
        plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="***", textsize=5)

#fold change
gdc %>% group_by(MLH1_low) %>% summarise(mean = mean(Mutations))

ggsave(filename="Images/TCGA_MLH1_low.tiff", dpi=600)

#MSH2 low vs rest
ggplot(data=gdc, aes(x=MSH2_low, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Normal", "TRUE"="Low (<8%)")) +
  labs(x="MSH2 Gene Expression", y="Mutation Count (log10)", title="Tumor Mutation Burden - TCGA") +
  scale_fill_discrete(guide="none", palette=c("grey", "yellow3")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2)),
        plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="n.s.", textsize=4)

ggsave(filename="Images/TCGA_MSH2_low.tiff", dpi=600)

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
ggplot(data=gdc, aes(x=HR, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()
ggplot(data=gdc, aes(x=pam50_f, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

ggplot(data=subset(gdc, MMR=="MLH1"), aes(x=HR, y=Mutations, fill=HR)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic()

##MSI
ggplot(data=gdc, aes(x=MMR, y=MANTIS, fill=MMR)) +
  geom_boxplot() + theme_classic() + labs(title="MANTIS Scores") +
  ylim(0, 0.5) + geom_hline(yintercept=0.4, linetype="dashed")

ggplot(data=gdc, aes(x=HR, y=MANTIS, fill=HR)) +
  geom_boxplot() + theme_classic() + labs(title="MANTIS Scores") +
  ylim(0, 0.6) + geom_hline(yintercept=0.4, linetype="dashed")

ggplot(data=gdc, aes(x=MMR, y=MANTIS, fill=HR)) +
  geom_boxplot() + theme_classic() + labs(title="MANTIS Scores") +
  ylim(0, 0.6) + geom_hline(yintercept=0.4, linetype="dashed")

#get basic stats for above plots
table(gdc$HR, gdc$MMR)
gdc %>% dplyr::group_by(MMR, HR) %>%
  dplyr::summarise(mean = mean(Mutations, na.rm=T), .groups="drop") %>%
  pivot_wider(names_from=MMR, values_from=mean)

gdc %>% group_by(MMR, HR) %>%
  dplyr::summarise(mean = mean(MANTIS, na.rm=T), .groups="drop") %>%
  pivot_wider(names_from=MMR, values_from=mean)

ggplot(data=gdc, aes(x=HR, y=MANTIS)) +
  geom_boxplot() + theme_classic() + labs(title="MANTIS Scores") +
  ylim(0.2, 0.5) + facet_grid(MMR ~ .) + geom_hline(yintercept=0.4, linetype="dashed")

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

