library(tidyverse)
library(ggsignif)
library(car)
library(survival)
library(survminer)
library(forestmodel)
library(labelled)

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
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Rest", "TRUE"="Low (<12%)")) +
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
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Rest", "TRUE"="Low (<8%)")) +
  labs(x="MSH2 Gene Expression", y="Mutation Count (log10)", title="Tumor Mutation Burden - TCGA") +
  scale_fill_discrete(guide="none", palette=c("grey", "#E3E418")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2)),
        plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="n.s.", textsize=4)

gdc %>% group_by(MSH2_low) %>% summarise(mean = mean(Mutations))

ggsave(filename="Images/TCGA_MSH2_low.tiff", dpi=600)

##PAM50 subtypes, MLH1 vs MSH2 vs rest
ggplot(data=gdc, aes(x=MMR, y=Mutations, fill=PAM50)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="Mutation Count per Subtype with Low Gene Expression")

ggplot(data=subset(gdc, PAM50 != "NA"), aes(x=PAM50, y=Mutations, fill=PAM50)) +
  geom_boxplot() + scale_y_continuous(trans="log10") + theme_classic() +
  labs(title="TCGA - Mutation Count per Subtype")

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
  ylim(0.2, 0.5) + facet_grid(MMR ~ .) +
  geom_hline(yintercept=0.4, linetype="dashed")

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

#relevel factor for legend order
gdc$PAM50 <- factor(gdc$PAM50, levels=c("Basal", "Normal", "Her2", "LumA", "LumB"))

cutoff <- quantile(gdc$MLH1, probs=0.1427)
ggplot(data=subset(gdc, PAM50 != "NA"), aes(x=PAM50, y=MLH1, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  labs(title="TCGA - MLH1 Expression by Subtype", x=NULL, y="mRNA z-score") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

cutoff <- quantile(gdc$MSH2, probs=0.05)
ggplot(data=subset(gdc, PAM50 != "NA"), aes(x=PAM50, y=MSH2, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  labs(title="TCGA - MSH2 Expression by Subtype", x=NULL, y="mRNA z-score") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))
ggsave(filename="Images/TCGA_enrichment_MSH2_supplement.tiff", dpi=600)

#enrichment plot
ggplot(data=subset(gdc, PAM50 != "NA"), aes(x=MMR, fill=PAM50)) +
  geom_bar(position="fill", stat="count", color="black") +
  # geom_text(stat="count", position=position_fill(vjust=0.5),
  #           aes(label=after_stat(count)), color="white") +
  theme_classic() +
  labs(title="", y="Percentage", x=NULL) +
  scale_x_discrete(labels=c("Low MLH1", "Low MSH2", "Rest")) +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  scale_y_continuous(expand=c(0,0), labels=c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

ggsave(filename="Images/TCGA_enrichment.tiff", dpi=600)

#stats for above plots
table(gdc$PAM50, gdc$MMR)
round(proportions(table(gdc$PAM50, gdc$MMR), margin=2)*100, digits=1)
chisq.test(table(gdc$PAM50, gdc$MMR), simulate.p.value=T)
chisq.test(table(gdc$PAM50 == "Basal", gdc$MMR == "MLH1"))
chisq.test(table(gdc$PAM50 == "LumA", gdc$MMR == "MSH2"))
chisq.test(table(gdc$PAM50 == "LumA" | gdc$PAM50 == "LumB", gdc$MMR == "MSH2"))

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


#Survival plots ----
gdc$Disease.Status.l <- ifelse(gdc$Disease.Free.Status == "0:DiseaseFree", 0, 1)
gdc$Disease.Free.Months <- gdc$Disease.Free.Months/12 #years for plotting

gdc$HR <- relevel(gdc$HR, ref='HR+/HER2-')

surv.df <- gdc
surv.df$time <- ifelse(surv.df$Disease.Free.Months >=10, 10, surv.df$Disease.Free.Months)
surv.df$event <- ifelse(surv.df$Disease.Free.Months >= 10, 0, surv.df$Disease.Status.l)


#factor for combination graph
surv.df$combo <- factor(NA, levels=c("Lum:MSH2", "Lum:None",
                                     "Basal:MLH1", "Basal:None"))
for(i in 1:length(surv.df$PAM50)){
  if(is.na(surv.df$PAM50[[i]])){
    surv.df$combo[[i]] <- NA
  }
  else if(surv.df$PAM50[[i]] == "LumA" | surv.df$PAM50[[i]] == "LumB"){
    surv.df$combo[[i]] <- ifelse(surv.df$MMR[[i]] == "MSH2",
                                 "Lum:MSH2", "Lum:None")
  }
  else if(surv.df$PAM50[[i]] == "Basal"){
    surv.df$combo[[i]] <- ifelse(surv.df$MMR[[i]] == "MLH1",
                                 "Basal:MLH1", "Basal:None")
  }
  else{
    surv.df$combo[[i]] <- NA
  }
}


surv.df$combo.mlh1 <- factor(NA, levels=c("Lum:None", "Lum:MLH1",
                                          "Basal:None", "Basal:MLH1"))
for(i in 1:length(surv.df$PAM50)){
  if(is.na(surv.df$PAM50[[i]])){
    surv.df$combo.mlh1[[i]] <- NA
  }
  else if(surv.df$PAM50[[i]] == "Basal"){
    surv.df$combo.mlh1[[i]] <- ifelse(surv.df$MMR[[i]] == "MLH1",
                                      "Basal:MLH1", "Basal:None")
  }
  else if(surv.df$PAM50[[i]] == "LumA" | surv.df$PAM50[[i]] == "LumB"){
    surv.df$combo.mlh1[[i]] <- ifelse(surv.df$MMR[[i]] == "MLH1",
                                      "Lum:MLH1", "Lum:None")
  }
  else{
    surv.df$combo.mlh1[[i]] <- NA
  }
}

surv.df$combo.msh2 <- factor(NA, levels=c("Lum:None", "Lum:MSH2",
                                          "Basal:None", "Basal:MSH2"))
for(i in 1:length(surv.df$PAM50)){
  if(is.na(surv.df$PAM50[[i]])){
    surv.df$combo.msh2[[i]] <- NA
  }
  else if(surv.df$PAM50[[i]] == "Basal"){
    surv.df$combo.msh2[[i]] <- ifelse(surv.df$MMR[[i]] == "MSH2",
                                      "Basal:MSH2", "Basal:None")
  }
  else if(surv.df$PAM50[[i]] == "LumA" | surv.df$PAM50[[i]] == "LumB"){
    surv.df$combo.msh2[[i]] <- ifelse(surv.df$MMR[[i]] == "MSH2",
                                      "Lum:MSH2", "Lum:None")
  }
  else{
    surv.df$combo.msh2[[i]] <- NA
  }
}


fit <- survfit(Surv(time, event) ~ MMR, data=surv.df)
ggsurvplot(fit, data=surv.df, pval=T, xlim=c(0,10), break.time.by=1,
           surv.scale="percent", title="TCGA (Overall)")

#MLH1
fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumA"))
ggsurvplot(fit, data=subset(surv.df, PAM50=="LumA"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Luminal A)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50=="LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Luminal B)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Luminal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Her2"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Her2"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (HER2)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Basal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Basal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Basal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))


fit <- survfit(Surv(time, event) ~ combo.mlh1,
               data=surv.df)
ggsurvplot(fit, data=surv.df,
           pval=F, xlim=c(0,10), break.time.by=1, xlab="Time (Years)",
           legend.labs=c("Luminal Rest\n(n=713)", "Luminal MLH1\n(n=53)",
                         "Basal Rest\n(n=107)", "Basal MLH1\n(n=82)"),
           linetype=c(1,1,3,3), censor.shape=124, censor.size=3,
           legend.title="", surv.scale="percent",
           subtitle="(TCGA)",
           title="Probability of Survival with MLH1 Loss",
           palette=c('black', '#31688E', 'black', '#31688E'))$plot +
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold", family="Arial"),
        plot.subtitle=element_text(hjust=0.5, size=14)) +
  annotate("text", label="Luminal HR = 2.60, p = 0.006\nBasal HR = 1.36, p = 0.43",
           x=0.1, y=0.15, hjust=0)

ggsave(filename="Images/TCGA_survival_MLH1_Combo_v2.svg", dpi=600)


#MSH2
fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumA"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Luminal A)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Luminal B)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Luminal)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "Her2"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Her2"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (HER2)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "Basal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Basal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSh2"), font.legend=c(13),
           legend.title="", surv.scale="percent",
           title="TCGA (Basal)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ combo.msh2,
               data=surv.df)
ggsurvplot(fit, data=surv.df,
           pval=F, xlim=c(0,10), break.time.by=1, xlab="Time (Years)",
           legend.labs=c("Luminal Rest\n(n=731)", "Luminal MSH2\n(n=35)",
                         "Basal Rest\n(n=188)", "Basal MSH2\n(n=1)"),
           linetype=c(1,1,3,3), censor.shape=124, censor.size=3,
           legend.title="", surv.scale="percent",
           subtitle="(TCGA)",
           title="Probability of Survival with MSH2 Loss",
           palette=c('black', '#CC7722', 'black', '#CC7722'))$plot +
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold", family="Arial"),
        plot.subtitle=element_text(hjust=0.5, size=14)) +
  annotate("text", label="Luminal HR = 0.97, p = 0.96\nBasal HR = 78.8, p < 0.001",
           x=0.1, y=0.15, hjust=0)



fit <- survfit(Surv(time, event) ~ combo,
               data=surv.df)
ggsurvplot(fit, data=surv.df,
           pval=F, break.time.by=1, xlab="Time (Years)",
           legend.labs=c("Low MSH2\n(n=35)", "Rest\n(n=731)",
                         "Low MLH1\n(n=82)", "Rest\n(n=107)"),
           linetype=c(2,1,2,1), censor.shape=124, censor=F,
           legend.title="Luminal                                            Basal",
           surv.scale="percent",
           subtitle="(TCGA)",
           title="Probability of Survival with MMR Loss",
           palette=c('#c29f1f', '#c29f1f', "#56B4E9", "#56B4E9"))$plot + 
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold", family="Arial"),
        plot.subtitle=element_text(hjust=0.5, size=14, face="bold"),
        legend.title=element_text(hjust=0.5, size=12),
        legend.key.size=unit(2, "line"), legend.location="plot",
        legend.position="bottom", legend.title.position="top") +
  annotate("text", label="Luminal HR = 0.97, p = 0.96\nBasal HR = 1.36, p = 0.43",
           x=0.1, y=0.15, hjust=0)

ggsave(filename="Images/TCGA_survival_Combo.svg", dpi=600)


#Forest plots ----

#set parameters for plot
panels <- list(
  list(width = 0.03),
  list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.1, display = ~level),
  list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.55, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(
    width = 0.05,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 3, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p-value"
  ),
  list(width = 0.03)
)

ggforest(coxph(Surv(time, event) ~ MLH1_low * PAM50,
               data=as.data.frame(surv.df)))

#MLH1
surv.df$pam50.f <- factor(surv.df$PAM50, levels=c("Luminal", "Her2", "Basal"))
surv.df$pam50.f[surv.df$PAM50 == "LumA" | surv.df$PAM50 == "LumB"] <-  "Luminal"

fit_list <- vector("list", 4)
i <- 1
for(hr in c("LumA", "LumB", "Her2", "Basal")){
  fit <- coxph(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == hr))
  print(hr)
  print(fit)
  fit_list[[i]] <- fit
  names(fit_list)[[i]] <- hr
  i <- i + 1
}

forest_model(model_list=fit_list, merge_models=T, panels=panels)
ggsave(filename="Images/TCGA_forest_MLH1_v2.tiff", dpi=600)

#MSH2
fit_list <- vector("list", 3)
i <- 1
for(hr in c("LumA", "LumB", "Basal")){
  fit <- coxph(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == hr))
  print(hr)
  print(fit)
  fit_list[[i]] <- fit
  names(fit_list)[[i]] <- hr
  i <- i + 1
}

forest_model(model_list=fit_list, merge_models=T, panels=panels)
ggsave(filename="Images/TCGA_forest_MSH2_v2.tiff", dpi=600)


#other forest plotting method
#MLH1
subgroup_hr <- function(data, label) {
  s <- summary(coxph(Surv(time, event) ~ MLH1_low, data = data))
  data.frame(
    subgroup = label,
    n        = nrow(data),
    events   = s$nevent,
    hr       = s$conf.int[1, "exp(coef)"],
    lower    = s$conf.int[1, "lower .95"],
    upper    = s$conf.int[1, "upper .95"]
  )
}

forest_df <- rbind(
  subgroup_hr(surv.df, label="Overall"),
  subgroup_hr(subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), label="Luminal"),
  subgroup_hr(subset(surv.df, PAM50 == "LumA"), label="LumA"),
  subgroup_hr(subset(surv.df, PAM50 == "LumB"), label="LumB"),
  subgroup_hr(subset(surv.df, PAM50 == "Her2"), label="Her2"),
  subgroup_hr(subset(surv.df, PAM50 == "Basal"), label="Basal")
)

forest_df

forest_df$lab <- sprintf("%.2f (%.2f, %.2f)", forest_df$hr, forest_df$lower,
                         forest_df$upper)
forest_df$row <- factor(forest_df$subgroup, levels = rev(forest_df$subgroup))



ggplot(forest_df, aes(x = hr, y = row)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0.22, orientation = "y", color = "#31688E") +
  geom_point(size = 2.9, color = "#31688E") +
  geom_text(aes(x = 14, label = lab), hjust = 0, size = 4) +
  scale_x_log10(breaks = c(0, 1, 2, 4, 6, 8, 10, 12)) +
  coord_cartesian(xlim = c(0.5, 12), clip = "off") +
  labs(x = "Hazard ratio (MLH1 Low vs Rest, log scale)", y = NULL,
       title = "TCGA Forest Plot") +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(6, 130, 6, 6),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=12))

ggsave(filename="Images/TCGA_forest_MLH1_v3.tiff", dpi=600)


#MSH2
subgroup_hr <- function(data, label) {
  s <- summary(coxph(Surv(time, event) ~ MSH2_low, data = data))
  data.frame(
    subgroup = label,
    n        = nrow(data),
    events   = s$nevent,
    hr       = s$conf.int[1, "exp(coef)"],
    lower    = s$conf.int[1, "lower .95"],
    upper    = s$conf.int[1, "upper .95"]
  )
}

forest_df <- rbind(
  subgroup_hr(surv.df, label="Overall"),
  subgroup_hr(subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), label="Luminal"),
  subgroup_hr(subset(surv.df, PAM50 == "LumA"), label="LumA"),
  subgroup_hr(subset(surv.df, PAM50 == "LumB"), label="LumB"),
  subgroup_hr(subset(surv.df, PAM50 == "Basal"), label="Basal")
)

forest_df

forest_df$lab <- sprintf("%.2f (%.2f, %.2f)", forest_df$hr, forest_df$lower,
                         forest_df$upper)
forest_df$row <- factor(forest_df$subgroup, levels = rev(forest_df$subgroup))



ggplot(forest_df, aes(x = hr, y = row)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0.22, orientation = "y", color = "#31688E") +
  geom_point(size = 2.9, color = "#31688E") +
  geom_text(aes(x = 14, label = lab), hjust = 0, size = 4) +
  scale_x_log10(breaks = c(0, 1, 2, 4, 6, 8, 10, 12)) +
  coord_cartesian(xlim = c(0.5, 12), clip = "off") +
  labs(x = "Hazard ratio (MSH2 Low vs Rest, log scale)", y = NULL,
       title = "TCGA Forest Plot") +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(6, 130, 6, 6),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=12))

ggsave(filename="Images/TCGA_forest_MLH1_v3.tiff", dpi=600)