library(tidyverse)
library(viridisLite)
library(ggsignif)
library(survival)
library(survminer)
library(forestmodel)

meta <- read.csv("metabric_prepped.csv")

color_blind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                         "#0072B2", "#D55E00", "#CC79A7", "#999999")



#Plots ----
##Boxplots ----
##MLH1 low vs rest
wilcox.test(Mutations ~ MLH1_low, data=meta)
wilcox.test(meta$Mutations[meta$MLH1_low], meta$Mutations[!meta$MLH1_low])

ggplot(data=meta, aes(x=MLH1_low, y=Mutations, fill=MLH1_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Normal", "TRUE"="Low (<12%)")) +
  labs(x="MLH1 Gene Expression", y="Mutation Count (log10)", title="Tumor Mutation Burden - Metabric") +
  scale_fill_discrete(guide="none", palette=c("grey", "steelblue")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2)),
        plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="*", textsize=5)

meta %>% group_by(MLH1_low) %>% summarise(mean = mean(Mutations, na.rm=T))

ggsave(filename="Images/Metabric_MLH1_low.tiff", dpi=600)

#MSH2 low vs rest
wilcox.test(meta$Mutations[meta$MSH2_low], meta$Mutations[!meta$MSH2_low])

ggplot(data=meta, aes(x=MSH2_low, y=Mutations, fill=MSH2_low)) +
  geom_boxplot() + scale_y_continuous(trans="log10") +
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Normal", "TRUE"="Low (<8%)")) +
  labs(x="MSH2 Gene Expression", y="Mutation Count (log10)", title="Tumor Mutation Burden - Metabric") +
  scale_fill_discrete(guide="none", palette=c("grey", "yellow3")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.2)),
        plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  geom_signif(comparisons=list(c("FALSE","TRUE")), annotations="n.s.", textsize=4)

meta %>% group_by(MSH2_low) %>% summarise(mean = mean(Mutations, na.rm=T))

ggsave(filename="Images/Metabric_MSH2_low.tiff", dpi=600)


#relevel factor for legend order
meta$PAM50 <- factor(meta$PAM50, levels=c("Basal", "Normal", "Her2", "LumA", "LumB"))

#stacked columns - enrichment plots
ggplot(data=subset(meta, PAM50 != "NC" & PAM50 != "claudin-low"), aes(x=MMR, fill=PAM50)) +
  geom_bar(position="fill", stat="count") +
  theme_classic() + labs(title="Metabric - PAM50 Subtype by Mismatch Repair Deficiency", y="Proportion") +
  scale_x_discrete(labels=c("Low MLH1", "Low MSH2", "Rest"))

percents <- data.frame(MMR=c(rep("MLH1", 5), rep("MSH2", 5), rep("None", 5)),
                       PAM50=c(rep(c("Basal", "Normal", "Her2", "LumA", "LumB"), 3)),
                       percent=c('37.2%', '2.94%', '19.6%', '10.3%', '29.9%', 
                                 '10.1%', '3.88%', '20.2%', '39.5%', '26.4%',
                                 '10.4%', '2.26%', '11.2%', '34.7%', '41.5%'))

#enrichment plots
cutoff <- quantile(meta$MLH1, probs=0.1283)
ggplot(data=subset(meta, PAM50 != "NA"), aes(x=PAM50, y=MLH1, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  labs(title="Metabric - MLH1 Expression by Subtype", x=NULL, y="mRNA z-score") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

cutoff <- quantile(meta$MSH2, probs=0.0106)
ggplot(data=subset(meta, PAM50 != "NA"), aes(x=PAM50, y=MSH2, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  labs(title="Metabric - MSH2 Expression by Subtype", x=NULL, y="mRNA z-score") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))
ggsave(filename="Images/Metabric_enrichment_MLH1_supplement.tiff", dpi=600)

ggplot(data=subset(meta, PAM50 != "NA"), aes(x=MMR, fill=PAM50)) +
  geom_bar(position="fill", stat="count", color="black") +
  # geom_text(stat="count", position=position_fill(vjust=0.5),
  #           aes(label=after_stat(count)), color="white") +
  theme_classic() + labs(title="Metabric",
                         y="Percentage", x=NULL) +
  scale_x_discrete(labels=c("Low MLH1", "Low MSH2", "Rest")) +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  scale_y_continuous(expand=c(0,0), labels=c("0%", "25%", "50%", "75%", "100%")) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

ggsave(filename="Images/Metabric_enrichment_v2.svg", dpi=600)

#stats for above plots
table(meta$PAM50, meta$MMR)
round(proportions(table(meta$PAM50, meta$MMR), margin=2)*100, digits=1)
chisq.test(table(meta$PAM50, meta$MMR), simulate.p.value=T)
chisq.test(table(meta$PAM50 == "Basal", meta$MMR == "MLH1"))
chisq.test(table(meta$PAM50 == "LumA", meta$MMR == "MSH2"))
chisq.test(table(meta$PAM50 == "LumA" | meta$PAM50 == "LumB", meta$MMR == "MSH2"))

#Survival plots ----
#colnames(meta)[c(61,62)] <- c("Surv.Status", "Relapse.Months")
meta$Relapse.Status.l <- ifelse(meta$Relapse.Status == "0:Not Recurred", 0, 1)
meta$Relapse.Months <- meta$Relapse.Months/12 #years for plotting

meta$HR <- relevel(meta$HR, ref='HR+/HER2-')
meta$MMR <- relevel(meta$MMR, ref="None")


surv.df <- meta
surv.df$time <- ifelse(surv.df$Relapse.Months >=10, 10, surv.df$Relapse.Months)
surv.df$event <- ifelse(surv.df$Relapse.Months >= 10, 0, surv.df$Relapse.Status.l)

# surv.df <- read.csv("Relapse_Free_METABRIC_combined_Sv_070626.csv", header=T)
# surv.df$time <- surv.df$time/12

fit <- survfit(Surv(time, event) ~ MMR,
               data=surv.df)
ggsurvplot(fit, data=surv.df, pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Overall)")

#MLH1

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumA"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal A)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal B)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))


fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Luminal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Luminal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Her2"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Her2"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (HER2)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Basal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Basal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Basal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

ggsave(filename="Images/Metabric_survival_MLH1_Basal.tiff", dpi=600)


#MSH2
fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumA"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal A)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal B)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Luminal)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "Her2"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Her2"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (HER2)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "Basal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Basal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Normal", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Relapse probability",
           title="Metabric (Basal)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

ggsave(filename="Images/Metabric_survival_MSH2_Basal.tiff", dpi=600)


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

#MLH1
surv.df$PAM50 <- factor(surv.df$PAM50, levels=c("Luminal", "Her2", "Basal"))

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

forest_model(model_list=fit_list, panels=panels, merge_models=T)
ggsave(filename="Images/Metabric_forest_MLH1.tiff", dpi=600)

#MSH2
#zero events in MSH2 subset for LumA and Basal
fit_list <- vector("list", 2)
i <- 1
for(hr in c("LumB", "Her2")){
  fit <- coxph(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == hr))
  print(hr)
  print(fit)
  fit_list[[i]] <- fit
  names(fit_list)[[i]] <- hr
  i <- i + 1
}

forest_model(model_list=fit_list, panels=panels, merge_models=T)
ggsave(filename="Images/Metabric_forest_MSH2.tiff", dpi=600)
