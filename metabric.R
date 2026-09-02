library(tidyverse)
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
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Rest", "TRUE"="Low (<12%)")) +
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
  theme_classic() + scale_x_discrete(labels=c("FALSE"="Rest", "TRUE"="Low (<8%)")) +
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
cutoff <- quantile(meta$MLH1, probs=0.12)
ggplot(data=subset(meta, PAM50 != "NA"), aes(x=PAM50, y=MLH1, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  labs(title="Metabric - MLH1 Expression by Subtype",
       x=NULL, y="mRNA Expression (Illumina microarray)") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))

cutoff <- quantile(meta$MSH2, probs=0.08)
ggplot(data=subset(meta, PAM50 != "NA"), aes(x=PAM50, y=MSH2, fill=PAM50)) +
  geom_violin() + stat_summary(fun=mean, geom="crossbar") +
  geom_hline(yintercept=cutoff, linetype="dashed", color="red") +
  scale_fill_manual(values=c('#E3E418', '#27AD81', '#31688E', '#443A83','#471164')) +
  theme_bw() +
  labs(title="Metabric - MSH2 Expression by Subtype",
       x=NULL, y="mRNA Expression (Illumina microarray)") +
  theme(legend.position="none",
        axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(1.5), hjust=0.3),
        legend.text=element_text(size=rel(1.2)),
        legend.title=element_text(size=rel(1.3)))
ggsave(filename="Images/Metabric_violin_MSH2.svg", dpi=600)

ggplot(data=subset(meta, PAM50 != "NA"), aes(x=MMR, fill=PAM50)) +
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

ggsave(filename="Images/Metabric_enrichment.tiff", dpi=600)

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


# surv.df <- meta
# surv.df$time <- ifelse(surv.df$Relapse.Months >=10, 10, surv.df$Relapse.Months)
# surv.df$event <- ifelse(surv.df$Relapse.Months >= 10, 0, surv.df$Relapse.Status.l)

# surv.df <- read.csv("Relapse_Free_METABRIC_combined_Sv_070626.csv", header=T)
# surv.df$time <- surv.df$time/12

#Disease specific survival
surv.df <- read.csv("Clinical_Overall_Survival_Data_from_METABRIC.csv", header=T)
surv.df <- left_join(surv.df, meta, by=c("sample_id"="Sample.ID"))
surv.df$event <- surv.df$dss
surv.df$time <- ifelse(surv.df$dss.time >=10, 10, surv.df$dss.time)
surv.df$event <- ifelse(surv.df$dss.time >= 10, 0, surv.df$dss)


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

fit <- survfit(Surv(time, event) ~ MMR,
               data=surv.df)
ggsurvplot(fit, data=surv.df, pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Overall)")

#MLH1

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumA"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Luminal A)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Luminal B)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Luminal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Her2"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Her2"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (HER2)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MLH1_low,
               data=subset(surv.df, PAM50 == "Basal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Basal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MLH1"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Basal)", palette=c("black", "#31688E"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ combo.mlh1,
               data=surv.df)
ggsurvplot(fit, data=surv.df,
           pval=F, xlim=c(0,10), break.time.by=1, xlab="Time (Years)",
           legend.labs=c("Luminal Rest \n(n=1157)", "Luminal MLH1 \n(n=76)",
                         "Basal Rest \n(n=212)", "Basal MLH1 \n(n=84)"),
           linetype=c(1,1,3,3), censor.shape=124, censor.size=3,
           legend.title="", surv.scale="percent",
           subtitle="(METABRIC)",
           title="Probability of Survival with MLH1 Loss",
           palette=c('black', '#31688E', 'black', '#31688E'))$plot +
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold", family="Arial"),
        plot.subtitle=element_text(hjust=0.5, size=14)) +
  annotate("text", label="Luminal HR = 1.36, p = 0.17\nBasal HR: 1.32, p = 0.19",
           x=0.1, y=0.15, hjust=0)


ggsave(filename="Images/Metabric_survival_MLH1_Combo_v2.tiff", dpi=600)


#MSH2
fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumA"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Luminal A)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Luminal B)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "LumA" | PAM50 == "LumB"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Luminal)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "Her2"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Her2"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (HER2)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == "Basal"))
ggsurvplot(fit, data=subset(surv.df, PAM50 == "Basal"), pval=T,
           xlim=c(0,10), xlab="Time (Years)", break.time.by=1,
           legend.labs=c("Rest", "Low MSH2"), font.legend=c(13),
           legend.title="", surv.scale="percent", ylab="Survival probability",
           title="Metabric (Basal)", palette=c("black", "#E3E418"))$plot +
  theme(plot.title = element_text(hjust=0.5, size=16))

fit <- survfit(Surv(time, event) ~ combo.msh2,
               data=surv.df)
ggsurvplot(fit, data=surv.df,
           pval=F, xlim=c(0,10), break.time.by=1, xlab="Time (Years)",
           legend.labs=c("Luminal Rest\n(n=1074)", "Luminal MSH2\n(n=83)",
                         "Basal Rest\n(n=192)", "Basal MSH2\n(n=20)"),
           linetype=c(1,1,3,3), censor.shape=124, censor.size=3,
           legend.title="", surv.scale="percent",
           subtitle="(METABRIC)",
           title="Probability of Survival with MSH2 Loss",
           palette=c('black', '#CC7722', 'black', '#CC7722'))$plot +
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold", family="Arial"),
        plot.subtitle=element_text(hjust=0.5, size=14)) + 
  annotate("text", label="Luminal HR = 1.17, p = 0.52\nBasal HR = 0.80, p = 0.63",
           x=0.1, y=0.15, hjust=0)
ggsave(filename="Images/Metabric_survival_MSH2_Combo_v2.svg", dpi=600)


fit <- survfit(Surv(time, event) ~ combo,
               data=surv.df)
ggsurvplot(fit, data=surv.df,
           pval=F, break.time.by=1, xlab="Time (Years)",
           legend.labs=c("Low MSH2\n(n=83)", "Rest\n(n=1150)",
                         "Low MLH1\n(n=84)", "Rest\n(n=212)"),
           linetype=c(2,1,2,1), censor.shape=124, censor=F,
           legend.title="Luminal                                            Basal",
           surv.scale="percent",
           subtitle="(METABRIC)",
           title="Probability of Survival with MMR Loss",
           palette=c('#c29f1f', '#c29f1f', "#56B4E9", "#56B4E9"))$plot + 
    theme(plot.title=element_text(hjust=0.5, size=16, face="bold", family="Arial"),
        plot.subtitle=element_text(hjust=0.5, size=14, face="bold"),
        legend.title=element_text(hjust=0.5, size=12),
        legend.key.size=unit(2, "line"), legend.location="plot",
        legend.position="bottom", legend.title.position="top") +
  annotate("text", label="Luminal HR = 1.17, p = 0.52\nBasal HR = 1.32, p = 0.19",
           x=0.1, y=0.15, hjust=0)

ggsave(filename="Images/Metabric_survival_Combo.svg", dpi=600)


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
ggsave(filename="Images/Metabric_forest_MLH1_v2.tiff", dpi=600)

#MSH2
#zero events in MSH2 subset for LumA and Basal
fit_list <- vector("list", 4)
i <- 1
for(hr in c("LumA", "LumB", "Her2", "Basal")){
  fit <- coxph(Surv(time, event) ~ MSH2_low,
               data=subset(surv.df, PAM50 == hr))
  print(hr)
  print(fit)
  fit_list[[i]] <- fit
  names(fit_list)[[i]] <- hr
  i <- i + 1
}

forest_model(model_list=fit_list, panels=panels, merge_models=T)
ggsave(filename="Images/Metabric_forest_MSH2_v2.svg", dpi=600)


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
  geom_text(aes(x = 6, label = lab), hjust = 0, size = 4) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.4, 5), clip = "off") +
  labs(x = "Hazard ratio (MLH1 Low vs Rest, log scale)", y = NULL,
       title = "Metabric Forest Plot") +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(6, 130, 6, 6),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=12))

ggsave(filename="Images/Metabric_forest_MLH1_v3.tiff", dpi=600)


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
                width = 0.22, orientation = "y", color = "black") +
  geom_point(size = 2.9, color = "black") +
  geom_text(aes(x = 5, label = lab), hjust = 0, size = 4) +
  scale_x_log10(breaks = c(0, 1, 2, 4, 6, 8, 10, 12)) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off") +
  labs(x = "Hazard ratio (MSH2 Low vs Rest, log scale)", y = NULL,
       title = "Metabric Forest Plot") +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(6, 130, 6, 6),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=12))

ggsave(filename="Images/Metabric_forest_MSH2_v3.tiff", dpi=600)
