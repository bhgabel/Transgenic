library(tidyverse)
library(ggsignif)

meta <- read.csv("metabric_prepped.csv")

#Plots ----
##Boxplots ----
##MLH1 low vs rest
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
