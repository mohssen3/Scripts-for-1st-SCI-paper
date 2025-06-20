library(tidyverse)
library(ggpubr)
library(rtatix)
library(cowplot)

setwd("~/alpha_diversity_and_gut_transit_time") #set to this directory

gut_transit<-read_csv("GI transit 102824 Red carmine.csv",skip = 2)

comparisons = list(c("Lam","T9"))
discrete_palettes <- list(
  c("#ef476f", "#049a8f"),
  c("#2a9d8f","#e9c46a","#e76f51"),
  RColorBrewer::brewer.pal(4, "Set2"),
  RColorBrewer::brewer.pal(6, "Accent")
)
df<-gut_transit%>%
    pivot_longer(-c(Mouse,Injury),names_to = "Timepoint",values_to = "Time") %>%
    mutate(Timepoint=factor(Timepoint,
                            levels=c("7dpi","15dpi","21dpi"),
                            labels=c("1 week","2 weeks","3 weeks")),
           Injury=factor(Injury,
                        levels=c("Lam","T9"),
                        labels=c("Lam", "T9")))
p1<-df %>%
  ggboxplot(., x = "Timepoint", y = "Time",ylab = "Time (mins)",
              color = "Injury", 
              palette = c("#2a9d8f","#e9c46a"),add = "jitter",title = "Gut Transit Time") + 
  ggpubr::stat_pvalue_manual(df %>%
                               select(c(Timepoint,Injury,Time)) %>%
                               group_by(Timepoint) %>%
                               pairwise_wilcox_test(Time ~ Injury, paired=FALSE,
                                                    ref.group = "Lam",
                                                    detailed=TRUE,
                                                    alternative="two.sided",
                                                    p.adjust.method ="fdr") %>%
                               ungroup() %>%
                               add_xy_position(group="Injury", x= "Timepoint",
                                               fun="max",step.increase = 0.03,
                                               dodge=0.8), hide.ns = TRUE, tip.length = 0.01)
p1  

ggsave("gut_transit_time.pdf", p1, device="pdf", width = 3, height = 4)
