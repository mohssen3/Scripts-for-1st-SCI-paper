library(tidyverse)
library(ggplot2)
library(cowplot)

setwd("~/BMS_plot") #set to this directory
metadata <- read_csv("./metadata_w_surgerydate.csv")

p1<- metadata %>%
  filter(Type=="Fecal") %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  select(Mouse,`6-month BMS score`,Group,Sex) %>%
  distinct()
title=paste0("Before Removing Outliers (n = ",
             length(p1$Mouse),
             " Mice)")
discrete_palettes <- list(
  c("#ef476f", "#049a8f"),
  c("#2a9d8f","#e9c46a","#e76f51"),
  RColorBrewer::brewer.pal(4, "Set2"),
  RColorBrewer::brewer.pal(6, "Accent")
)
p2_before<- p1%>%
  ggplot(.,aes(x=Group,y=`6-month BMS score`,color=Sex)) +
  theme(plot.title = element_text(hjust = 0.5))+
#  geom_violin() +
  geom_point(position=position_dodge2(width = .5)) +
  labs(title=str_wrap(title,60)) +
  scale_color_discrete(type=discrete_palettes)+
  ylim(0,9)+
  scale_y_continuous(breaks=seq(0,9,by=1)) +
  theme_classic()

p2_before
#ggsave(paste0("Figure1A2",".png"),p2_before,device="png",width=5,height = 3)
#ggsave(paste0("Figure1A2",".pdf"),p2_before,device="pdf",width=5,height = 3)

p1<- metadata %>%
  filter(Type=="Fecal"&
           Sex_SurgeryDate!="M_Day1"&
           Sex_SurgeryDate!="F_Day4") %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  select(Mouse,`6-month BMS score`,Group,Sex) %>%
  distinct() 
title=paste0("After Removing Outliers (n = ",
             length(p1$Mouse),
             " Mice)")
p2_after<- p1 %>%
  ggplot(.,aes(x=Group,y=`6-month BMS score`,color=Sex)) +
  #  geom_violin() +
  geom_point(position=position_dodge2(width = .5)) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title=str_wrap(title,60))+
  theme_classic()+
  ylim(0,9)+
  scale_y_continuous(breaks=seq(0,9,by=1))
  

p2_after
#ggsave(paste0("Figure1A2_filtered",".png"),p2_after,device="png",width=5,height = 3)
#ggsave(paste0("Figure1A2_filtered",".pdf"),p2_after,device="pdf",width=5,height = 3)



cow <- plot_grid(p2_before,p2_after, ncol=2, align = "v", axis="1")
cow #view the multi-panel figure 
ggsave("cowplot_bms.pdf",cow,device = "pdf",width = 9,height = 3)
