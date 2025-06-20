library(vegan)
library(ggplot2)
library(ggpubr)
library(rtatix)
library(tidyverse)
library(cowplot)


setwd("~/bray_relative_to0dpi") #set to this directory

trimmed_mean_tibble <- read_csv("Normalized_Trimmed_Mean.csv")

metadata <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" )

abundance <- inner_join(metadata,trimmed_mean_tibble,by="Sample")[-c(2:13)] %>%
  column_to_rownames("Sample")


calc_dist <- function(OTU_Table,Normalization_method) {
  if(missing(Normalization_method)) {
    bray <- vegdist(OTU_Table,method = "bray")
  } else {
    transformed<-decostand(OTU_Table,Normalization_method)
    bray <- vegdist(transformed,method = "bray")
  }
}

bray<- calc_dist(abundance)
discrete_palettes <- list(
  c("#ef476f", "#049a8f"),
  c("#2a9d8f","#e9c46a","#e76f51"))


bray_plot<-
  bray %>%
  as.matrix() %>%
  as_tibble(rownames= "Sample") %>%
  pivot_longer(-Sample) %>%
  filter(Sample < name) %>%
  inner_join(.,metadata,by=c("Sample" = "Sample")) %>%
  inner_join(.,metadata,by=c("name" = "Sample")) %>%
  mutate(day.x=as.numeric(str_replace_all(Timepoint.x, c("6_months"= "180","dpi"=""))),
         day.y=as.numeric(str_replace_all(Timepoint.y, c("6_months"= "180","dpi"="")))
  ) %>%
  filter(Mouse.x == Mouse.y & day.x==0 & day.y!=0 & Type.y=="Fecal" & Type.x=="Fecal" )# %>%#& day.y!=180) %>%

bray_plot$Timepoint.y = factor(bray_plot$Timepoint.y, 
                               levels = c("0dpi", "7dpi", "21dpi",
                                          "35dpi","63dpi","6_months"))

summary_data <- bray_plot %>%
  select(c(day.y,value,Mouse.x,Group.y,Timepoint.y)) %>%
  group_by(Timepoint.y, Group.y) %>%
  summarize(median_y = median(value), sd_y = sd(value)) %>%
  ungroup() %>%
  mutate(y_min = median_y - 2*sd_y,
         y_max = median_y + 2*sd_y)

p1<- ggplot(bray_plot,aes(x=Timepoint.y, y=value,alpha=0.01)) +
  theme_classic()+
  geom_boxplot(aes(fill=Group.y),outlier.color = NA,size=0.3,alpha=0.2) +
  scale_color_discrete(type = discrete_palettes)+
  scale_fill_discrete(type = discrete_palettes)+
  geom_line(data=summary_data,aes(group=Group.y,x = Timepoint.y, y = median_y,color=Group.y),
            position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
            size=4,alpha=0.9)+
  geom_point(aes(color=Group.y),size=1,
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1,seed=1),
             alpha=0.5)+
  geom_point(data=summary_data, aes(x = Timepoint.y, y = median_y,fill=Group.y),
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
             size = 5,alpha=1,colour="black",pch=21)+
  stat_pvalue_manual(bray_plot %>%
                       select(c(Timepoint.y,Group.y,value)) %>%
                       group_by(Timepoint.y) %>%
                       wilcox_test(value ~ Group.y, paired=FALSE, 
                                   ref.group = "Lam",
                                   detailed=TRUE,
                                   alternative="less",
                                   p.adjust.method ="fdr") %>%
                       ungroup() %>%
                       add_xy_position(group="Group.y", x= "Timepoint.y",
                                       fun="max",
                                       dodge=0.8))
p1

########### repeat for male and female after removing outliers (or write a function for it)
metadata_filtered_m <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4" & Sex=="M")

abundance_filtered_m <- inner_join(metadata_filtered_m,trimmed_mean_tibble,by="Sample")[-c(2:13)] %>%
  column_to_rownames("Sample")



calc_dist <- function(OTU_Table,Normalization_method) {
  if(missing(Normalization_method)) {
    bray <- vegdist(OTU_Table,method = "bray")
  } else {
    transformed<-decostand(OTU_Table,Normalization_method)
    bray <- vegdist(transformed,method = "bray")
  }
}


bray<- calc_dist(abundance_filtered_m)
discrete_palettes <- list(
  c("#ef476f", "#049a8f"),
  c("#2a9d8f","#e9c46a","#e76f51"))


bray_plot_m<-
  bray %>%
  as.matrix() %>%
  as_tibble(rownames= "Sample") %>%
  pivot_longer(-Sample) %>%
  filter(Sample < name) %>%
  inner_join(.,metadata_filtered_m,by=c("Sample" = "Sample")) %>%
  inner_join(.,metadata_filtered_m,by=c("name" = "Sample")) %>%
  mutate(day.x=as.numeric(str_replace_all(Timepoint.x, c("6_months"= "180","dpi"=""))),
         day.y=as.numeric(str_replace_all(Timepoint.y, c("6_months"= "180","dpi"="")))
  ) %>%
  filter(Mouse.x == Mouse.y & day.x==0 & day.y!=0 & Type.y=="Fecal" & Type.x=="Fecal" )# day.y!=180)
#write.csv(bray_plot, "bray_plot.csv")

bray_plot_m$Timepoint.y = factor(bray_plot_m$Timepoint.y, 
                               levels = c("0dpi", "7dpi", "21dpi",
                                          "35dpi","63dpi","6_months"))

summary_data <- bray_plot_m %>%
  select(c(day.y,value,Mouse.x,Group.y,Timepoint.y)) %>%
  group_by(Timepoint.y, Group.y) %>%
  summarize(median_y = median(value), sd_y = sd(value)) %>%
  ungroup() %>%
  mutate(y_min = median_y - 2*sd_y,
         y_max = median_y + 2*sd_y)

p2<- ggplot(bray_plot_m,aes(x=Timepoint.y, y=value,alpha=0.01))+
  theme_classic()+
  geom_boxplot(aes(fill=Group.y),outlier.color = NA,size=0.3,alpha=0.2) +
  scale_color_discrete(type = discrete_palettes)+
  scale_fill_discrete(type = discrete_palettes)+
  geom_line(data=summary_data,aes(group=Group.y,x = Timepoint.y, y = median_y,color=Group.y),
            position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
            size=4,alpha=0.9)+
  geom_point(aes(color=Group.y),size=1,
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1,seed=1),
             alpha=0.5)+
  geom_point(data=summary_data, aes(x = Timepoint.y, y = median_y,fill=Group.y),
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
             size = 5,alpha=1,colour="black",pch=21)+
  stat_pvalue_manual(bray_plot_m %>%
                       select(c(Timepoint.y,Group.y,value)) %>%
                       group_by(Timepoint.y) %>%
                       wilcox_test(value ~ Group.y, paired=FALSE, 
                                   ref.group = "Lam",
                                   detailed=TRUE,
                                   alternative="less",
                                   p.adjust.method ="fdr") %>%
                       # ungroup() %>%
                       add_xy_position(group="Group.y", x= "Timepoint.y",
                                       fun="max",
                                       dodge=0.8))
p2


########### Now female
metadata_filtered_f <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4" & Sex=="F")

abundance_filtered_f <- inner_join(metadata_filtered_f,trimmed_mean_tibble,by="Sample")[-c(2:13)] %>%
  column_to_rownames("Sample")



calc_dist <- function(OTU_Table,Normalization_method) {
  if(missing(Normalization_method)) {
    bray <- vegdist(OTU_Table,method = "bray")
  } else {
    transformed<-decostand(OTU_Table,Normalization_method)
    bray <- vegdist(transformed,method = "bray")
  }
}

bray<- calc_dist(abundance_filtered_f)
discrete_palettes <- list(
  c("#ef476f", "#049a8f"),
  c("#2a9d8f","#e9c46a","#e76f51"))

bray_plot_f<-
  bray %>%
  as.matrix() %>%
  as_tibble(rownames= "Sample") %>%
  pivot_longer(-Sample) %>%
  filter(Sample < name) %>%
  inner_join(.,metadata_filtered_f,by=c("Sample" = "Sample")) %>%
  inner_join(.,metadata_filtered_f,by=c("name" = "Sample")) %>%
  mutate(day.x=as.numeric(str_replace_all(Timepoint.x, c("6_months"= "180","dpi"=""))),
         day.y=as.numeric(str_replace_all(Timepoint.y, c("6_months"= "180","dpi"="")))
  ) %>%
  filter(Mouse.x == Mouse.y & day.x==0 & day.y!=0 & Type.y=="Fecal" & Type.x=="Fecal")# & day.y==180)

bray_plot_f$Timepoint.y = factor(bray_plot_f$Timepoint.y, 
                                 levels = c("0dpi", "7dpi", "21dpi",
                                            "35dpi","63dpi","6_months"))

summary_data <- bray_plot_f %>%
  select(c(day.y,value,Mouse.x,Group.y,Timepoint.y)) %>%
  group_by(Timepoint.y, Group.y) %>%
  summarize(median_y = median(value), sd_y = sd(value)) %>%
  ungroup() %>%
  mutate(y_min = median_y - 2*sd_y,
         y_max = median_y + 2*sd_y)

p3<- ggplot(bray_plot_f,aes(x=Timepoint.y, y=value,alpha=0.01))+
  theme_classic()+
  geom_boxplot(aes(fill=Group.y),outlier.color = NA,size=0.3,alpha=0.2) +
  scale_color_discrete(type = discrete_palettes)+
  scale_fill_discrete(type = discrete_palettes)+
  geom_line(data=summary_data,aes(group=Group.y,x = Timepoint.y, y = median_y,color=Group.y),
            position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
            size=4,alpha=0.9)+
  geom_point(aes(color=Group.y),size=1,
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1,seed=1),
             alpha=0.5)+
  geom_point(data=summary_data, aes(x = Timepoint.y, y = median_y,fill=Group.y),
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
             size = 5,alpha=1,colour="black",pch=21)+
  stat_pvalue_manual(bray_plot_f %>%
                       select(c(Timepoint.y,Group.y,value)) %>%
                       group_by(Timepoint.y) %>%
                       wilcox_test(value ~ Group.y, paired=FALSE, 
                                   ref.group = "Lam",
                                   detailed=TRUE,
                                   alternative="less",
                                   p.adjust.method ="fdr") %>%
                       
                       ungroup() %>%
                       add_xy_position(group="Group.y", x= "Timepoint.y",
                                       fun="max",
                                       dodge=0.8))
p3

#### put all 3 together

cow <- plot_grid(p1,p2,p3, ncol=1, align = "v", axis="1")
cow #view the multi-panel figure 
ggsave("cowplot_3brayRelativeto0.pdf",cow,device = "pdf",width = 7,height = 8)


