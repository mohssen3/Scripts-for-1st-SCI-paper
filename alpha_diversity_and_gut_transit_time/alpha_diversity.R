library(tidyverse)
library(vegan)
library(ggpubr)
library(cowplot)

setwd("~/alpha_diversity_and_gut_transit_time") #change directory to this directory

trimmed_mean_tibble <- read_csv("Normalized_Trimmed_Mean.csv")


metadata <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  #  filter(Type=="Cecum"& Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4")#
  filter(Type=="Fecal") #& Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4")
abundance <- inner_join(metadata,trimmed_mean_tibble,by="Sample")[-c(2:13)] %>%
  #  column_to_rownames("Sample") %>%
  pivot_longer(-Sample, names_to = "mag", values_to="tm") %>%
  dplyr::group_by(Sample) %>%
  mutate(relabund=tm/sum(tm)) %>%
  ungroup() %>%
  select(-tm) %>%
  #  decostand(., "hellinger") %>%
  #  rownames_to_column("Sample") %>%
  as_tibble()


# 0dpi
m_metadata_0dpi_filtered <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex=="M" & Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Timepoint=="0dpi")
f_metadata_0dpi_filtered <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex=="F" & Type=="Fecal" & Sex_SurgeryDate!="F_Day4" & Timepoint=="0dpi")


m_metadata_after0_filtered <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex=="M" & Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Timepoint!="0dpi")
f_metadata_after0_filtered <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex=="F" & Type=="Fecal" & Sex_SurgeryDate!="F_Day4" & Timepoint!="0dpi")

m_abundances_0dpi_filtered <- inner_join(m_metadata_0dpi_filtered,abundance,by="Sample")[-c(2:13)]
f_abundances_0dpi_filtered <- inner_join(f_metadata_0dpi_filtered,abundance,by="Sample")[-c(2:13)]

m_abundances_after0_filtered <- inner_join(m_metadata_after0_filtered,abundance,by="Sample")[-c(2:13)]
f_abundances_after0_filtered <- inner_join(f_metadata_after0_filtered,abundance,by="Sample")[-c(2:13)]



comparisons = list(c("Lam","T10"), c("T10","T4"), c("Lam","T4"))
discrete_palettes <- list(
  c("#ef476f", "#049a8f"),
  c("#2a9d8f","#e9c46a","#e76f51"),
  RColorBrewer::brewer.pal(4, "Set2"),
  RColorBrewer::brewer.pal(6, "Accent")
)
box_plot <- function(OTU_Table,Title){
  comparisons = list(c("Lam","T10"), c("T10","T4"), c("Lam","T4"))
  p1<- OTU_Table%>%
    dplyr::group_by(Sample) %>%
    summarize(richness=specnumber(relabund),
              shannon=diversity(relabund,index="shannon"),
              simpson=diversity(relabund,index="simpson"),
              invsimpson=1/simpson) %>%
    ungroup() %>%
    pivot_longer(-Sample,names_to = "metric",values_to = "value") %>%
    inner_join(.,metadata) %>%
    #   ggplot(aes(x=Group,y=value,color=Group)) +
    #   geom_point() +
    # #  geom_smooth() +
    #   facet_wrap(~metric,nrow=4,scales="free_y") +
    filter(metric=="shannon") %>%
    mutate(Timepoint=factor(Timepoint,
                            levels=c("0dpi","7dpi","21dpi","35dpi","63dpi","6_months"),
                            labels=c("0dpi", "7dpi","21dpi","35dpi","63dpi", "180dpi")),
           Group=factor(Group,
                        levels=c("Lam","T10","T4")))%>%
    ggboxplot(., x = "Group", y = "value",ylab = "Shannon Index",
              color = "Group", 
              palette = c("#2a9d8f","#e9c46a","#e76f51"),add = "jitter",title = Title) + 
    stat_compare_means(comparisons = comparisons,method = "wilcox.test",
                       paired=FALSE,
                       method.args = list(exact=FALSE),
                       label = "p.signif",hide.ns =TRUE,
                       ref.group = ".all") +
    facet_wrap(~Sex,nrow=1, scales = "free_x")
  p1
  return(p1)

}
#without facet_wrap
box_plot(m_abundances_0dpi_filtered,"Alpha Diversity - Excluding Outliers - 0dpi - Male Mice")
box_plot(m_abundances_after0_filtered,"Alpha Diversity - Excluding Outliers - After Injury - Male Mice")

box_plot(f_abundances_0dpi_filtered,"Alpha Diversity - Excluding Outliers - 0dpi - Female Mice")
box_plot(f_abundances_after0_filtered,"Alpha Diversity - Excluding Outliers - After Injury - Female Mice")

#with facet_wrap
metadata_after0<-read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Timepoint!="0dpi")
abundance_after0<- inner_join(metadata_after0,abundance,by="Sample")[-c(2:13)]

box_plot(abundance_after0,"Alpha Diversity for Samples Collected - Male and female_relabund_beforeoutliers")
#box_plot(abundance,"Alpha Diversity for Samples Collected - Male and female_relabund_beforeoutliers_cecal")
#box_plot(abundance,"Alpha Diversity for Samples Collected - Male and female_relabund_afteroutliers_cecal")

metadata_after0_filtered<-read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex_SurgeryDate!= "M_Day1" & Sex_SurgeryDate!= "F_Day4" ) %>%
  filter(Type=="Fecal" & Timepoint!="0dpi")
abundance_after0_filtered<- inner_join(metadata_after0_filtered,abundance,by="Sample")[-c(2:13)]

box_plot(abundance_after0_filtered,"Alpha Diversity for Samples Collected - Male and female_relabund_after_removing_outliers")


metadata_0dpi_filtered<-read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Timepoint=="0dpi") %>%
  filter(Sex_SurgeryDate!= "M_Day1" & Sex_SurgeryDate!= "F_Day4" )
abundance_0dpi_filtered<- inner_join(metadata_0dpi_filtered,abundance,by="Sample")[-c(2:13)]

p<-box_plot(abundance_0dpi_filtered,"Alpha Diversity - 0dpi - Male and female_relabund_after_filtering_outliers")
p
ggsave(paste0("Alpha Diversity - 0dpi - Male and female_relabund_after_filtering_outliers",
              ".pdf"),
       p,
       device="pdf", 
       width = 5,
       height = 5)

box_plot_time <- function(OTU_Table,Title){
  comparisons = list( c("0dpi", "7dpi"), c("0dpi","21dpi"), c("0dpi","35dpi"), c("0dpi","63dpi"), c("0dpi","6_months"),
                      c("7dpi", "21dpi"), c("7dpi","35dpi"),c("7dpi","63dpi"), c("7dpi","6_months"),
                      c("21dpi", "35dpi"), c("21dpi","63dpi"), c("21dpi","6_months"),
                      c("35dpi","63dpi"), c("35dpi","6_months")
  )
  p1<- OTU_Table %>%
    # pivot_longer(-Sample, names_to = "mag", values_to="tm") %>%
    # dplyr::group_by(Sample) %>%
    # mutate(relabund=tm/sum(tm)) %>%
    # ungroup() %>%
    # select(-tm) %>%
    dplyr::group_by(Sample) %>%
    summarize(richness=specnumber(relabund),
              shannon=diversity(relabund,index="shannon"),
              simpson=diversity(relabund,index="simpson"),
              invsimpson=1/simpson) %>%
    ungroup() %>%
    pivot_longer(-Sample,names_to = "metric",values_to = "value") %>%
    inner_join(.,metadata) %>%
    #   ggplot(aes(x=Group,y=value,color=Group)) +
    #   geom_point() +
    # #  geom_smooth() +
    #   facet_wrap(~metric,nrow=4,scales="free_y") +
    filter(metric=="shannon") %>%
    mutate(Timepoint=factor(Timepoint,
                            levels=c("0dpi","7dpi","21dpi","35dpi","63dpi","180dpi"),
                            labels=c("Day0", "1 week","3 weeks","5 weeks","7 weeks", "6 months")),
           Group=factor(Group,
                        levels=c("Lam","T10","T4"))) 
  summary_data<-  p1 %>%
    select(c(Timepoint,value,Group)) %>%
    group_by(Timepoint, Group) %>%
    summarize(median_y = median(value)) %>%
    ungroup()
  plot1<- ggplot(p1,aes(x=Timepoint, y=value)) +
    geom_point(aes(color=Group),size=1,
               position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1,seed=1),
               alpha=0.5)+
    geom_boxplot(aes(fill=Group),width=0.75,alpha=0.5,outliers = FALSE,
                 linewidth=0.2)+
    scale_color_manual(name=NULL,
                       breaks=c("Lam",
                                "T10",
                                "T4"),
                       values=c("#2a9d8f","#e9c46a","#e76f51")) +
    scale_fill_manual(name=NULL,
                      breaks=c("Lam",
                               "T10",
                               "T4"),
                      values=c("#2a9d8f","#e9c46a","#e76f51"))+
    theme_classic()+
    
    # ggboxplot(., x = "Timepoint", y = "value",title = Title,ylab = "Shannon Index",
    #           color = "Group", palette = c("#2a9d8f","#e9c46a","#e76f51"),add = "jitter") + 
    ggpubr::stat_pvalue_manual(p1 %>%
                                 group_by(Timepoint) %>%
                                 wilcox_test(data=.,value ~ Group,
                                             paired=FALSE,
                                             ref.group = "Lam",
                                             detailed=TRUE,
                                             alternative="less",
                                             p.adjust.method ="fdr") %>%
                                 ungroup() %>%
                                 add_xy_position(group="Group", x= "Timepoint",
                                                 fun="median_iqr",dodge=0.75)) +
    geom_point(data= summary_data,aes(x = Timepoint, y = median_y,fill=Group),
               position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
               size = 4,alpha=1,colour="black",pch=21)+
    geom_line(data=summary_data,aes(group=Group,x = Timepoint, y = median_y,color=Group),
              position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0,seed=1),
              size=1,alpha=0.9)+
    

    labs(title="Title",
         x="Timepoint",
         y="Shannon Index")+
    theme(
      legend.background = element_rect(color="black", fill = NA),
      legend.margin = margin(l=5,t=-5, r=5, b=5),
      legend.position=c(0.9,0.9),
      #for male: legend.position=c(0.8,0.2)
      #for female: legend.position = c(0.9, 0.9),
      plot.title = element_text(hjust = 0, vjust = 0,face = "bold"),
      axis.title.x = element_text(face = "bold",margin = margin(r = 3)),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_markdown()
    )
  plot1

}

m_metadata_alltimepoints_filtered<-read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex_SurgeryDate!= "M_Day1" & Sex_SurgeryDate!= "F_Day4" ) %>%
  filter(Type=="Fecal" & Sex=="M")
m_abundance_alltimepoints_filtered<- inner_join(m_metadata_alltimepoints_filtered,abundance,by="Sample")[-c(2:13)]

male<-box_plot_time(m_abundance_alltimepoints_filtered,"Alpha Diversity for Samples - Male Mice")
male

f_metadata_alltimepoints_filtered<-read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex_SurgeryDate!= "M_Day1" & Sex_SurgeryDate!= "F_Day4" ) %>%
  filter(Type=="Fecal" & Sex=="F")
f_abundance_alltimepoints_filtered<- inner_join(f_metadata_alltimepoints_filtered,abundance,by="Sample")[-c(2:13)]

female<- box_plot_time(f_abundance_alltimepoints_filtered,"Alpha Diversity for Samples - Female Mice")
female

#### put the 2 together

cow <- plot_grid(male,female, ncol=1, align = "v", axis="1")
cow #view the multi-panel figure 
#ggsave("cowplot_alpha_diversity_v1.pdf",cow,device = "pdf",width = 7,height = 8)
ggsave("cowplot_alpha_diversity_v2.pdf",cow,device = "pdf",width = 6,height = 6)
