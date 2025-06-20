library(tidyverse)
library(readxl)
library(broom)
library(vegan)
library(glue)
library(ggplot2)
library(plotly)
library(patchwork)
library(expm)
require(matrixcalc)
#remotes::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

setwd("~/pcoas/") #set to this directory

trimmed_mean_tibble <- read_csv("Normalized_Trimmed_Mean.csv")
trimmed_mean <- read_csv("Normalized_Trimmed_Mean.csv") %>%
  column_to_rownames("Sample")


metadata <- read_csv("metadata_w_surgerydate.csv")


calc_dist <- function(OTU_Table,Normalization_method) {
  if(missing(Normalization_method)) {
    bray <- vegdist(OTU_Table,method = "bray")
  } else {
    transformed<-decostand(OTU_Table,Normalization_method)
    bray <- vegdist(transformed,method = "bray")
  }
}


pcoa_plot <- function(OTU_Table, Title,Factor) {
  dist_matrix <- vegdist(OTU_Table,method = "bray")
  pcoa <- cmdscale(dist_matrix, eig=TRUE,add=TRUE)
  positions= pcoa$points
  
  colnames(positions) <- c("pcoa1","pcoa2")
  percent_explained <- 100*pcoa$eig /sum(pcoa$eig)
  pretty_pe <- format(round(percent_explained[1:2], digits=1),nsmall=1,trim=TRUE)
  labs <- c(glue("PCo 1 ({pretty_pe[1]}%)"),glue("PCo 2 ({pretty_pe[2]}%)"))
  par(oma=c(0,0,2,0))
  pcoa_metadata<- positions %>%
    as_tibble(rownames = "Sample") %>%
    inner_join(.,metadata)
  
  distance_tibble<- dist_matrix %>% as.matrix() #%>% as_tibble() 
  rownames(distance_tibble) <- colnames(distance_tibble)
  distance_tibble <- distance_tibble %>% as_tibble(rownames="Sample")
  distance_metadata<- inner_join(distance_tibble,metadata)
  
  test <-adonis2(dist_matrix~get(Factor),data=distance_metadata,permutations=9999)
  test
  
  pcoa_centroid <- pcoa_metadata %>% 
    group_by(get(Factor)) %>%
    summarize(pcoa1=mean(pcoa1),pcoa2=mean(pcoa2)) %>%
    rename_with(.cols = 1, ~toString(Factor))
  discrete_palettes <- list(
    c("#ef476f", "#049a8f"),
    c("#2a9d8f","#e9c46a","#e76f51"),
    RColorBrewer::brewer.pal(4, "Set2"),
    RColorBrewer::brewer.pal(6, "Accent")
  )
  p1 <- pcoa_metadata %>%
    ggplot(aes(x=pcoa1,y=pcoa2,color=get(Factor)))+
    geom_point(aes(text=paste0(Sample,"\n",Mouse,
                               "\nInjury Level: ",Group,
                               "\nTimepoint: ",Timepoint,
                               "\nSex: ",Sex,
                               "\nSequencing Depth: ", Reads,
                               " MB\nType: ",Type,
                               "\nSurgery Date: ", Surgery_date,
                               "\n6-month BMS Score: ",`6_month_BMS_score`))) +
    stat_ellipse(type='t',level=0.75,show.legend = FALSE)+
    labs(x=labs[1],y=labs[2],
         caption=paste("adonis (permanova) pvalue of: ",test$'Pr(>F)'[1])) +
    geom_point(data=pcoa_centroid,size=5,shape=21,color="black",
               aes(fill=get(Factor)) , show.legend = FALSE) +
    scale_color_discrete(name=Factor,type = discrete_palettes) +
    theme_classic() +
    theme(plot.caption=element_text(hjust=0))

  
  p2<- p1 +theme(plot.margin = margin(0.5, 1, 0.5, 1, "cm")) + 
    ggtitle(Title) + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(Title,".png"),p2,device="png",width=6,height = 4)
  ggsave(paste0(Title,".pdf"),p2,device="pdf",width=6,height = 4)
  p3<-ggplotly(p2, tooltip="text")
  htmlwidgets::saveWidget(p3, paste0(Title,".html"))
  ggplotly(p2, tooltip="text")
}


joined <- inner_join(metadata,trimmed_mean_tibble,by="Sample")


#############
#####0dpi
Fecal_no_repeats_0dpi<- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")

distance_Fecal_no_repeats_0dpi_mat <- calc_dist(Fecal_no_repeats_0dpi)
Factor="Sex"
pcoa_plot(distance_Fecal_no_repeats_0dpi_mat,"PCoA of Fecal Samples only at 0dpi Before Excluding Outliers","Sex")

####
#0dpi after excluding day1M,day4F
day0_excl<- joined %>% 
  filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi") %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  column_to_rownames("Sample")

distance_Fecal_no_repeats_0dpi_mat <- calc_dist(day0_excl)
Factor="Sex"
pcoa_plot(distance_Fecal_no_repeats_0dpi_mat,"PCoA of Fecal Samples only at 0dpi - Excluding Outliers","Sex")


#########
##after0dpi
Fecal_no_repeats_after0 <- joined %>% 
  filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")

distance_Fecal_no_repeats_after0_mat <- calc_dist(Fecal_no_repeats_after0)

Factor="Sex"
pcoa_plot(distance_Fecal_no_repeats_after0_mat,"PCoA of Fecal Samples only AFTER 0dpi Before Excluding Outliers","Sex")
# Factor="Group"
# pcoa_plot(distance_Fecal_no_repeats_after0_mat,"Bray_Fecal_No_repeats_after0_colored_by_Injury_Level","Group")


####
#After0, excluding mday1,fday4
Fecal_no_repeats_after0 <- joined %>% 
  filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi") %>%
  filter(Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_after0_mat <- calc_dist(Fecal_no_repeats_after0)
# Factor="Sex_SurgeryDate"
# pcoa_plot(distance_Fecal_no_repeats_after0_mat,"Bray_excl_Fecal_No_repeats_after0_colored_by_Sex_Date","Sex_SurgeryDate")
Factor="Sex"
pcoa_plot(distance_Fecal_no_repeats_after0_mat,"PCoA of Fecal Samples only AFTER 0dpi - Excluding Outliers","Sex")
# Factor="Group"
# pcoa_plot(distance_Fecal_no_repeats_after0_mat,"Bray_excl_Fecal_No_repeats_after0_colored_by_Injury_Level","Group")

##after0dpi male day1excluded
Fecal_no_repeats_after0_m <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi" & Sex=="M" & Sex_SurgeryDate!="M_Day1") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_after0_m_mat <- calc_dist(Fecal_no_repeats_after0_m)
Factor="Group"
pcoa_plot(distance_Fecal_no_repeats_after0_m_mat,"PCoA of Fecal Samples only AFTER 0dpi - Excluding Outliers\nMale Mice only","Group")

##0dpi male day1excluded
Fecal_no_repeats_at0_m <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi" & Sex=="M" & Sex_SurgeryDate!="M_Day1") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_at0_m_mat <- calc_dist(Fecal_no_repeats_at0_m)
Factor="Group"
pcoa_plot(distance_Fecal_no_repeats_at0_m_mat,"PCoA of Fecal Samples only at 0dpi - Excluding Outliers\nMale Mice only","Group")


##after0dpi female day4excluded
Fecal_no_repeats_after0_f <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi" & Sex=="F" & Sex_SurgeryDate!="F_Day4") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_after0_f_mat <- calc_dist(Fecal_no_repeats_after0_f)

Factor="Group"
pcoa_plot(distance_Fecal_no_repeats_after0_f_mat,"PCoA of Fecal Samples only AFTER 0dpi - Excluding Outliers\nFemale Mice only","Group")

##0dpi female day4excluded
Fecal_no_repeats_at0_f <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi" & Sex=="F" & Sex_SurgeryDate!="F_Day4") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_at0_f_mat <- calc_dist(Fecal_no_repeats_at0_f)

Factor="Group"
pcoa_plot(distance_Fecal_no_repeats_at0_f_mat,"PCoA of Fecal Samples only at 0dpi - Excluding Outliers\nFemale Mice only","Group")

##0dpi male day1excluded and female day4excluded
Fecal_no_repeats_at0 <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_at0_mat <- calc_dist(Fecal_no_repeats_at0)

Factor="Group"
pcoa_plot(distance_Fecal_no_repeats_at0_mat,"PCoA of Fecal Samples only at 0dpi - Excluding Outliers\n All Mice","Group")
