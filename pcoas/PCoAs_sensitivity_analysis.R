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
###Ordination
setwd("~/pcoas") #set to this directory

library(ggpubr)
##tutorial

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

pcoa_plot <- function(dist_matrix, Title,Factor) {
  #  bray <- vegdist(OTU_Table,method = "bray")
  pcoa <- cmdscale(dist_matrix, eig=TRUE,add=TRUE)
  positions= pcoa$points
  
  colnames(positions) <- c("pcoa1","pcoa2")
  percent_explained <- 100*pcoa$eig /sum(pcoa$eig)
  pretty_pe <- format(round(percent_explained[1:2], digits=1),nsmall=1,trim=TRUE)
  labs <- c(glue("PCO 1 ({pretty_pe[1]}%)"),glue("PCO 2 ({pretty_pe[2]}%)"))
  par(oma=c(0,0,2,0))
  pcoa_metadata<- positions %>%
    as_tibble(rownames = "Sample") %>%
    inner_join(.,metadata)
  pcoa_metadata$Timepoint <-factor(pcoa_metadata$Timepoint,
                                   levels=c("0dpi","7dpi","21dpi",
                                            "35dpi","63dpi", "6_months"))
  pcoa_centroid <- pcoa_metadata %>% 
    group_by(get(Factor)) %>%
    summarize(pcoa1=mean(pcoa1),pcoa2=mean(pcoa2)) %>%
    rename_with(.cols = 1, ~toString(Factor))
  p1 <- pcoa_metadata %>%
    ggplot(aes(x=pcoa1,y=pcoa2,color=get(Factor)))+
    geom_point(aes(text=paste0(Sample,"\n",Mouse,
                               "\nInjury Level: ",Group,
                               "\nTimepoint: ",Timepoint,
                               "\nSex: ",Sex,
                               "\nSequencing Depth: ", Reads,
                               " MB\nType: ",Type,
                               "\nSurgery Date: ", Surgery_date,
                               "\n6_month BMS Score: ",`6_month_BMS_score`))) +
    stat_ellipse(type='t',level=0.7,show.legend = FALSE)+
    labs(x=labs[1],y=labs[2]) +
    geom_point(data=pcoa_centroid,size=5,shape=21,color="black",
               aes(fill=get(Factor)) , show.legend = FALSE) +
    #    guides(fill=guide_legend(title="New Legend Title"))
    scale_color_discrete(name=Factor) +
    theme_classic() +
    theme(plot.caption=element_text(hjust=0))
  
  #    scale_color_discrete(name=toString(Factor),
  #                      labels=pcoa_metadata[as.name(Factor)])
  
  
  p2<- p1 +theme(plot.margin = margin(0.5, 1, 0.5, 1, "cm")) + 
    ggtitle(Title) + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(Title,".png"),p2,device="png",width=9,height = 7)
  p3<-ggplotly(p2, tooltip="text")
  htmlwidgets::saveWidget(p3, paste0(Title,".html"))
  ggplotly(p2, tooltip="text")
  #ggsave(pcoa)
}



joined <- inner_join(metadata,trimmed_mean_tibble,by="Sample")

###fecal no repeats
Fecal_no_repeats<- joined %>% filter(Type=="Fecal") %>% 
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_mat <- calc_dist(Fecal_no_repeats)
pcoa_plot(distance_Fecal_no_repeats_mat,"Fecal_No_repeats_colored_by_Sex_Date","Sex_SurgeryDate")
pcoa_plot(distance_Fecal_no_repeats_mat,"Fecal_No_repeats_colored_by_Sex","Sex")

distance_tibble<- distance_Fecal_no_repeats_mat %>% as.matrix() #%>% as_tibble()
rownames(distance_tibble) <- colnames(distance_tibble)
distance_tibble <- distance_tibble %>% as_tibble(rownames="Sample")
distance_metadata<- inner_join(distance_tibble,metadata)

all_dist<-adonis.pair(distance_Fecal_no_repeats_mat,as.factor(distance_metadata$Sex_SurgeryDate),nper=9999,corr.method="BH")
write.csv(all_dist,"Sex_and_SurgeryDate_adonis_BH_corrected.csv")

##################
## 2 mice (partial injury) excluded
Fecal_no_repeats_mice_excluded<- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Mouse != "Mouse54" & Mouse != "Mouse9") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_mice_excluded_mat <- calc_dist(Fecal_no_repeats_mice_excluded)
pcoa_plot(distance_Fecal_no_repeats_mice_excluded_mat,"Fecal_No_repeats_mice_excluded_colored_by_Sex_Date","Sex_SurgeryDate")
pcoa_plot(distance_Fecal_no_repeats_mice_excluded_mat,"Fecal_No_repeats_mice_excluded_colored_by_Sex","Sex")
#############
#####0dpi
Fecal_no_repeats_0dpi<- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_0dpi_mat <- calc_dist(Fecal_no_repeats_0dpi)
distance_tibble<- distance_Fecal_no_repeats_0dpi_mat %>% as.matrix() #%>% as_tibble()
rownames(distance_tibble) <- colnames(distance_tibble)
distance_tibble <- distance_tibble %>% as_tibble(rownames="Sample")
distance_metadata<- inner_join(distance_tibble,metadata)

pcoa_plot(distance_Fecal_no_repeats_0dpi_mat,"Fecal_No_repeats_0dpi_colored_by_Sex_Date","Sex_SurgeryDate")
pcoa_plot(distance_Fecal_no_repeats_0dpi_mat,"Fecal_No_repeats_0dpi_colored_by_Sex","Sex")


a<-adonis.pair(distance_Fecal_no_repeats_0dpi_mat,as.factor(distance_metadata$Sex_SurgeryDate),nper=9999,corr.method="BH")
write.csv(a,"Sex_and_SurgeryDate_adonis_BH_corrected_0dpi.csv")

#########
##after0dpi
Fecal_no_repeats_after0 <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_after0_mat <- calc_dist(Fecal_no_repeats_after0)

pcoa_plot(distance_Fecal_no_repeats_after0_mat,"Fecal_No_repeats_after0_colored_by_Sex_Date","Sex_SurgeryDate")
pcoa_plot(distance_Fecal_no_repeats_after0_mat,"Fecal_No_repeats_after0_colored_by_Sex","Sex")
pcoa_plot(distance_Fecal_no_repeats_after0_mat,"Fecal_No_repeats_after0_colored_by_Injury_Level","Group")

######

##after0dpi male day1excluded
Fecal_no_repeats_after0_m <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi" & Sex=="M" & Sex_SurgeryDate!="M_Day1") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_after0_m_mat <- calc_dist(Fecal_no_repeats_after0_m)

pcoa_plot(distance_Fecal_no_repeats_after0_m_mat,"Fecal_No_repeats_after0_male_excluding_day1_colored_by_Injury_Level","Group")

##after0dpi female day1excluded
Fecal_no_repeats_after0_f <- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint!="0dpi" & Sex=="F" & Sex_SurgeryDate!="F_Day4") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_after0_f_mat <- calc_dist(Fecal_no_repeats_after0_f)

pcoa_plot(distance_Fecal_no_repeats_after0_f_mat,"Fecal_No_repeats_after0_female_excluding_day4_colored_by_Injury_Level","Group")
pcoa_plot(distance_Fecal_no_repeats_after0_f_mat,"Fecal_No_repeats_after0_female_excluding_day4_colored_by_Injury","disease_status")

#########################
Fecal_no_repeats_0dpi<- joined %>% filter(Type=="Fecal" ) %>% 
  filter(Timepoint=="0dpi") %>%
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_0dpi_mat <- calc_dist(Fecal_no_repeats_0dpi)

distance_tibble<- distance_Fecal_no_repeats_0dpi_mat %>% as.matrix() #%>% as_tibble()
rownames(distance_tibble) <- colnames(distance_tibble)
distance_tibble <- distance_tibble %>% as_tibble(rownames="Sample")
distance_metadata<- inner_join(distance_tibble,metadata)

a<-adonis.pair(distance_Fecal_no_repeats_0dpi_mat,as.factor(distance_metadata$Sex_SurgeryDate),nper=9999,corr.method="bonferroni")
write.csv(a,"Sex_and_SurgeryDate_adonis_bonferroni_corrected_0dpi.csv")

##
Fecal_no_repeats<- joined %>% filter(Type=="Fecal") %>% 
  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  column_to_rownames("Sample")
distance_Fecal_no_repeats_mat <- calc_dist(Fecal_no_repeats)

distance_tibble<- distance_Fecal_no_repeats_mat %>% as.matrix() #%>% as_tibble()
rownames(distance_tibble) <- colnames(distance_tibble)
distance_tibble <- distance_tibble %>% as_tibble(rownames="Sample")
distance_metadata<- inner_join(distance_tibble,metadata)

all_dist<-adonis.pair(distance_Fecal_no_repeats_mat,as.factor(distance_metadata$Sex_SurgeryDate),nper=9999,corr.method="bonferroni")
write.csv(all_dist,"Sex_and_SurgeryDate_adonis_bonferroni_corrected_alldata.csv")
