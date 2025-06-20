library(ggplot2)
library(ggpubr)
library(viridis)
library(tidyverse)

setwd("~/checkM_and_perc_mapped/") #set to this directory
meta<- read.csv("Meta.csv")
percent_mapped<- read.csv("percentage_mapped.csv")
a<-inner_join(percent_mapped,meta)

meta<- read.csv("Meta.csv") %>%
  filter(Type=="Fecal") %>%
  filter(Sex_SurgeryDate!="F_Day4"&Sex_SurgeryDate!="M_Day1") %>%
  filter(!str_detect(Sample,"Repeat"))


CheckM<- read.csv("CheckM_galah95.txt",sep="\t") %>%
  rename("MAG"="Bin.Id") %>%
  select(c(MAG,Completeness,Contamination))
taxonomy<- read.csv("95_Clean_classification_dereplicated_taxonomy.tsv",sep="\t")

merged_tax_qual<- inner_join(CheckM,taxonomy) %>%
  mutate(phylum = str_remove(phylum, "p__")) %>%
  mutate(family = str_remove(family, "f__")) %>%
  mutate(order = str_remove(order, "o__")) %>%
  mutate(Phylum_Family = paste0(phylum," - ",family)) %>%
  rename(Phylum=phylum) %>%
  mutate(Contamination= case_when(
    Contamination < 5 ~ 'Less than 5%',
    Contamination >= 5 & Contamination < 10 ~ '5%-10%')) %>% 
  arrange(Phylum_Family)
  

# merged_tax_qual %>%
# #  select(c(MAG,family,Completeness,Contamination)) %>%
#   ggviolin(.,x = "Contamination",y="Completeness",colour="phylum")

# ggboxplot(merged_tax_qual,x="phylum",y="Completeness",color="Contamination")+
#   geom_point()
# 
# p1<-ggboxplot(merged_tax_qual, x = "phylum", y = "Completeness",
#           color = "Contamination",add = "jitter") +
#   scale_color_viridis()+
#   labs(title = "MAG Quality (n=263)")+    #  scale_fill_brewer(palette = "Set2") +
#   theme(plot.title = element_text(hjust = 0.5))
# p1
# ggsave("MAG Quality.png",p1,device = "png",width = 11,height = 8)
overall_median <- median(merged_tax_qual$Completeness)

# Create a boxplot with a dashed line at the overall median
checkm_p<- ggplot(merged_tax_qual, aes(y = Phylum_Family, x = Completeness,
                                       color=factor(Phylum,
                                                    levels = c("Verrucomicrobiota",
                                                               "Proteobacteria",
                                                               "Firmicutes_B",
                                                               "Firmicutes_A",
                                                               "Firmicutes",
                                                               "Bacteroidota",
                                                               "Actinobacteriota"),
                                                    label=c("Verrucomicrobiota",
                                                            "Proteobacteria",
                                                            "Firmicutes_B",
                                                            "Firmicutes_A",
                                                            "Firmicutes",
                                                            "Bacteroidota",
                                                            "Actinobacteriota"))))+#,fill=Phylum)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 0.2)+
#  scale_color_manual(values=c("#F9C22E","#124559"))+
  scale_color_manual(values= c("#80B1D3","#EA9C00","#FB6250","#349BB7","#AA6A84","#068270","#FFD56B"),
                     name="Phylum")+
#  scale_fill_manual(values=c("#F9C22E","#124559"))+
#  scale_color_brewer(palette = "Set3")+
#  scale_x_continuous(breaks = c(60,80, 100)) +
  scale_x_continuous(#expand = c(0, 60),
                     limits = c(0,101),
                     breaks = c(0,20,40,60,80, 100)) +
  geom_vline(xintercept = overall_median, linetype = "dashed", color = "red", size = 0.5) +
  theme_classic()
checkm_p
ggsave("checkm_p_2.png",checkm_p,device = "png",width = 9,height = 5)
ggsave("checkm_p_2.pdf",checkm_p,device = "pdf",width = 10,height = 6)


meta<- read.csv("Meta.csv")
mapping <- c("1" = "Lam", "2" = "T10", "3" = "T4")
# Map values using indexing
meta$Injury <- mapping[meta$Group]
mapping <- c("1" = "#2a9d8f", "2" = "#e9c46a", "3" = "#e76f51")
meta$color<-mapping[meta$Group]

percent_mapped<- read.csv("percentage_mapped.csv")
perc_mapped<-inner_join(percent_mapped,meta) %>%
  mutate(Injury=factor(Injury,levels=c("Lam", "T10","T4")))
mean_perc_mapped<-mean(perc_mapped$Percentage.Mapped)
p2<-ggboxplot(perc_mapped, x = "Injury", y = "Percentage.Mapped",
          color="Injury",palette=c("#2a9d8f", "#e9c46a","#e76f51"),
          add = "jitter",add.params = list(size=0.4),
          facet.by = "Sex") +
#  scale_color_viridis()+
  labs(title = "Percent Reads Mapped to Dereplicated MAGs (n=263)")+    #  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(#expand = c(0, 60),
                     breaks = c(0,20,40,60,80, 100),
                     limits = c(0,100)) +
#  scale_x_discrete(labels = c("Lam", "T10", "T4")) +
  ylab("Reads mapped to MAGs (%)")+
  geom_hline(yintercept = mean_perc_mapped, linetype = "dashed", color = "red", size = 0.5)+
  geom_hline(yintercept = 54.7, linetype = "dashed", color = "#349BB7", size = 0.5)

p2

ggsave("Reads_to_MAG_mapped.png",p2,device = "png",width = 5.5,height = 3.5)
ggsave("Reads_to_MAG_mapped.pdf",p2,device = "pdf",width = 5.5,height = 3.5)

#  scale_color_manual(labels = c("Label", "T10", "T4"))

