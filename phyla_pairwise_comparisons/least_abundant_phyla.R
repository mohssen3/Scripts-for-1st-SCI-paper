library(tidyverse)
#library(readxl)
library(broom)
library(vegan)
library(ggplot2)
library(ggtext)
library(stringr)
library(svglite)
library(ggtext)
library(ggpubr)
library(rtatix)
library(cowplot)

setwd("~/phyla_pairwise_comparisons") #set to this directory

taxonomy<- read.csv("95_Clean_classification_dereplicated_taxonomy.tsv",sep="\t") %>%
  mutate(domain = str_remove(domain, "d__")) %>%
  mutate(phylum = str_remove(phylum, "p__")) %>%
  mutate(class = str_remove(class, "c__")) %>%
  mutate(order = str_remove(order, "o__")) %>%
  mutate(family = str_remove(family, "f__")) %>%
  mutate(genus = str_remove(genus, "g__")) %>%
  mutate(species = str_remove(species, "s__")) %>%
  #  select(c(MAG,phylum,family,species))  %>%
  rename(mag=MAG)



metadata <- read_csv("metadata_w_surgerydate.csv")

trimmed_mean_tibble <- read_csv("Normalized_Trimmed_Mean.csv") %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  pivot_longer(-Sample, names_to = "mag", values_to = "tm") %>%
  group_by(Sample) %>%
  mutate(rel_abund= 100* tm/sum(tm)) %>%
  ungroup() %>%
  select(-tm) %>%
  pivot_wider(names_from = "mag",values_from = "rel_abund")




metadata <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  mutate(Group=factor(Group,levels=c("Lam", "T10","T4")),
         Timepoint=factor(Timepoint,
                          levels=c("0dpi", "7dpi","21dpi","35dpi","63dpi","6_months"),
                          labels=c("0dpi","7dpi","21dpi", "35dpi","63dpi","180dpi")))# %>%
set.seed(1995310)

per_timepoint<- function(sex,timepoint,taxonomy_level){
  metadata_time_sex<- metadata %>%
    filter(!!timepoint & Sex==sex)
  abundance_time_sex <- inner_join(metadata_time_sex,trimmed_mean_tibble,by="Sample")[-c(2:13)]
  
  joined<-abundance_time_sex %>%
    pivot_longer(-Sample, names_to = "mag", values_to = "rel_abund") %>%
    inner_join(.,metadata) %>%
    inner_join(.,taxonomy) %>%
    pivot_longer(c(domain,phylum,class,order,family,genus,species,fam_species),names_to = "level",values_to = "taxon")
  
  taxon_relabund<- joined %>% 
    filter(level==taxonomy_level) %>% 
    select(Sample,Group, mag, rel_abund,taxon) %>%
    group_by(Sample,Group,taxon) %>%
    summarize(rel_abund=sum(rel_abund)) %>%
    ungroup()# %>%
  kruskal<-taxon_relabund %>%
    nest(data = -taxon) %>%
    mutate(experimental_tests=map(.x=data,~kruskal.test(rel_abund~Group,data=.x) %>% tidy())) %>%
    unnest(experimental_tests) %>%
    mutate(p.experiment=p.adjust(p.value,method="BH")) %>%
    select(taxon,data,p.experiment) %>%
    filter(p.experiment<0.05)
  #  print(n=200) %>%
  pairwise<- kruskal %>%
    mutate(pairwise_tests = map(.x=data,
                                ~pairwise.wilcox.test(x=.x$rel_abund,
                                                      g=.x$Group, paired = FALSE,
                                                      p.adjust.method = "BH",exact=FALSE) %>%
                                  tidy())) %>%
    unnest(pairwise_tests) %>%
    filter(p.value < 0.05) %>%
    select(taxon, group1, group2, p.value) %>%
    ungroup()
  
  return(pairwise)
}

s0<- per_timepoint("M",quote(Timepoint=="0dpi"|Timepoint=="7dpi"),"phylum") %>%
  filter(group2=="Lam")
unique(s0$taxon)

s1<- per_timepoint("M",quote(Timepoint=="7dpi"|Timepoint=="21dpi"),"phylum") %>%
  filter(group2=="Lam")
unique(s1$taxon)

s2<- per_timepoint("M",quote(Timepoint=="21dpi"|Timepoint=="35dpi"),"phylum") %>%
  filter(group2=="Lam")
unique(s2$taxon)

s3<- per_timepoint("M",quote(Timepoint=="35dpi"|Timepoint=="63dpi"),"phylum") %>%
  filter(group2=="Lam")
unique(s3$taxon)

s4<- per_timepoint("M",quote(Timepoint=="63dpi"|Timepoint=="180dpi"),"phylum") %>%
  filter(group2=="Lam")
unique(s4$taxon)

all_sig_taxa<-unique(bind_rows(s1,s2,s3,s4)$taxon) #for phylum

to_plot<- function(sex,taxonomy_level,Title){
  metadata_time_sex<- metadata %>%
    filter(Sex==sex)
  abundance_time_sex <- inner_join(metadata_time_sex,trimmed_mean_tibble,by="Sample")[-c(2:13)]
  
  joined<-abundance_time_sex %>%
    pivot_longer(-Sample, names_to = "mag", values_to = "rel_abund") %>%
    inner_join(.,metadata) %>%
    inner_join(.,taxonomy) %>%
    pivot_longer(c(domain,phylum,class,order,family,genus,species,fam_species),names_to = "level",values_to = "taxon")
  
  taxon_relabund<- joined %>% 
    filter(level==taxonomy_level) %>% 
    select(Sample,Group, mag, rel_abund,taxon) %>%
    group_by(Sample,Group,taxon) %>%
    summarize(rel_abund=sum(rel_abund)) %>%
    ungroup()

  short_metadata<- metadata %>% filter(Sex==sex) %>% select(Sample,Timepoint)
  to_be_plotted<- taxon_relabund %>%
    filter(taxon_relabund$taxon %in% all_sig_taxa) %>%##
    group_by(Group,taxon) %>%
    summarize(median=median(rel_abund)) %>%
    ungroup() %>%
    group_by(taxon) %>%
    #for family and phylum
    summarize(pool=max(median) <0.25,
              median=max(median)) %>%

    ungroup() %>%
    filter(pool==TRUE) %>%
    inner_join(.,taxon_relabund, by="taxon",multiple="all") %>%
    inner_join(short_metadata) %>%
        filter(Timepoint!="0dpi") #change here ! to = as directed below
  to_be_plotted$Group <- as.character(to_be_plotted$Group)
  to_be_plotted$taxon = with(to_be_plotted, reorder(taxon,median))
  
  
  p1<- 
    ggplot(to_be_plotted,aes(x=taxon, y=rel_abund)) +
    scale_fill_manual(name=NULL,
                      breaks=c("Lam",
                               "T10",
                               "T4"),
                      values=c("#2a9d8f","#e9c46a","#e76f51")) +
    theme_classic()+
    geom_boxplot(aes(fill=Group),width=0.7,alpha=0.5,outliers = FALSE)+
    ggpubr::stat_pvalue_manual(to_be_plotted %>%
                                 group_by(taxon) %>%
                                 wilcox_test(data=.,rel_abund ~ Group,
                                             paired=FALSE,
                                             ref.group = "Lam",
                                             detailed=TRUE,
                                             p.adjust.method ="fdr") %>%
                                 ungroup() %>%
                                 rstatix::add_xy_position(group="Group",
                                                          x= "taxon",
                                                          fun="max",
                                                          dodge=0.9))+
    geom_jitter(aes(color=Group),size=0.4,position=position_jitterdodge(jitter.width = 0.03,
                                                                        dodge.width = 0.7)) +
    scale_color_manual(name=NULL,
                       breaks=c("Lam",
                                "T10",
                                "T4"),
                       values=c("black","black","black")) +
    
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1,2,3))+
    labs(title=Title,
         x="Phylum",
         y="Relative Abundance (%)")#+
  p1
  return(p1)
}
###
p1<-to_plot("M","phylum","Male Post-surgery")
p1

p3<-to_plot("F","phylum","Female Post-surgery")
p3

####change ! to = in filter(Timepoint=="0dpi")

p2<-to_plot("M","phylum","Male Pre-surgery") #change ! to =
p2

p4<-to_plot("F","phylum","Female Pre-surgery") #change ! to =
p4

##put plots together
cow1<-plot_grid(p2,p1, ncol=2, align = "v", axis="1")
cow1 #view the multi-panel figure 
ggsave("cowplot_phyla_male_leastabundant.pdf",cow1,device = "pdf",width = 6,height = 3)



cow2<-plot_grid(p4,p3, ncol=2, align = "v", axis="1")
cow2 #view the multi-panel figure 
ggsave("cowplot_phyla_female_leastabundant.pdf",cow2,device = "pdf",width = 6,height = 3)

