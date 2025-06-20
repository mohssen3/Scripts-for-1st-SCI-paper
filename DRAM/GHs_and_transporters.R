library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd("~/DRAM") #set to this directory
file_path <- "metabolism_summary.xlsx" ##put in the full path if you get an error
sheets <- excel_sheets(file_path)

data_list <- map(sheets, ~ read_excel(file_path, sheet = .x))

final_table <- bind_rows(data_list, .id = "sheet_name")

taxonomy<- read_tsv("95_Clean_classification_dereplicated_taxonomy.tsv") %>%
  mutate(across(everything(), ~ str_replace_all(.x, "d__|p__|o__|c__|f__|g__|s__", "")))


#############
###GHs
#############

carbs_groups<- read.csv("carbs_groups.csv")
taxa_dram<- carbs_groups %>% inner_join(taxonomy) %>%
  select(MAG, species, family, change_direction,Group) %>%
  filter(change_direction!="No Change")
dram_of_taxa<-final_table %>%
  pivot_longer(-c(sheet_name,gene_id,gene_description,module,header,subheader), names_to = "MAG", values_to = "gene_freq") %>%
  inner_join(taxa_dram)

heatmap_tibble<-dram_of_taxa %>%
  separate_rows(subheader, sep = ", ") %>%
  filter(str_starts(gene_id, "GH")) %>%
  mutate(subheader = str_remove_all(subheader, " Oligo| Backbone")) %>%
  filter(!is.na(subheader)) %>%
  group_by(species, family, subheader, change_direction,Group) %>%
  summarize(sum_of_gene_freq=sum(gene_freq)) %>%
  ungroup() %>%
  mutate(subheader = str_remove_all(subheader, regex(" Cleavage*", ignore_case=TRUE))) %>%
  pivot_wider(id_cols = c(species,family,change_direction,Group), names_from = subheader, values_from = sum_of_gene_freq) %>%
  select(species, where(~ any(. != 0))) %>%
  select(-c(Fucose,Rhamnose,Arabinose))

heatmap_matrix <- heatmap_tibble[, -1:-4] %>% as.matrix()
rownames(heatmap_matrix) <- heatmap_tibble$species



ha = rowAnnotation(Family=heatmap_tibble$family,
                   col = list("Direction" = c("Increasing" ="#fb8500", "Decreasing" = "#219ebc",
                                              "Decreasing (rare)" = "#21bc8d","No Change"="ivory"),
                              "Family"= family_colors),
                   annotation_name_side = "top",
                   annotation_name_rot = 75,
                   Direction=heatmap_tibble$change_direction)


heatmap_obj<- Heatmap(
  heatmap_matrix,right_annotation = ha,cluster_columns = TRUE,cluster_rows = TRUE,
  row_split = heatmap_tibble$change_direction,
  col = colorRamp2(c(min(heatmap_matrix)/3, max(heatmap_matrix)/3), c("ivory", "#84a59d")),#modified the scales here. Real numbers are 3x (change in illustrator)
  heatmap_legend_param = list(
    title = "Number of genes"
  ),
  row_title = "MAGs", row_names_side = "left",
  show_row_names = TRUE,
  show_column_names = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_names_rot = 75,
  row_names_max_width = max_text_width(rownames(heatmap_matrix)),
  column_names_max_height = unit(10, "cm"),  
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 9),
  column_names_side = "top"
)

max(heatmap_matrix)
heatmap_obj

pdf("polysaccharides_heatmap.pdf", width = 8, height = 7)  # Specify the dimensions
draw(heatmap_obj)
dev.off()

#####################
####GH boxplots
#####################

library(ggtext)
library(ggpubr)
library(rstatix)

carbs_groups<- read.csv("carbs_groups.csv")

taxa_dram<- carbs_groups %>% inner_join(taxonomy) %>%
  select(MAG, species, family, change_direction) %>%
  filter(change_direction!="No Change") %>%
  filter(change_direction!="Decreasing (rare)")


dram_of_taxa<-final_table %>%
  pivot_longer(-c(sheet_name,gene_id,gene_description,module,header,subheader), names_to = "MAG", values_to = "gene_freq") %>%
  inner_join(taxa_dram)

heatmap_tibble<-dram_of_taxa %>%
  separate_rows(subheader, sep = ", ") %>%
  filter(str_starts(gene_id, "GH")) %>%
  mutate(subheader = str_remove_all(subheader, " Oligo| Backbone")) %>%
  filter(!is.na(subheader)) %>%
  group_by(species, family, subheader, change_direction) %>%
  summarize(sum_of_gene_freq=sum(gene_freq)) %>%
  ungroup() %>%
  group_by(subheader) %>%
  filter(!all(sum_of_gene_freq == 0)) %>%
  ungroup() %>%
  mutate(sqrts=sqrt(sum_of_gene_freq)) %>%
  mutate(subheader = str_remove_all(subheader, regex(" Cleavage*", ignore_case=TRUE))) %>%
  filter(!str_detect(subheader, "Fucose|Rhamnose|Arabinose"))



discrete_palettes <- list(
  c("#219ebc","#fb8500"),
  c("#2a9d8f","#21bc8d","#e9c46a","#e76f51"))

carbs_boxplot_obj<-ggplot(heatmap_tibble,aes(x=subheader, y=sum_of_gene_freq,alpha=1)) +
  theme_classic()+
  geom_boxplot(aes(fill=change_direction),outlier.color = NA, outlier.shape = NA,size=0.3,alpha=0.4) +
  scale_color_discrete(type = discrete_palettes)+
  scale_fill_discrete(type = discrete_palettes)+
  geom_point(aes(color=change_direction),size=1,
             position=ggplot2::position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1,seed=1),
             alpha=0.5) +
  ggpubr::stat_pvalue_manual(heatmap_tibble %>%
                               select(c(subheader,change_direction,sum_of_gene_freq)) %>%
                               group_by(subheader) %>%
                               pairwise_wilcox_test(sum_of_gene_freq ~ change_direction,
                                                    paired=FALSE,
                                                    detailed=TRUE, exact=FALSE,
                                                    alternative="two.sided",
                                                    p.adjust.method ="fdr") %>%
                               ungroup() %>%
                               add_xy_position(group="change_direction", x= "subheader",
                                               fun="max",step.increase = 0.03,
                                               dodge=0.8), hide.ns = FALSE, tip.length = 0.01)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(
    title = "Glycoside Hydrolases",
    #        subtitle = "Data grouped by treatment",
    x = "Substrate",      # Custom X-axis label
    y = "Sum of MAG's Gene Counts"  # Custom Y-axis label
    #        color = "Treatment Group"  # Legend title (if not hiding)
  )
carbs_boxplot_obj

ggsave("polysaccharides_boxplots.pdf",carbs_boxplot_obj,device = "pdf",width = 8,height = 5)


#####################
###transporters
#####################
transporters_groups<- read.csv("transporters_group.csv")

taxa_dram<- transporters_groups %>% inner_join(taxonomy) %>%
  select(MAG, species, family, change_direction,Group)


dram_of_taxa<-final_table %>%
  pivot_longer(-c(sheet_name,gene_id,gene_description,module,header,subheader), names_to = "MAG", values_to = "gene_freq") %>%
  inner_join(taxa_dram)

###These are all "sugar" modules in this sheet

modules_of_interest<- paste("Arabinosaccharide transport system","Putative multiple sugar transport system",
                            "D-Allose transport system","D-Xylose transport system",
                            "Erythritol transport system","Fructose transport system",
                            "Glucose/arabinose transport system","Glucose/mannose transport system",
                            "Inositol transport system", "Inositol-phosphate transport system",
                            "L-Arabinose transport system", "L-Arabinose/lactose transport system",
                            "Maltose/maltodextrin transport system", "PTS system",
                            "Multiple sugar transport system","Oligogalacturonide transport system",
                            "Putative simple sugar transport system", "Putative sorbitol/mannitol transport system",
                            "Putative xylitol transport system", "Raffinose/stachyose/melibiose transport system",
                            "Rhamnose transport system", "Ribose transport system",
                            "Trehalose transport system", "Trehalose/maltose transport system",
                            "Putative chitobiose transport system", "Putative fructooligosaccharide transport system",
                            sep="|")
transporters<-dram_of_taxa %>%
  filter(str_detect(module, modules_of_interest)) %>%
  select(c(module, MAG, species, family, change_direction,gene_freq,Group)) %>%
  group_by(module,MAG, species, family, change_direction,Group) %>%
  summarize(sum_gene_freq=sum(gene_freq)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(MAG,species,family,change_direction,Group), names_from = module, values_from = sum_gene_freq)

heatmap_tibble<-transporters%>%
  select(-MAG) %>%
  select(species, where(~ any(. != 0))) #%>%
#  filter(change_direction!="No Change")

heatmap_matrix <- heatmap_tibble[, -1:-4] %>% as.matrix()
rownames(heatmap_matrix) <- heatmap_tibble$species

colors_9 <- brewer.pal(9, "Set3")
family_level <- sort(unique(heatmap_tibble$family))
family_colors <- setNames(colors_9[seq_along(family_level)], family_level)

ha = rowAnnotation(Family=heatmap_tibble$family,
  col = list("Direction" = c("Increasing" ="#fb8500", "Decreasing" = "#219ebc",
                             "Decreasing (rare)" = "#21bc8d","No Change"="ivory"),
             "Family"= family_colors,
             "Group" = c("Scarce" ="#ef476f", "Rich" = "#049a8f")),
  annotation_name_side = "top",
  annotation_name_rot = 75,
  Direction=heatmap_tibble$change_direction,
  Group=heatmap_tibble$Group)


max(heatmap_matrix)

heatmap_obj<-Heatmap(
  heatmap_matrix,right_annotation = ha,cluster_columns = TRUE,cluster_rows = TRUE,
  row_split=paste(heatmap_tibble$Group,heatmap_tibble$change_direction, sep="_"),
  col = colorRamp2(c(min(heatmap_matrix)/3, max(heatmap_matrix)/3), c("ivory", "#84a59d")),#modified the scales here. Real numbers are 3x (change in illustrator)
  row_title = "MAGs", row_names_side = "left",
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "Number of genes"
    ),
  rect_gp = gpar(col = "white", lwd = 1),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_names_rot = 75,
  row_names_max_width = max_text_width(rownames(heatmap_matrix)),
  column_names_max_height = unit(14, "cm"),  
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 9),
  column_names_side = "top"
)


heatmap_obj

pdf("transporters_heatmap.pdf", width = 12, height = 10)  # Specify the dimensions
draw(heatmap_obj)
dev.off()


