library(tidyverse)
#library(VennDiagram)
#library(UpSetR)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd("~/database_comparisons") #set to this directory

input_file <- "All_databases_representative_genomes.tsv"
database_name<- "All_databases"

database_clustering_output <- read.csv(input_file,
                                       header = FALSE,
                                       stringsAsFactors = FALSE,
                                       sep= "\t")
colnames(database_clustering_output) <- c("Representative", "Member")

first_pattern<- paste0("./",database_name,"/")
second_pattern<- ".fa"

clusters <- database_clustering_output %>%
  mutate(Representative = str_replace_all(Representative,
                                          setNames(c("", ""), c(first_pattern, second_pattern))
  ),
  Member = str_replace_all(Member,
                           setNames(c("", ""), c(first_pattern, second_pattern))
  )
  ) %>%
  group_by(Representative) %>%
  summarise(Members = paste(Member, collapse = ","), .groups = "drop") %>%
  separate_longer_delim(Members, delim = ",") %>%
  group_by(Representative) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = row, values_from = Members,names_prefix = "col") %>%
  ungroup()# %>%
#  filter(str_detect(Representative, "Mouse|maxbin"))


genome_origins <- clusters %>%
  pivot_longer(cols = -Representative, names_to = "Type", values_to = "Genome") %>%
  filter(Genome != "") %>%  # Remove empty entries
  mutate(
    Database = case_when(
      grepl("CBAJ", Genome) ~ "CBAJDB",
      grepl("MGG", Genome) ~ "CMMG",
      grepl("iMGMC", Genome) ~ "iMGMC",
      grepl("MGBC", Genome) ~ "MGBC",
      grepl("Mouse|maxbin", Genome) ~ "Mohssen_etal"
    )
  ) %>%
  select(Representative, Genome, Database)# %>%

# Create presence/absence matrix
presence_absence <- genome_origins %>%
  group_by(Representative, Database) %>%
  summarise(Present = 1, .groups = "drop") %>%
  pivot_wider(names_from = Database, values_from = Present, values_fill = list(Present = 0)) %>%
  mutate(common_Mohssen_winner= ifelse(grepl("Mouse|maxbin", Representative), 1,0))


presence_absence_tibble<-presence_absence %>%
  filter(common_Mohssen_winner ==1) %>%
  select(-common_Mohssen_winner)

presence_absence_heatmap <- presence_absence_tibble[, -1] %>% as.matrix()
rownames(presence_absence_heatmap) <- presence_absence_tibble$Representative

taxonomy<- read.csv("95_Clean_classification_dereplicated_taxonomy.tsv",
                    sep="\t") %>%
  rename(Representative=MAG) %>%
  inner_join(presence_absence_tibble,.) %>%
  mutate(across(everything(), ~ str_replace_all(.x, "p__|o__|c__|f__|g__|s__", ""))) %>%
  mutate(across(everything(), ~ ifelse(.x == "", paste0("Unclassified ", cur_column()), .x)))


colors_9 <- brewer.pal(9, "Set3")

phylum_level <- sort(unique(taxonomy$phylum))
class_level <- sort(unique(taxonomy$class))
order_level <- sort(unique(taxonomy$order))
family_level <- sort(unique(taxonomy$family))

phylum_colors <- setNames(colors_9[seq_along(phylum_level)], phylum_level)
class_colors <- setNames(colors_9[seq_along(class_level)], class_level)
order_colors <- setNames(colors_9[seq_along(order_level)], order_level)
family_colors <- setNames(colors_9[seq_along(family_level)], family_level)

ha = rowAnnotation(Phylum=taxonomy$phylum,
                   Class=taxonomy$class,
                   Order=taxonomy$order,
                   Family=taxonomy$family, 
                   annotation_name_side = "top",
                   annotation_name_rot = 75,
                   col = list(Phylum = phylum_colors,
                              Class= class_colors,
                              Order = order_colors,
                              Family= family_colors),
                   Genus= anno_text(
                     taxonomy$genus,
                     gp = gpar(fontsize = 7.5, col = "black"), 
                     just = "left" 
                   ),
                   Species = anno_text(
                     taxonomy$species,
                     gp = gpar(fontsize = 7.5, col = "black"),
                     just = "left"
                   )
                   )

heatmap_obj<- Heatmap(
  presence_absence_heatmap,right_annotation = ha,
  row_split = interaction(taxonomy$phylum, taxonomy$order, sep = " & "),
  col = colorRamp2(c(0, 1), c("ivory", "#84a59d")),
  width = unit(2, "cm"), 
  height = unit(18, "cm"),
  row_title = "MAGs", row_names_side = "left",
  show_row_names = FALSE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "MAG in Dataset", 
    at = c(0, 1), 
    labels = c("Absent", "Present")
  ),
  rect_gp = gpar(col = "white", lwd = 1),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_names_rot = 75,
  row_names_max_width = max_text_width(rownames(presence_absence_heatmap)),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 9),
  column_names_side = "top"
)

heatmap_obj

pdf("All_databases_heatmap_for_53_improved_or_unique_MAGs_newcolors.pdf", width = 6, height = 8.5)  # Specify the dimensions
draw(heatmap_obj)
dev.off()
