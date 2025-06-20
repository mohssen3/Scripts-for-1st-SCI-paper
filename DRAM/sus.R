library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd("~/DRAM") #set to this directory

taxonomy<- read_tsv("95_Clean_classification_dereplicated_taxonomy.tsv") %>%
  mutate(across(everything(), ~ str_replace_all(.x, "d__|p__|o__|c__|f__|g__|s__", "")))
DRAM<- read_tsv("merged_annotations.tsv") %>%
  rename(MAG=fasta) #%>%
#  select(MAG,rank,kegg_genes_id,ko_id,kegg_hit,uniref_id,uniref_hit,pfam_hits)
  
carbs_groups<- read.csv("carbs_groups.csv")

taxa_dram<- carbs_groups %>% inner_join(taxonomy) %>%
  select(MAG, species, family, change_direction,Group)# %>%
#  filter(change_direction!="No Change")
dram_of_taxa<-DRAM %>%
  inner_join(taxa_dram,.)

sus_systems<- c("susR","susA", "susB", "susC","susD","susE", "susF","susG")
sus_dram<-dram_of_taxa %>%
  filter(if_any(everything(), ~ str_detect(., regex(paste(sus_systems, collapse = "|"), ignore_case = TRUE))))

sus_dram<-dram_of_taxa %>%
  mutate(
    matches = pmap_chr(across(everything()), function(...) {
      row_values <- tolower(c(...))
      found <- sus_systems[str_detect(
        paste(row_values, collapse = " "),  # combine all values into a single string
        regex(sus_systems, ignore_case = TRUE)
      )]
      if (length(found) > 0) str_c(found, collapse = ", ") else NA_character_
    })
  ) %>%
  filter(!is.na(matches))  %>%
  select(MAG,species,family,change_direction,Group,matches)

summary_sus<- sus_dram %>%
#  mutate(row_id = row_number()) %>%  # add unique row id
#  separate_rows(matches, sep = ",") %>%
  group_by(MAG,species,family,change_direction,Group,matches) %>%
  summarize(total_occurrences = n(),
#            unique_row_count = n_distinct(row_id),
            .groups = "drop")


p_sus<-ggplot(summary_sus, aes(x = reorder(species, total_occurrences), y = total_occurrences)) +
  theme_classic()+
  geom_bar(stat = "identity", aes(fill = matches)) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title = "Total Occurrences of Each Starch Utilization System",
    x = "Species",
    y = "Total Occurrences",
    fill="Sus Component"
  ) +
  theme_minimal() +
  coord_flip()
p_sus

ggsave("sus_barplots.pdf",p_sus,device = "pdf",width = 9,height = 6)
