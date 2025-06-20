library(tidyverse)
library(VennDiagram)
#library(UpSetR)
#library(forcats)

setwd("~/database_comparisons") #set to this directory


# input_file <- "CBAJ_representative_genomes.tsv"
# input_file <- "CMMG_representative_genomes.tsv"
# input_file <- "MGBC_representative_genomes.tsv"
# input_file <- "iMGMC_representative_genomes.tsv"
# input_file <- "All_databases_representative_genomes.tsv"
# database_name<- "All_databases"

#database_name<- str_split(input_file, "_",n=2)[[1]][1]


our_dataset_clusters <- function(database_name) {
  input_file=paste0(database_name,"_representative_genomes.tsv")
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
    ungroup()
  
  genome_origins <- clusters %>%
    pivot_longer(cols = -Representative, names_to = "Type", values_to = "Genome") %>%
    filter(Genome != "") %>%  # Remove empty entries
    mutate(Database = ifelse(grepl("Mouse|maxbin", Genome), "Mohssen_etal", database_name)) %>%
    select(Representative, Genome, Database)%>%
    filter(str_detect(Representative, "Mouse|maxbin"))%>%
    group_by(Representative, Database) %>%
    summarise(Present = 1, .groups = "drop") %>%
    pivot_wider(names_from = Database, values_from = Present, values_fill = list(Present = 0))
  return(genome_origins)
  
}

CMMG<- our_dataset_clusters("CMMG") 
MGBC<- our_dataset_clusters("MGBC")
iMGMC<- our_dataset_clusters("iMGMC")


union_integrated_dbs<-full_join(CMMG,MGBC) %>%
  full_join(iMGMC) %>%
  select(-Mohssen_etal) %>%
  pivot_longer(cols=-Representative,names_to = "Database", values_to = "Presence") %>%
  group_by(Database,Presence) %>%
  summarise(Count = n(), .groups = "drop") %>%
  filter(!is.na(Presence)) %>%
  mutate(Presence = case_when(
    Presence == 1 ~ "Improved MAGs in our database",
    Presence == 0 ~ "Unique MAGs to our database",
  ))
discrete_palettes <- list(c("#81b29a",
                            "#e07a5f"))


db_size <- tibble(
  Database = c("CMMG", "iMGMC", "MGBC"),
  size = c(1573, 1296 , 1091))
text_positions <- union_integrated_dbs %>%
  group_by(Database) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  inner_join(db_size) %>%
  mutate(percent=paste0(round(Total/size*100, 2), "%"))

barplot<- ggplot(union_integrated_dbs, 
       aes(x = fct_reorder(Database, Count, .desc = TRUE),
                                 y = Count, 
                                 fill = Presence)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size = 4) +
  geom_text(data = text_positions, aes(x = Database, y = Total, label = percent), 
             vjust = -0.5, size = 4, color= "#3d405b", fontface= "bold",inherit.aes = FALSE) +
  scale_fill_discrete(type = discrete_palettes)+
  scale_y_continuous(breaks = seq(0, max(union_integrated_dbs$Count) + 80, by = 30)) +
  theme_classic() +
  labs(title = "Pairwise Clustering with Integrated Databases",
       x = "Database",
       y = "Number of MAGs",
       fill = NULL) +
  theme(
    legend.position = c(0.7, 0.9),
    legend.justification = "center"  # Centers it
  )+
  guides(fill = guide_legend(ncol = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  
barplot
pdf("stacked_barplot_pairwise_clustering_integrated_dbs_only.pdf", width = 5, height = 5)
print(barplot)
dev.off()
