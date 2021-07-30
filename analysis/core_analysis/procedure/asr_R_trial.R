setwd("~/Documents/OSUDocs/Projects/Disease_LHS/GCMP_Global_Disease/")
library(tidyr)
library(phytools)

map <- read.delim("metadata/input/GCMP_EMP_map_r28_no_empty_samples.txt")
coral_tree <- read.tree("analysis/coral_disease_vs_adiv/input/huang_roy_molecular.newick") 

map <- map %>% filter(sample_type_EMP != "control blank" )
map <- map %>% filter(sample_type_EMP != "" )
map <- map %>% filter(host_species_id != "Missing: Not collected" )
map <- map %>% filter(host_species_id != "Not applicable" )

#head(map)
map <- map %>% separate(host_scientific_name, into = c("Genus", "Species"), sep = " ", remove = FALSE, extra = "merge")
#View(map$Genus)

unique_genera <- unique(map$Genus)

for (unique_genus in unique_genera) {
  print(unique_genus)
  genus_map <- map %>% filter(Genus == unique_genus)
  species_in_genus <- unique(genus_map$Huang_Roy_tree_name)
  if("none" %in% species_in_genus) {
    print("Skipping!")
    next
  }
  print(paste("Species in genus:", species_in_genus))
  print(paste("Tree label:", coral_tree$tip.label))
  genus_MRCA <- findMRCA(coral_tree, tips = species_in_genus)
  print(paste("Genus MRCA:", genus_MRCA))
}


trial_MRCA <- findMRCA(coral_tree, c("POC_Pocillopora"))
coral_tree$tip.label$POC_Pocillopora
View(coral_tree$tip.label)
