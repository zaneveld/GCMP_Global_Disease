
#install.packages("tidyr",repos='http://cran.us.r-project.org')
#install.packages("phytools",repos='http://cran.us.r-project.org')
#install.packages("dplyr",repos='http://cran.us.r-project.org')
#install.packages("stringr", repos = 'http://cran.us.r-project.org')

library(dplyr)
library(tidyr)
library(phytools)
library(stringr)

map <- read.table("../input/GCMP_EMP_map_r29.txt",sep="\t",header=TRUE,fill=TRUE)
coral_tree <- read.tree("../input/huang_roy_molecular.newick") 

mrca_full <- mrca(coral_tree, full = FALSE)
View(mrca_full)
write.csv(mrca_full, "../output/huang_roy_mrca.csv")

map <- map %>% filter(sample_type_EMP != "control blank" )
map <- map %>% filter(sample_type_EMP != "" )
map <- map %>% filter(host_species_id != "Missing: Not collected" )
map <- map %>% filter(host_species_id != "Not applicable" )
map <- map %>% filter(str_count(Huang_Roy_tree_name, "_") > 1)

print(map)
View(map$Huang_Roy_tree_name)
map <- map %>% separate("host_scientific_name", into = c("Genus", "Species"), sep = " ", remove = FALSE, extra = "merge")
#View(map$Genus)


unique_genera <- unique(map$Genus)

for (unique_genus in unique_genera) {
  print(paste("Unique genus:", unique_genus))
  genus_map <- map %>% filter(Genus == unique_genus)
  species_in_genus <- unique(genus_map$Huang_Roy_tree_name)
  if("none" %in% species_in_genus) {
    print("Skipping!")
    next
  }
  #print(paste("Tree label:", coral_tree$tip.label))
  print(paste("Species in genus:", species_in_genus))
  genus_MRCA <- findMRCA(coral_tree, tips = species_in_genus)
  print(paste("Genus MRCA:", genus_MRCA))
}

findMRCA(coral_tree, tips = c("FAV_Cyphastrea_serailia", "FAV_Cyphastrea_micropthalma"))


trial_MRCA <- findMRCA(coral_tree, c("POC_Pocillopora"))
coral_tree$tip.label$POC_Pocillopora
View(coral_tree$tip.label)

