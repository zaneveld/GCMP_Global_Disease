#install.packages("tidyr",repos='http://cran.us.r-project.org')
#install.packages("phytools",repos='http://cran.us.r-project.org')
#install.packages("dplyr",repos='http://cran.us.r-project.org')

library(dplyr)
library(tidyr)
library(phytools)

map <- read.table("../input/GCMP_EMP_map_r29.txt",sep="\t",header=TRUE,fill=TRUE)
coral_tree <- read.tree("../input/huang_roy_molecular.newick") 



map <- map %>% filter(sample_type_EMP != "control blank" )
map <- map %>% filter(sample_type_EMP != "" )
map <- map %>% filter(host_species_id != "Missing: Not collected" )
map <- map %>% filter(host_species_id != "Not applicable" )

print(map)
map <- map %>% separate("host_scientific_name", into = c("Genus", "Species"), sep = " ", remove = FALSE, extra = "merge")
#View(map$Genus)

mrca_full <- mrca(coral_tree, full = FALSE)
write.csv(mrca_full, "../output/huang_roy_full_mrca_table.csv")

unique_genera <- unique(map$Genus)
#make a vector mapping genus to node number
genus_mrcas = c()
genus_names = c()
for (unique_genus in unique_genera) {
  print("************") 
  print(unique_genus)
  genus_map <- map %>% filter(Genus == unique_genus)
  species_in_genus <- unique(genus_map$Huang_Roy_tree_name)
  if("none" %in% species_in_genus) {
    print("Skipping!")
    next
  }
  print(length(species_in_genus))
  if (length(species_in_genus) == 1){
    print("only 1 species in this genus!")
    print(paste("Species in genus:",species_in_genus))
    genus_MRCA <- coral_tree$tip.label[match(species_in_genus,coral_tree$tip.label)]
    print(paste("Genus MRCA:", genus_MRCA))
    
  }
  else{
    genus_MRCA <- findMRCA(coral_tree, tips = species_in_genus,type="node")
    print(paste("Genus MRCA:", genus_MRCA))
  }
  genus_mrcas <- c(genus_mrcas,genus_MRCA)
  genus_names <- c(genus_names,unique_genus)
}

i <- 1
for (genus_name in genus_names) {
    print(paste("Genus:",genus_name))
    print(paste("Genus MRCA node name:",genus_mrcas[i]))
    i <- i + 1
}

df <- data.frame(genus_mrcas,row.names=genus_names)
write.csv(df,'../output/genus_to_mrca_mapping.csv')
