#Goal: This script identifies the most recent common ancestor (MRCA) of each
#genus in the analysis to help with data merging at the genus level.

#Uncomment to install packages
#install.packages("tidyr",repos='http://cran.us.r-project.org')
#install.packages("phytools",repos='http://cran.us.r-project.org')
#install.packages("dplyr",repos='http://cran.us.r-project.org')

library(dplyr)
library(tidyr)
library(phytools)

map <- read.table("../input/GCMP_EMP_map_r29.txt",sep="\t",header=TRUE,fill=TRUE)
coral_tree <- read.tree("../input/huang_roy_molecular.newick") 


#Filter out samples with no binomial name
map <- map %>% filter(sample_type_EMP != "control blank" )
map <- map %>% filter(sample_type_EMP != "" )
map <- map %>% filter(host_species_id != "Missing: Not collected" )
map <- map %>% filter(host_species_id != "Not applicable" )

#Split scientific name into genus and species, filling missing with NA
map <- map %>% separate("host_scientific_name", into = c("Genus", "Species"), sep = " ", remove = FALSE, extra = "merge")

#Collect the unique genus names to iterate over
unique_genera <- unique(map$Genus)

#make a vector mapping genus to node number
genus_mrcas = c()
genus_names = c()
for (unique_genus in unique_genera) {
  print("************") 
  print(unique_genus)
  genus_map <- map %>% filter(Genus == unique_genus)
  species_in_genus <- unique(genus_map$Huang_Roy_tree_name)
  
  #skip genera that are in the GCMP mapping file but 
  #not in the Huang Roy tree
  if("none" %in% species_in_genus) {
    print("Skipping!")
    next
  }

  #We need special handling of genera with one species
  #in the dataset. These cause an error in findMRCA.
  #To avoid that we set the species as representative of the 
  #genus in these cases. 
  print(paste("Species in this genus:",length(species_in_genus)))
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
  
  #Add the results to our vectors of genera and their names
  genus_mrcas <- c(genus_mrcas,genus_MRCA)
  genus_names <- c(genus_names,unique_genus)
}

i <- 1
for (genus_name in genus_names) {
    print(paste("Genus:",genus_name))
    print(paste("Genus MRCA node name:",genus_mrcas[i]))
    i <- i + 1
}

#Write results to file
df <- data.frame(genus_mrcas,row.names=genus_names)
write.csv(df,'../output/genus_to_mrca_mapping.csv')
