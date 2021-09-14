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
map <- map %>% filter(Huang_Roy_tree_name != "Missing: Not collected" )
map <- map %>% filter(Huang_Roy_tree_name != "Not applicable" )
map <- map %>% filter(Huang_Roy_tree_name != "" )
map <- map %>% filter(Huang_Roy_tree_name != "none" )

#Split scientific name into genus and species, filling missing with NA
map <- map %>% separate("Huang_Roy_tree_name", into = c("Family_Code","Genus", "Species"), sep = "_", remove = FALSE, extra = "merge")

map <- map %>% filter(!is.na(map$Genus))
print(paste("Genus column:",map$Genus))
print(paste("Species column:",map$Species))
#Collect the unique genus names in the GCMP metadata to iterate over
unique_genera <- unique(map$Genus)

#Collect unique genus names on the tree

coral_tree_tip_names <- coral_tree$tip.label
coral_tree_tip_df <- as.data.frame(coral_tree_tip_names) 
coral_tree_tip_df <- coral_tree_tip_df %>% separate("coral_tree_tip_names", into = c("Family_Code","Genus", "Species"), sep = "_",remove = FALSE, extra = "merge")
unique_genera_on_tree <- unique(coral_tree_tip_df$Genus)

#NOTE: there are some GCMP samples identified only to the genus level.
#At the same time, there are Huang Roy tree nodes identified only to the genus level (POR_Porites)
#Currently these are retained and mapped together. This is a caveat to the weighting of species
#within genera BUT will almost always be better than simple averaging across the tree, and retains
#maximum 16S data (e.g. we expect that having 6 vs. 9 samples is worse for the final estimate
#than retaining 9 and risking that 3 of them are misplaced on the tree, especially given that
#tree relationships are inexact). This has no effect in genera with just one species in the dataset

#make a vector mapping genus to node number
genus_mrcas = c()
genus_names = c()
gcmp_species = c()
for (unique_genus in unique_genera_on_tree) {
  print("************") 
  print(unique_genus)
  genus_map <- map %>% filter(Genus == unique_genus)
  genus_map_tree <- coral_tree_tip_df %>% filter(Genus == unique_genus)
  
  species_in_genus <- unique(genus_map$Huang_Roy_tree_name)
  species_in_genus_on_tree <- unique(genus_map_tree$coral_tree_tip_names)

  if (length(species_in_genus) == 0) {
      print("Skipping ... no data for this genus")
      next
  }
  print(paste("Species in genus:",species_in_genus))  
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
  #print(paste("Species in this genus:",length(species_in_genus)))
  if (length(species_in_genus) == 1){
    print("only 1 species in this genus!")
    print(paste("Species in genus:",species_in_genus))
    
    genus_MRCA <- coral_tree$tip.label[match(species_in_genus,coral_tree$tip.label)]
    print(paste("Genus MRCA:", genus_MRCA))
    #Alternative strategy if we got an NA
    #In cases where we have e.g. POR_Goniopora (incomplete id) and Huang Roy has
    # two specific species of Goniopora, we find all tip labels that have Goniopora in
    # them and then take their MRCA?
    if (is.na(genus_MRCA)){
        genus_MRCA <- findMRCA(coral_tree,tips = species_in_genus_on_tree)
    }
  }
  else{
      print(paste("Species in genus on tree:",species_in_genus_on_tree))
      genus_MRCA <- findMRCA(coral_tree, tips = species_in_genus_on_tree,type="node")
      print(paste("Genus MRCA:", genus_MRCA))
  }
  #Add the results to our vectors of genera and their names
  genus_mrcas <- c(genus_mrcas,genus_MRCA)
  genus_names <- c(genus_names,unique_genus)
  gcmp_species <- c(gcmp_species,length(species_in_genus))
  print(genus_mrcas)
  print(length(genus_mrcas))
  print(genus_names)
  print(length(genus_names))
  
}
print("Summarizing results...")
i <- 1
for (genus_name in genus_names) {
    print(paste("Genus:",genus_name))
    print(paste("Genus MRCA node name:",genus_mrcas[i]))
    print(paste("GCMP species in this genus:",gcmp_species[i]))
    i <- i + 1
}

#Write results to file
df <- data.frame(genus_mrcas,row.names=genus_names)
write.csv(df,'../output/genus_to_mrca_mapping.csv')
