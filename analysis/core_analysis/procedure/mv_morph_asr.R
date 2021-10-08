library(mvMORPH)
library(tidyr)
library(dplyr)
library(phytools)

#Set all filepaths we'll need
tree.filepath <- "../input/Huang_Roy_molecular_r2.newick"
adiv.trait.table.filepath <- "../output/adiv_trait_table_tissue.tsv"
disease.trait.table.filepath <- "../input/merged_disease_table_r2_names_updated.txt"
adiv.trait.name <- "observed_features_tissue"
adiv.trait2.name<- "Total_diseased_FRRP"

adiv.tip.label.col <- "Huang_Roy_tree_name"
#Load tree
tree <- read.tree(tree.filepath)

#Check tree structure
print(paste("Binary tree:", is.binary(tree)))
min_branch_length = 1e-2
tree$edge.length[tree$edge.length==0] <- min_branch_length
na_edges <- tree$edge[is.na(tree$edge.length[tree$edge[,2]])]
tree$edge.length[na_edges] <- min_branch_length
tree$edge.length[is.na(tree$edge.length[tree$edge[,2]])] <- min_branch_length
tree$edge_length <- replace_na(tree$edge.length,min_branch_length)

#Load alpha diversity trait table
adiv.trait.table <- read.csv(adiv.trait.table.filepath,sep="\t")


#Load disease trait table
disease.trait.table <- read.csv(disease.trait.table.filepath,sep="\t")
print(colnames(adiv.trait.table))

#Produce merged disease and adiv table
adiv.trait.table <- merge(x=adiv.trait.table,y=disease.trait.table,by="Huang_Roy_tree_name",all=TRUE)

#Extract our trait of interest from the trait table
adiv.trait.data <- data.frame(adiv.trait.table[["Huang_Roy_tree_name"]],adiv.trait.table[[adiv.trait.name]],adiv.trait.table[[adiv.trait2.name]])
colnames(adiv.trait.data) <- c("Huang_Roy_tree_name",adiv.trait.name,adiv.trait2.name)

#Write the reduced table to disk at this step
write.table(adiv.trait.data, file='../output/disease_and_adiv_table.tsv', quote=FALSE, sep='\t', col.names = NA)

#Before we fit a model, we need _some_ value - even if it's just NA - for all tips
#Since the GCMP sampled only some tips, we need to map those onto a vector of all the tips on the tree
tip.labels <- data.frame(Huang_Roy_tree_name=tree$tip.label)
row.names(tip.labels) <- tip.labels$Huang_Roy_tree_name

full.trait.table <- merge(x=tip.labels,y=adiv.trait.data,
   by="Huang_Roy_tree_name",all.x=TRUE)
trait.data <- full.trait.table
row.names(trait.data) <- tip.labels$Huang_Roy_tree_name
write.table(trait.data, file='../output/disease_and_adiv_table_full_trait_table.tsv', quote=FALSE, sep='\t', col.names = NA)


#Drop the Huang Roy tree name column
trait.data <- trait.data[ , !(names(trait.data) %in% c("Huang_Roy_tree_name"))]

#print("Assigning a random extra trait column")
#trait.data[[adiv.trait.name]] <- runif(476,1,100)
#trait.data[[adiv.trait2.name]] <- runif(476,1,100)



#Filter table to tree tips
trait.data <- trait.data[rownames(trait.data) %in% tree$tip.label,]
#Drop NA rows for our metric
trait.data <- trait.data[!is.na(trait.data[,adiv.trait.name]) | !is.na(trait.data[,adiv.trait2.name]),]

write.table(trait.data, file='../output/disease_and_adiv_table_no_double_NA_values_FINAL.tsv', quote=FALSE, sep='\t', col.names = NA)


#Filter tree to table
tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% row.names(trait.data)])




print("Final data frame...")
print(trait.data)
print("About to fit the model....")

fit<-mvBM(tree,trait.data,model="BM1")
#fit<-mvOU(tree,trait.data,model="OU1",method="inverse",optimization="subplex",param = list(decomp="cholesky"))

tip.label.length <- length(tree$tip.label)
trait.data.length <- length(trait.data)
print(paste("Tip label length:",tip.label.length))
print(paste("Trait data length:",trait.data.length))

#fit<-phylo.impute(tree,trait.data)
print(fit)
print("Length of fit:")
print(length(fit))
print("Attributes of fit:")
print(attributes(fit))

print("About to run estim")
imp <- estim(tree,trait.data,fit,asr=TRUE)
print("Done with imputation")
print(imp$estimate)
print("I can has trait estimates?")


map <- read.table("../input/GCMP_EMP_map_r29.txt",sep="\t",header=TRUE,fill=TRUE)
coral_tree <- tree

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
