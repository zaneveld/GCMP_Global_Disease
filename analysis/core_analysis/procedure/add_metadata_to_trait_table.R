###Adding host taxonomy and other host-associated metadata to trait table for easy subsetting downstream###

#Set working directory (unnecessary if working from command line in the core analysis/procedure directory)
setwd("~/Documents/OSUDocs/Projects/Disease_LHS/GCMP_Global_Disease/analysis/core_analysis/procedure")

#Load libraries
library(plyr)
library(dplyr)
library(tidyr)

##Load the original mapping file
map <- read.delim("../input/GCMP_EMP_map_r28_no_empty_samples.txt")

#Identify the variables that we want to add as metadata & make sure they are factors
map$country <- as.factor(map$country)
map$ocean <- as.factor(map$ocean)
map$ocean_area <- as.factor(map$ocean_area)
map$shelf_location <- as.factor(map$shelf_location)
map$binary_turf_contact <- as.factor(map$binary_turf_contact)

#Take those factors and make them binary, so that we can apply them to genus level (not just species)
map$genus_in_australia <- ifelse(map$country == "Australia", 1, 0)
map$genus_in_colombia <- ifelse(map$country == "Colombia", 1, 0)
map$genus_in_curacao <- ifelse(map$country == "Curacao", 1, 0)
map$genus_in_france <- ifelse(map$country == "France", 1, 0)
map$genus_in_panama <- ifelse(map$country == "Panama", 1, 0)
map$genus_in_saudiarabia <- ifelse(map$country == "Saudi Arabia", 1, 0)
map$genus_in_singapore <- ifelse(map$country == "Singapore", 1, 0)

map$genus_in_atlantic <- ifelse(map$ocean == "Atlantic", 1, 0)
map$genus_in_indian <- ifelse(map$ocean == "Indian", 1, 0)
map$genus_in_pacific <- ifelse(map$ocean == "Pacific", 1, 0)

map$genus_in_caribbean <- ifelse(map$ocean_area == "Caribbean", 1, 0)
map$genus_in_coralsea <- ifelse(map$ocean_area == "Coral Sea", 1, 0)
map$genus_in_easternindian <- ifelse(map$ocean_area == "Eastern Indian", 1, 0)
map$genus_in_easternpacific <- ifelse(map$ocean_area == "Eastern Pacific", 1, 0)
map$genus_in_redsea <- ifelse(map$ocean_area == "Red Sea", 1, 0)
map$genus_in_southchinasea <- ifelse(map$ocean_area == "South China Sea", 1, 0)
map$genus_in_tasmansea <- ifelse(map$ocean_area == "Tasman Sea", 1, 0)
map$genus_in_westernindian <- ifelse(map$ocean_area == "Western Indian", 1, 0)

map$genus_in_aquariumdeep <- ifelse(map$shelf_location == "aquarium-deepsea", 1, 0)
map$genus_in_inshore <- ifelse(map$shelf_location == "inshore", 1, 0)
map$genus_in_midshelf <- ifelse(map$shelf_location == "midshelf", 1, 0)
map$genus_in_offshore <- ifelse(map$shelf_location == "offshore", 1, 0)

map$genus_in_turfcontact <- ifelse(map$binary_turf_contact == c("y", "y "), 1, 0 ) 
#apparently there is one factor level that is a "y " with an added space^. Add this is by using the concatenate function


#Extract host genus ID & taxonomy columns & other metadata necessary for subset
map.meta <- map %>% select("host_genus_id", 
                          "host_clade_sensu_fukami_numeric","host_clade_sensu_fukami",
                          "taxonomy_string_to_family", "taxonomy_string_to_order", "taxonomy_string_to_subclass", 
                          "taxonomy_string_to_class", "complex_robust", "functional_group_sensu_darling",
                          "genus_in_australia", "genus_in_curacao", "genus_in_france", "genus_in_panama",
                          "genus_in_saudiarabia", "genus_in_singapore", "genus_in_atlantic", "genus_in_indian",
                          "genus_in_pacific", "genus_in_caribbean", "genus_in_coralsea", "genus_in_easternindian",
                          "genus_in_easternpacific", "genus_in_redsea", "genus_in_southchinasea", "genus_in_tasmansea",
                          "genus_in_westernindian", "genus_in_aquariumdeep", "genus_in_inshore", "genus_in_midshelf",
                          "genus_in_offshore", "genus_in_turfcontact")


##subset to just the unique genera
map.genus.unique <- map.meta[!duplicated(map.meta$host_genus_id), ]  

#Load the trait data
trait.table <- read.delim("../output/GCMP_trait_table_genus.tsv")
#Rename the variable "host_genus" to "host_genus_id" to match the mapping file above:
trait.table <- rename(trait.table, host_genus_id = host_genus) 

#R-Join the two datasets so that taxonomic info is added according to Huang Roy tree name
trait.table.meta <- right_join(map.genus.unique, trait.table, by = "host_genus_id")
#Write a CSV
write.csv(trait.table.meta, "../output/GCMP_trait_table_meta.csv")
#done!

