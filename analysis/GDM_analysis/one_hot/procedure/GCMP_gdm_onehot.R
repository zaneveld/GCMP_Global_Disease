setwd("~/Documents/OSUDocs/Projects/Disease_LHS/GCMP_Global_Disease/")

#Load libraries
library(phyloseq)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(gdm)
library(qiime2R)



#Combine the two dataframes - one with the dummy variables (one_hot_encoded) and the full metadata

one_hot <- read_tsv("metadata/output/one_hot_encoded_metadata.tsv")
all_meta <- read_tsv("analysis/decontamination/input/GCMP_EMP_map.txt")
meta_1000 <- read_tsv("analysis/decontamination/output/GCMP_decontaminated_1000_sample-metadata.txt")
meta.comb <- (merge(one_hot, all_meta, on='SampleID'))
meta.comb.1000 <- merge(one_hot, meta_1000, on="SampleID")
meta.comb.1000 <- meta.comb.1000[!duplicated(meta.comb.1000$SampleID), ]
#write it as a tsv so that we can then reload it as a phyloseq object
write_tsv(meta.comb, "analysis/GDM_analysis/one_hot/output/GCMP_metadata_onehot_combined.tsv")
write_tsv(meta.comb.1000, "analysis/GDM_analysis/one_hot/output/GCMP_metadata_onehot_combined_1000.tsv")


#Create phyloseq object with decontaminated qza files
#Something wrong with the re-exported taxonomy file from decontam. 
#Using original tax reference file here because I don't think it matters
##Using all rarefied to 1000 files for this
physeq <- qza_to_phyloseq("analysis/decontamination/output/GCMP_decontaminated_1000_feature-table.qza", #feature table
                          "analysis/decontamination/output/GCMP_decontaminated_1000_tree-rooted.qza", #tree
                          "analysis/organelle_removal/output/silva_metaxa2_reference_taxonomy.qza", #taxonomy reference
                          "analysis/GDM_analysis/one_hot/output/GCMP_metadata_onehot_combined_1000.tsv") #mapping file

#Filter into just tissue
physeq.tissue <- subset_samples(physeq, tissue_compartment_T == 1)
#Could do the same for every compartment? Could loop it over or write out:
#physeq.skeleton <- subset_samples(physeq, tissue_compartment_S == 1)

##Prep for GDM##
#Create a bray curtis dissimilarity matrix & make sure the row names and column names are sample IDs
bray_dm <- as.matrix(phyloseq::distance(physeq.tissue, method="bray", type = "samples"))
View(bray_dm)
#create a metadata object and make sure the rownames are sampleIDs
metadata <- as(sample_data(physeq.tissue), "data.frame")
View(metadata)

#Check that the metadata row names and distance matrix row names are the same
#All should answer "TRUE"
nrow(bray_dm) == ncol(bray_dm)
nrow(bray_dm) == nrow(metadata)
all(rownames(bray_dm) == colnames(bray_dm))
all(rownames(bray_dm) == rownames(metadata))


#re-order the metadata to match the distance matrix order
match_metadata <- match(rownames(metadata), rownames(bray_dm))
View(match_metadata)
metadata_ord  <- metadata[match_metadata, ]
#check that the rownames were not lost
View(metadata_ord)

#Check that the metadata row names and distance matrix row names are the same
#All should answer "TRUE"
nrow(bray_dm) == ncol(bray_dm)
nrow(bray_dm) == nrow(metadata_ord)
all(rownames(bray_dm) == colnames(bray_dm))
all(rownames(bray_dm) == rownames(metadata_ord))

#All true! Woohooo!!

#Now make sure you only include data with no NAs in your predictor variables
#Here we will work with reef name, sst, depth and coordinates to trial

#remove all missing data from these categories
metadata_ord_filt <- metadata_ord %>% filter(temperature != "Missing: Not collected" & temperature != "NA" &
                                               latitude != "Missing: Not collected" & longitude != "Missing: Not collected")
metadata_ord_filt <- metadata_ord_filt %>% filter(temperature != "None")
View(metadata_ord_filt)

#Make all variables numeric - predictor variables MUST be numeric classes 
metadata_ord_filt$temperature <- as.numeric(metadata_ord_filt$temperature) #For some reason this brings in NAs, but can't see them on df
metadata_ord_filt$depth <- as.numeric(metadata_ord_filt$depth) #For some reason this brings in NAs, but can't see them on df
metadata_ord_filt$latitude <- as.numeric(metadata_ord_filt$latitude)
metadata_ord_filt$longitude <- as.numeric(metadata_ord_filt$longitude)

View(metadata_ord_filt)

#Make sure the distance matrix only includes samples in metadata
bray_dm_filt <- subset(bray_dm, rownames(bray_dm) %in% rownames(metadata_ord_filt), colnames(bray_dm) 
                       %in% rownames(metadata_ord_filt), drop = FALSE)


#Double check the number of columns and rows match and all rownames are the same for metadata and matrix
nrow(bray_dm_filt) == ncol(bray_dm_filt)
nrow(bray_dm_filt) == nrow(metadata_ord_filt)
all(rownames(bray_dm_filt) == colnames(bray_dm_filt))
all(rownames(bray_dm_filt) == rownames(metadata_ord_filt))

#Next make site-pair table
# we'll need the helper scripts provided from epiphyte workflow
source("analysis/GDM_analysis/one_hot/procedure/jld_gdm_helpers.r")
# and we'll make a list of the predictors we want to use.
# note that location data are a list, with labels Lat and Lon. Must have those
# names if location is included (but it doesn't have to be).
# NOTE - distance matrix inputs to GDM just go in this list.

predictor_list <- list(
  depth=metadata_ord_filt$depth,
  temp=metadata_ord_filt$temperature,
  complex=metadata_ord_filt$complex_robust_complex,
  robust=metadata_ord_filt$complex_robust_robust,
  width=metadata_ord_filt$colony_width1,
  location= list(Lat =metadata_ord_filt$latitude, Lon=metadata_ord_filt$longitude)
)


# 4. make site-pair table from predictor_list
spt <- site_pair_from_list(responseMat=bray_dm_filt, predList=predictor_list )

# 5. run GDM
# make geo=F if you didn't include location in predictor_list
fit1 <- gdm(spt, geo=T)
# you can also use gdm.varImp to do backward elimination, TAKES FOREVER
#variables_importance <- gdm.varImp(spt, geo=T, fullModelOnly=FALSE, parallel=TRUE, cores=4)

# 6. plot GDM
plot_gdm_jld(fit1, pred_colors="auto")

# 7. look at how much beta div our model explained (like an R^2)
fit1$explained 
