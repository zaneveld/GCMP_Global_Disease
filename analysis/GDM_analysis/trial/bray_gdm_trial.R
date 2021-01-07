##Trial Generalized Dissimilarity Model (GDM) 
#Using bray curtis distance matrix from .qza output from QIIME2

#Set working directory

setwd("~/Documents/OSU_Documents/Projects/Disease_LHS/gdm_analysis/gdm_disease/")

#Load libraries
library(phyloseq)
library(qiime2R)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(gdm)


#Create phyloseq object from feature table
metadata <- read_tsv("GCMP_EMP_map_r28_no_empty_samples.txt")
View(metadata)
bray <- read_qza("bdiv_greengenes_metaxa2_tissue_braycurtis_1000.qza") 
View(bray)

#create a distance matrix & make sure the names are correct across columns and rows
bray_dm <- as.matrix(bray$data)
View(bray_dm)
rownames(bray_dm) <- colnames(bray_dm)
rownames(metadata) <- metadata$SampleID
View(metadata)

#subset the metadata to match the otu table
metadata_tissue <- subset(metadata, SampleID %in% rownames(bray_dm))
View(metadata_tissue)
#rename the rownames as sample ID
rownames(metadata_tissue) <- metadata_tissue$SampleID


#re-order the metadata to match the distance matrix order
match_metadata <- match(rownames(bray_dm), rownames(metadata_tissue))
match_metadata
metadata_ord  <- metadata_tissue[match_metadata, ]
View(metadata_ord)
#check that the rownames were not lost
rownames(metadata_ord) <- metadata_ord$SampleID

#Check that the metadata row names and distance matrix row names are the same
#All should answer "TRUE"
nrow(bray_dm) == ncol(bray_dm)
nrow(bray_dm) == nrow(metadata_ord)
all(rownames(bray_dm) == colnames(bray_dm))
all(rownames(bray_dm) == rownames(metadata_ord))

#All true!

#Now make sure you only include data with no NAs in your predictor variables
#Here we will work with reef name, sst, depth and coordinates to trial

#remove all missing data from these categories
metadata_ord_filt <- metadata_ord %>% filter(reef_name != "Missing: Not collected" & surface_temperature != "Missing: Not collected" &
                                               latitude != "Missing: Not collected" & longitude != "Missing: Not collected")
View(metadata_ord_filt)

#Make all variables numeric - predictor variables MUST be numeric classes 
metadata_ord_filt$surface_temperature <- as.numeric(metadata_ord_filt$surface_temperature)
metadata_ord_filt$depth <- as.numeric(metadata_ord_filt$depth)
metadata_ord_filt$latitude <- as.numeric(metadata_ord_filt$latitude)
metadata_ord_filt$longitude <- as.numeric(metadata_ord_filt$longitude)

#Make sure the distance matrix only includes samples in metadata
bray_dm_filt <- subset(bray_dm, rownames(bray_dm) %in% rownames(metadata_ord_filt))
match_bray <- match(rownames(bray_dm_filt), colnames(bray_dm_filt))
match_bray
bray_dm_filt_2  <- bray_dm_filt[,match_bray]

#Double check the number of columns and rows match and all rownames are the same for metadata and matrix
nrow(bray_dm_filt_2) == ncol(bray_dm_filt_2)
nrow(bray_dm_filt_2) == nrow(metadata_ord_filt)
all(rownames(bray_dm_filt_2) == colnames(bray_dm_filt_2))
all(rownames(bray_dm_filt_2) == rownames(metadata_ord_filt))

#Next make site-pair table
# we'll need the helper scripts provided from epiphyte workflow
source("jld_gdm_helpers.r")
# and we'll make a list of the predictors we want to use.
# note that location data are a list, with labels Lat and Lon. Must have those
# names if location is included (but it doesn't have to be).
# NOTE - distance matrix inputs to GDM just go in this list.

predictor_list <- list(
  #host_genus=metadata_ord_filt$host_genus,
  #reef=metadata_ord_filt$reef_name, 
  temp=metadata_ord_filt$surface_temperature,
  depth=metadata_ord_filt$depth,
  location= list(Lat =metadata_ord_filt$latitude, Lon=metadata_ord_filt$longitude)
  )


# 4. make site-pair table from predictor_list
spt <- site_pair_from_list(responseMat=bray_dm_filt_2, predList=predictor_list )

# 5. run GDM
# make geo=F if you didn't include location in predictor_list
fit1 <- gdm(spt, geo=T)
# you can also use gdm.varImp to do backward elimination, TAKES FOREVER
# variables_importance <- gdm.varImp(spt, geo=T, fullModelOnly=FALSE, parallel=TRUE, cores=4)

# 6. plot GDM
plot_gdm_jld(fit1, pred_colors="auto")


# 7. look at how much beta div our model explained (like an R^2)
fit1$explained
# wow it works pretty well when you make fake data!

