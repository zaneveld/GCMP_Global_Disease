### R Script for Pre-Processing Seq Data for GCMP Disease/LHS Manuscript ###

#load libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(decontam)

#Set working directory
setwd("~/Documents/OSUDocs/Projects/Disease_LHS/GCMP_Global_Disease/analysis")

#use qiime2R to upload data into a phyloseq object
#mapping file has "#" in it - I replaced all "#" with "." in order to upload
#this new map is called "GCMP_EMP_map.txt" found in input file 
#For trial use the following code


#For final use the following code
physeq <- qza_to_phyloseq("organelle_removal/output/effects_of_rarefaction_analysis/feature_table_silva_metaxa2_all.qza", #feature table
                          "phylogeny_insertion/output/insertion-tree_silva_GCMP.qza", #tree
                          "organelle_removal/output/silva_metaxa2_reference_taxonomy.qza", #taxonomy reference
                          "decontamination/input/GCMP_EMP_map.txt") #mapping file


#Check the taxonomic classification is correct
rank_names(physeq)
tax_table(physeq)


#make a sample data frame
sample.data <- as(sample_data(physeq), "data.frame") #423 observations, 170 variables

#Using the package Decontam remove any potential contaminants found in negative/blank controls
#First inspect library size 
sample.data$LibrarySize <- sample_sums(physeq)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = sample_type_EMP)) +
  geom_point()
ggsave("output/SILVA_library_size.pdf", plot= last_plot())

#identify the blank samples by creating a new column called is.neg (provides TRUE/FALSE)
sample_data(physeq)$is.neg <- sample_data(physeq)$sample_type_EMP == "control blank"
sample_data(physeq)$is.neg
#Next check for contaminants using prevalence and threshold 0.5
#The 0.5 probability threshold will identify all sequences that are more prevalent in negative controls 
#than in positive samples. This is more aggressive and conservative. Default is 0.1, which is less aggressive.
contamdf.prev <- isContaminant(physeq, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant) # 662 contaminant ASVs at a level of 0.5, 238 at 0.1 
head(which(contamdf.prev$contaminant))

#Look at the # of times several of these taxa were observed in positive and negative samples
ps <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.neg <- prune_samples(sample_data(ps)$is.neg == "TRUE", ps)
ps.pos <- prune_samples(sample_data(ps)$is.neg == "FALSE", ps)
# Make data.frame of prevalence in positive and negative samples
df.ps <- data.frame(ps.pos=taxa_sums(ps.pos), ps.neg=taxa_sums(ps.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.ps, aes(x=ps.neg, y=ps.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("output/SILVA_prevalence_plot_0.5.pdf", plot = last_plot())

#Remove contaminants
physeq.noncont <- prune_taxa(!contamdf.prev$contaminant, physeq)

#Check final numbers
tax_table(physeq.noncont) 
sample_data(physeq.noncont) #1383 samples total, 171 variables
#Count Sequences
sum(sample_sums(physeq.noncont)) #37,469,008

saveRDS(physeq.noncont, "output/physeq_noncont.RDS")
physeq.noncont <- readRDS("output/physeq_noncont.RDS")
#Using function phyloseq2qiime2 written by Dr. Christian Edwardson
#To go from a phyloseq object to biom, newick and text files
#This function uses the following packages that must be installed or you will get an error:
#(phyloseq)
#(biomformat)
#(ape)
#(Biostrings)
#(dada2)
source("phyloseq2QIIME2.R")

#could set your working directory to your output folder for this
setwd("~/Documents/OSUDocs/Projects/Disease_LHS/Decontamination/output/")
phyloseq2qiime2(physeq.noncont)
#Saves output into working directory. 
#this code is currently having issues writing a newick file (getting a C stack error)
#Try the following code instead from the castor package:
library(castor)
write_tree(phy_tree(physeq.noncont),"physeq.noncont_tree-rooted.newick", append = FALSE)
#worked!
#Now just need to use QIIME2 to convert biom and newick back to qza files
#See biom_tre-to-qza.sh 
