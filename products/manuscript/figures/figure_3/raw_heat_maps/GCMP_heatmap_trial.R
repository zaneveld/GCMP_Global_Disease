###GCMP Heat map for publication###

#The following code is for running a heat map against the following:
# 1) %disease (perc_dis)
#2) growth rate (growth_rate_mm_per_year)
#3) microbial alpha div (observed_features_tissue)
#4) microbial dominance (dominance_tissue)
#5) endozoicomonas dominance (endo_relabund)
#6) alpha vs gamma proteobacteria domination

#First load libraries

library(phytools)
library(RColorBrewer)
library(dplyr)

#Load in a tree file and trait table and filter the trait table into the columns intended for the heat map

#We are going to do everything across the 3 compartments: coral mucous, tissue, & skeleton

tree <- read.tree('input/huang_roy_genus_tree.newick')
meta <- read.csv('input/GCMP_trait_table.csv', header = TRUE)

#Bring in the growth data
growth <- read.csv("input/GCMP_wGrowth.csv", header = TRUE)
traits <- merge(meta, growth, all = TRUE, by = "host_genus_id")

#Bring in endozoicomonas data
endo <- read.csv("input/GCMP_endo.csv", header = TRUE)
traits <- merge(traits, endo, all = TRUE, by = "host_genus_id")


heatmap.traits <- traits %>% select("host_genus_id", "perc_dis", "growth_rate_mm_per_year", "observed_features_tissue", "dominance_tissue", "tissue_endo", 
                                    "observed_features_mucus", "dominance_mucus", "mucus_endo", 
                                    "observed_features_skeleton", "dominance_skeleton", "skeleton_endo",
                                    "most_abundant_class_tissue", "most_abundant_class_mucus", "most_abundant_class_skeleton")
rownames(heatmap.traits) <- heatmap.traits$host_genus_id
heatmap.traits <- heatmap.traits[,c(2:15)]

#Remove rows that have all NAs
heatmap.traits <- heatmap.traits[rowSums(is.na(heatmap.traits)) != ncol(heatmap.traits), ]
head(heatmap.traits)

#drop all rows in which microbiome data is unavailable
library(tidyr)
heatmap.traits.micro.only <- heatmap.traits %>% drop_na(observed_features_tissue)

#Run a min-max normalization to standardize all columns

my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
heatmap.traits.micro.only.st <- heatmap.traits.micro.only
heatmap.traits.micro.only.st$perc_dis <- heatmap.traits.micro.only.st$perc_dis/my.max(heatmap.traits.micro.only.st$perc_dis)
heatmap.traits.micro.only.st$growth_rate_mm_per_year <- heatmap.traits.micro.only.st$growth_rate_mm_per_year/my.max(heatmap.traits.micro.only.st$growth_rate_mm_per_year)
heatmap.traits.micro.only.st$observed_features_tissue <- heatmap.traits.micro.only.st$observed_features_tissue/my.max(heatmap.traits.micro.only.st$observed_features_tissue)
heatmap.traits.micro.only.st$observed_features_mucus <- heatmap.traits.micro.only.st$observed_features_mucus/my.max(heatmap.traits.micro.only.st$observed_features_mucus)
heatmap.traits.micro.only.st$observed_features_skeleton <- heatmap.traits.micro.only.st$observed_features_skeleton/my.max(heatmap.traits.micro.only.st$observed_features_skeleton)
heatmap.traits.micro.only.st$dominance_tissue <- heatmap.traits.micro.only.st$dominance_tissue/my.max(heatmap.traits.micro.only.st$dominance_tissue)
heatmap.traits.micro.only.st$dominance_mucus <- heatmap.traits.micro.only.st$dominance_mucus/my.max(heatmap.traits.micro.only.st$dominance_mucus)
heatmap.traits.micro.only.st$dominance_skeleton <- heatmap.traits.micro.only.st$dominance_skeleton/my.max(heatmap.traits.micro.only.st$dominance_skeleton)
heatmap.traits.micro.only.st$tissue_endo <- heatmap.traits.micro.only.st$tissue_endo/my.max(heatmap.traits.micro.only.st$tissue_endo)
heatmap.traits.micro.only.st$mucus_endo <- heatmap.traits.micro.only.st$mucus_endo/my.max(heatmap.traits.micro.only.st$mucus_endo)
heatmap.traits.micro.only.st$skeleton_endo <- heatmap.traits.micro.only.st$skeleton_endo/my.max(heatmap.traits.micro.only.st$skeleton_endo)
summary(heatmap.traits.micro.only.st)
#Check that Max is 1 for every column

#Prune the tree to micro only taxa
genera.to.keep <- rownames(heatmap.traits.micro.only.st)
#Prune tree to match trait table
pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, genera.to.keep))

#Make a colour scheme for the heta map for each comp
#disease_ramp_color_palette <- colorRampPalette(c("lemonchiffon1", "red"))(100)
disease_ramp_color_palette <- colorRampPalette(c("#E7E1E1", "#601D1F"))(100)

mucus_ramp_color_palette <- colorRampPalette(c("#C7F2F2", "#098F94"))(100)
tissue_ramp_color_palette <- colorRampPalette(c("#F7D8CB", "#E98A5C")) (100)
skeleton_ramp_color_palette <- colorRampPalette(c("#E5C5EA", "#8D3F95")) (100)

#Run the heatmap code on tissue

#pull out standardized tissue traits
heatmap.traits.tissue <- heatmap.traits.micro.only.st[,c(1:5)]
#pull out gamma and alpha classes (yes/no = 1/0) - can change colour scheme in illustrator for these
heatmap.traits.tissue$gamma_class <- ifelse(heatmap.traits.micro.only.st$most_abundant_class_tissue == "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria", 1, 0)
heatmap.traits.tissue$alpha_class <- ifelse(heatmap.traits.micro.only.st$most_abundant_class_tissue == "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria", 1, 0)

#Run the following code all together to save
pdf(file="../output/heatmap_tissue.pdf")
#phylo.heatmap(pruned.tree, heatmap.traits, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
phylo.heatmap(pruned.tree, heatmap.traits.tissue, split = c(1,0.5), standardize=FALSE, colors = tissue_ramp_color_palette, fsize=0.3)
dev.off()

##pull out standardized mucus traits

heatmap.traits.mucus <- heatmap.traits.micro.only.st[,c(1:2,6:8)]
#pull out gamma and alpha classes (yes/no = 1/0) - can change colour scheme in illustrator for these
heatmap.traits.mucus$gamma_class <- ifelse(heatmap.traits.micro.only.st$most_abundant_class_mucus == "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria", 1, 0)
heatmap.traits.mucus$alpha_class <- ifelse(heatmap.traits.micro.only.st$most_abundant_class_mucus == "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria", 1, 0)

pdf(file="../output/heatmap_mucus.pdf")
#phylo.heatmap(pruned.tree, heatmap.traits, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
phylo.heatmap(pruned.tree, heatmap.traits.mucus, split = c(1,0.5), standardize=FALSE, colors = mucus_ramp_color_palette, fsize=0.3)
dev.off()


##Skeleton
heatmap.traits.skeleton <- heatmap.traits.micro.only.st[,c(1:2,9:11)]
#pull out gamma and alpha classes (yes/no = 1/0) - can change colour scheme in illustrator for these
heatmap.traits.skeleton$gamma_class <- ifelse(heatmap.traits.micro.only.st$most_abundant_class_skeleton == "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria", 1, 0)
heatmap.traits.skeleton$alpha_class <- ifelse(heatmap.traits.micro.only.st$most_abundant_class_skeleton == "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria", 1, 0)

pdf(file="../output/heatmap_skeleton.pdf")
#phylo.heatmap(pruned.tree, heatmap.traits, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
phylo.heatmap(pruned.tree, heatmap.traits.skeleton, split = c(1,0.5), standardize=FALSE, colors = skeleton_ramp_color_palette, fsize=0.3)
dev.off()


#In case we want to pull out the disease/growth rate data as a different colour
heatmap.traits.disgr <- heatmap.traits.micro.only.st[,c(1:2)]
pdf(file="output/heatmap_disease_growth_only.pdf")
#phylo.heatmap(pruned.tree, heatmap.traits, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
phylo.heatmap(pruned.tree, heatmap.traits.disgr, split = c(1,0.5), standardize=FALSE, colors = disease_ramp_color_palette, fsize=0.3)
dev.off()


