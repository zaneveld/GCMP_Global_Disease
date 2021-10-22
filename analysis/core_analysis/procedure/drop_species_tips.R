library(ape)
library(phytools)
library(phangorn)

sp_tree <- read.tree('../input/huang_roy_molecular.newick')

# main

tips<-sp_tree$tip.label
genera<-unique(sapply(strsplit(tips,"_"), function(x) x[2]))
ii <- sapply(genera, function(x,y) grep(x,y)[1], y = tips)
genus_tree <- drop.tip(sp_tree, setdiff(sp_tree$tip.label, tips[ii]))
genus_tree$tip.label <- sapply(strsplit(genus_tree$tip.label,"_"), function(x) x[2])
plotTree(genus_tree,ftype='i')
write.tree(genus_tree, file="../output/huang_roy_genus_tree.newick")