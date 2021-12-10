# Ancestral State Reconstruction of alpha diversity 
library(phytools)

ctree <- read.tree('../../input/pruned_tree.newick')
trait_table <- read.table('../../input/trait_table_r11_by_disease_adiv_rank_zeroed.txt', sep="\t", header = TRUE, row.names = 1)

adiv_vector <- as.matrix(trait_table$adiv)[,1]
adiv_vector
fit <- fastAnc(ctree, adiv_vector, vars=TRUE, CI=TRUE)
fit
