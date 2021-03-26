library(phytools)
library(RColorBrewer)


ctree <- read.tree('../../input/pruned_tree.newick')
preY <- read.table('../../input/trait_table_r11_by_disease_adiv_rank_zeroed.txt', sep="\t", header = TRUE, row.names = 1)
Y <- preY[,c(1:6)] 
disease_color <- c("lightcoral","deepskyblue")
disease_ramp_color_palette <- colorRampPalette(c("lemonchiffon1", "red"))(100)


# All diseases - HICORDIS

### standardized

##### pdf
pdf(file="../../output/disease_heatmap_standardized.pdf")
phylo.heatmap(ctree, Y, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

##### png
png(file="../../output/disease_heatmap_standardized.png")
phylo.heatmap(ctree, Y, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

### nonstandardized

##### pdf
pdf(file="../../output/disease_heatmap_nonstandardized.pdf")
phylo.heatmap(ctree, Y, standardize=FALSE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

##### png
png(file="../../output/disease_heatmap_nonstandardized.png")
phylo.heatmap(ctree, Y, standardize=FALSE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

# SEB and adiv

### pdf
pdf(file="../../output/SEB_vs_adiv_standardized.pdf")
phylo.heatmap(ctree, preY[,c(4,8)], standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

### png
png(file="../../output/SEB_vs_adiv_standardized.png")
phylo.heatmap(ctree, preY[,c(4,8)], standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

# Phylomorphospace - adiv and rank_disase_inverse

ctree <- paintSubTree(ctree, node=38, state="3")
ctree <- paintSubTree(ctree, node=61, state="1")
ctree <- paintSubTree(ctree, node=39, state="2")
cols <- c("deepskyblue", "sandybrown","black")
names(cols)<- 1:3

plot(ctree, cols, fsize=0.5)

### pdf
pdf(file="../../output/adiv_rank_disease_inverse_phylomorphospace.pdf")
phylomorphospace(ctree, preY[,c(8,9)], colors=cols, node.size=1.25, node.by.map=TRUE, fsize=0.01, label="off")
dev.off()

### png
png(file="../../output/adiv_rank_disease_inverse_phylomorphospace.png")
phylomorphospace(ctree, preY[,c(8,9)], colors=cols, node.size=1.25, node.by.map=TRUE, fsize=0.01, label="off")
dev.off()

# PIC

PositivizeContrasts <- function(x, y) {
  #Cache the sign of x so it doesn't
  #change as we reflect points about x-axis!
  sign_of_x <- sign(x)
  x.positivized <- x * sign_of_x
  y.positivized <- y * sign_of_x
  return(cbind(x.positivized, y.positivized))
}

phy <- preY[,c(8,4)]
rawpicX <- pic(phy[,1], ctree)
rawpicY <- pic(phy[,2], ctree)

pos_result <- PositivizeContrasts(rawpicX, rawpicY)

picX <- pos_result[,1]
picY <- pos_result[,2]

pic_fitXY <- lm(picY ~ picX)


sink('../../output/adiv_SEB_pic_results.txt')
print(summary(pic_fitXY))
sink()
