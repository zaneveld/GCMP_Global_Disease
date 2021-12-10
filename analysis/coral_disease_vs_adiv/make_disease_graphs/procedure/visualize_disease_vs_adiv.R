library(phytools)
library(RColorBrewer)


ctree <- read.tree('../../input/pruned_tree.newick')
trait_table <- read.table('../../input/trait_table_r11_by_disease_adiv_rank_zeroed.txt', sep="\t", header = TRUE, row.names = 1)
disease_table <- trait_table[,c(1:6)] 
disease_color <- c("lightcoral","deepskyblue")
disease_ramp_color_palette <- colorRampPalette(c("lemonchiffon1", "red"))(100)


# All diseases - HICORDIS

### standardized

##### pdf
pdf(file="../../output/disease_heatmap_standardized.pdf")
phylo.heatmap(ctree, disease_table, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

##### png
par(mar=c(1,1,1,1))
png(file="../../output/disease_heatmap_standardized.png", res=300)
phylo.heatmap(ctree, disease_table, standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

### nonstandardized

##### pdf
pdf(file="../../output/disease_heatmap_nonstandardized.pdf")
phylo.heatmap(ctree, disease_table, standardize=FALSE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

##### png
png(file="../../output/disease_heatmap_nonstandardized.png", res=300)
phylo.heatmap(ctree, disease_table, standardize=FALSE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

# SEB and adiv

### pdf
pdf(file="../../output/SEB_vs_adiv_standardized.pdf")
phylo.heatmap(ctree, trait_table[,c(4,8)], standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
dev.off()

### png
png(file="../../output/SEB_vs_adiv_standardized.png", res=300)
phylo.heatmap(ctree, trait_table[,c(4,8)], standardize=TRUE, fsize=0.6, colors=disease_ramp_color_palette,length=1)
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
phylomorphospace(ctree, trait_table[,c(trait_table$adiv,trait_table$rank_disease_inverse)], colors=cols, node.size=1.25, node.by.map=TRUE, fsize=0.01, label="off")
dev.off()

### png
png(file="../../output/adiv_rank_disease_inverse_phylomorphospace.png", res=300)
phylomorphospace(ctree, trait_table[,c(trait_table$adiv,trait_table$rank_disease_inverse)], colors=cols, node.size=1.25, node.by.map=TRUE, fsize=0.01, label="off")
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

phy <- trait_table[,c(8,4)]
rawpicX <- pic(phy[,1], ctree)
rawpicY <- pic(phy[,2], ctree)

pos_result <- PositivizeContrasts(rawpicX, rawpicY)

picX <- pos_result[,1]
picY <- pos_result[,2]

pic_fitXY <- lm(picY ~ picX - 1)


sink('../../output/adiv_SEB_pic_results.txt')
print(summary(pic_fitXY))
sink()

# SEB adiv PIC scatterplot

### pdf
pdf(file="../../output/adiv_SEB_pic_scatterplot.pdf")
plot(picX,picY,xlab="Microbiome alpha diversity",ylab="Skeletal Eroding Band",bg='gray',pch=16)
abline(pic_fitXY, col="red")
dev.off()

# FIX THIIIISSSS
### png
png(file="../../output/adiv_SEB_pic_scatterplot.png", res=300)
par(mar=c(1,1,1,1))
plot(picX,picY,xlab="Microbiome alpha diversity",ylab="Skeletal Eroding Band",bg='gray',pch=16)
abline(pic_fitXY, col="red")
dev.off()

# Total Disease barplot

#pdf
pdf(file="../../output/total_disease_barplot.pdf")
plotTree.barplot(ctree, trait_table[,c(7,11)],fsize=0.5, col=disease_color,xlab="Total Disease")
dev.off()