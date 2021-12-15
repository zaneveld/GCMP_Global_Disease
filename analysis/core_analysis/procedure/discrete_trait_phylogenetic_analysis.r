#Output results to a log file as well as to the screen
print("Usage: asr_viz_phytools.r <treefile> <trait_table>")

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
trait_table_fp <-args[1]
tree_fp <- args[2]
x_trait <-args[3]
y_trait <-args[4]

sink(paste("PIC_results_log_",x_trait,"_",y_trait,".txt",sep=""),append=FALSE,split=TRUE)
library(phytools)

print(paste("Analyzing",x_trait,"vs.",y_trait))

trait_table <- read.table(trait_table_fp,comment.char="",header=T,row.names=1,as.is=T,sep="\t")
tree <- read.tree(tree_fp)


#Ensure tree has been filtered so all tips are names in X and Y
filtered_tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% names(Y)])
tree <- drop.tip(filtered_tree,tree$tip.label[!tree$tip.label %in% names(X)])

#Dichotomize tree
tree <- multi2di(tree)

#Plot dominant microbe tree
microbial_trait = 'most_abundant_order_tissue'
pdf(paste("microbial_",microbial_trait,"_reconstruction_rightwards.pdf",sep=""))
#Drop NA rows for our x and y traits
trait_table[trait_table==''] <- NA
trait_table<- trait_table[!is.na(trait_table[,microbial_trait]),]

x <- as.factor(trait_table[,microbial_trait])
names(x) <- row.names(trait_table)
#Drop tree tips without microbiome data
tree <- drop.tip(filtered_tree,tree$tip.label[!tree$tip.label %in% names(x)])

#Dichotomize tree
tree <- multi2di(tree)
setEnv = TRUE
fitER <- rerootingMethod(tree, x, model = "ER")


palette <- c(
  "dodgerblue2",
  "gray70",
  "#FF7F00", # orange
  "#FB9A99", # lt pink
  "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "black", "gold1",
  "skyblue2",
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
cols <- setNames(palette[1:length(unique(x))], sort(unique(x)))
print(cols)

print("simulating trees!")
n_simulations = 100
map.trees<-make.simmap(tree,x,nsim=n_simulations,model="ARD",opt.method="optim")
plotTree(tree,ftype="i",fsize=0.5,offset=0.7,
    lwd=6)
par(fg="transparent",lend=1)
plotTree(tree,ftype="i",fsize=0.5,offset=0.7,
    lwd=4,color="white",add=TRUE)

## now plot our 100 stochastic map trees
## with 99% transparency
for(i in 1:length(map.trees))
    plot(map.trees[[i]],
      colors=sapply(cols,make.transparent,alpha=1.0/n_simulations),
      add=TRUE,lwd=4,ftype="i",fsize=0.5,offset=0.5)

    par(fg="black")
    nodelabels(pie=summary(map.trees)$ace,piecol=cols,cex=0.5)
    legend(x="bottomleft",levels(x),pch=22,
    pt.bg=cols,pt.cex=1.5,bty="n",cex=0.7)


#dev.off()

#plot(map.trees[[1]],colors=cols,lwd=4,ftype="i",fsize=0.5,offset=0.5,add=TRUE)
#nodelabels(pie=summary(map.trees)$ace,piecol=cols,cex=0.5)
#legend(x="bottomleft",levels(x),pch=22,pt.bg=cols,pt.cex=1.5,bty="n",cex=0.7)

dev.off()

## estimate ancestral states under a ER model
#print(lapply(fitER, round, digits = 3))


