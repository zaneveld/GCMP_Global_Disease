
#A script for running phytools PIC and phylomorphospace analysis from
#the command line

#See below for usage. 
print("Usage: Rscript phylomorphospace_r14.r <path_to_trait_table> <path_to_tree> <x_trait_column> <y_trait_column>")

library(ggplot2)
library(phytools)

args <- commandArgs(trailingOnly=TRUE)
trait_table_fp <-args[1]
tree_fp <- args[2]
x_trait <-args[3]
y_trait <-args[4]

#Filtering is not accessible from the commandline, but can be set in the script
#Example of how to filter results to just genera
#with Gammaproteobacteria as the dominant class

#filter_column <- "most_abundant_class_tissue"
#filter_value <- "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria"
filter_column <- NA
filter_value <- NA

filter_column <- args[5]
filter_value <- args[6]

output_dir_for_all_pics = "../output/PIC_results/"
dir.create(output_dir_for_all_pics, showWarnings = FALSE)

output_dir <- paste0(output_dir_for_all_pics,"PIC_",x_trait,"_vs_",y_trait)
if (!is.na(filter_column)){
    output_dir <- paste0(output_dir,"_",filter_column,"_is_",gsub("\\;","_",filter_value))
}
output_dir <- paste0(output_dir,"/")

print(paste("Outputting results to:",output_dir))
dir.create(output_dir)

#Output results to a log file as well as to the screen
sink(paste0(output_dir,x_trait,"_vs_",y_trait,"_results_log.txt"),append=FALSE,split=TRUE)

print("Phylogenetic Independent Contrast Analysis Report")
print("--------------------------------------------------")
print(paste("Analyzing",x_trait,"vs.",y_trait))
print(paste("Trait table filepath:",trait_table_fp))
print(paste("Tree filepath:",tree_fp))
print(paste("Filtering data based on column:",filter_column))
print(paste("Including data only if filter column value is:",filter_value))
print("--------------------------------------------------")

trait_table <- read.table(trait_table_fp,comment.char="",header=T,row.names=1,as.is=T,sep="\t")
tree <- read.tree(tree_fp)
print(tree)
#Filter table to tree tips
trait_table <- trait_table[rownames(trait_table) %in% tree$tip.label,]

if (!is.na(filter_column)){
    trait_table <- subset(trait_table,trait_table[[filter_column]] == filter_value)
}

#Coerce trait columns to numeric, possibly inducing NAs
trait_table[,x_trait]<-as.numeric(trait_table[,x_trait])
trait_table[,y_trait]<-as.numeric(trait_table[,y_trait])

#Drop NA rows for our x and y traits
trait_table<- trait_table[!is.na(trait_table[,x_trait]),]
trait_table <- trait_table[!is.na(trait_table[,y_trait]),]

#Restrict to just our x column of interest
X <- trait_table[,x_trait]
names(X) <- rownames(trait_table)

#Get our y-value column of interest
Y <- trait_table[,y_trait]
names(Y) <- rownames(trait_table)

print(X)
print(Y)

#Ensure tree has been filtered so all tips are names in X and Y
filtered_tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% names(Y)])
tree <- drop.tip(filtered_tree,tree$tip.label[!tree$tip.label %in% names(X)])

#Dichotomize tree
tree <- multi2di(tree)

#Record the filtered tree and trait table

write.tree(tree,paste0(output_dir,x_trait,"_vs_",y_trait,"_filtered_tree.newick"))
write.table(trait_table,paste0(output_dir,x_trait,"_vs_",y_trait,"_filtered_table.csv"))

print("Calculating PICs")
raw.pic.X <- pic(X,tree)
raw.pic.Y <- pic(Y,tree)

#****POSITIVIZE PIC VALUES (negative signs are arbitrary)
#Note: this is not just the absolute value since the difference x=-1 y=2
#has a different interpretation than x=1 y=2

#Code snippet from:https://github.com/bomeara/ComparativeMethodsInR/blob/master/ContinuousTrait_Answers.R
#for contrasts, you should positivize them, since the order doesn't matter. This is NOT taking absolute value.

PositivizeContrasts <- function(x, y) {
    #Cache the sign of x so it doesn't
    #change as we reflect points about x-axis!
    sign_of_x <- sign(x)
    x.positivized <- x * sign_of_x
    y.positivized <- y * sign_of_x
    return(cbind(x.positivized, y.positivized))
}
positivized.results <- PositivizeContrasts(raw.pic.X, raw.pic.Y)
pic.X <-positivized.results[,1]
pic.Y <-positivized.results[,2]

pic_df <- data.frame(pic.X,pic.Y)

# regress through origin (e.g. expect that 0 change in one trait is on average correlated with zero change in the other)

fitYX <- lm(pic.Y ~ pic.X -1 )
print(paste("Summary lm pic.Y ~ pic.X -1 for" ,x_trait,"(x) vs. ",y_trait,"(y)"))
print(summary(fitYX))

## this is a projection of the tree into morphospace
##This code snippit is adapted from a phytools tutorial (http://www.phytools.org/Cordoba2017/ex/3/PICs.html)

pdf(paste0(output_dir,x_trait,"_vs_",y_trait,"_phylomorphospace.pdf"))
phylomorphospace(tree,cbind(X,Y),xlab=x_trait,ylab=y_trait,label="off",node.size=c(0,0))
points(X,Y,pch=21,bg="firebrick",cex=1.4)
dev.off()

pdf(paste0(output_dir,x_trait,"_vs_",y_trait,"_rank_phylomorphospace.pdf",sep=""))
phylomorphospace(tree,cbind(rank(X),rank(Y)),xlab=paste("rank ",x_trait),ylab=paste("rank ",y_trait),label="off",node.size=c(0,0))
points(rank(X),rank(Y),pch=21,bg= "firebrick",cex=1.4)
dev.off()


#Save raw PIC contrasts as a pdf
pdf(paste0(output_dir,x_trait,"_vs_",y_trait,"_pic_scatter_YX.pdf"))

print("Dataframe for PIC analysis:")
print(pic_df)

#Do some quick calculations to figure out the y-axis limits
#We want them symmetrical above and below the null 0 slope expectation

largest_limit = max(abs(min(pic.Y)),abs(max(pic.Y)))

#largest limit must be positive, so set the lower bound
#to be that far below 0
lower_bound = largest_limit * -1
upper_bound = largest_limit 

#Set the null model label 75% of the way along the axis
null_model_label_x_pos = max(pic.X)*0.75 
null_model_label_y_pos = 0 + 0.05*upper_bound

fontsize = 16

ggplot(pic_df, aes(pic.X,pic.Y)) + 
    geom_smooth(method = "lm",se=TRUE,level=0.95,col = "black",fill="cyan",formula = y ~ x -1) +
    geom_point(size = 3, shape=21, fill = "lightgrey",colour = "black",alpha=0.5) + 
    labs(x = paste("Contrast in ",x_trait), y = paste("Contrast in ",y_trait))+
    scale_y_continuous(limits=c(lower_bound,upper_bound))+ 
    geom_hline(yintercept=0, color = "firebrick", size=1, linetype="dotted")+
    geom_label(
        label="Null expectation", 
        x=null_model_label_x_pos,
        y=null_model_label_y_pos,
        label.size = NA,
        size = 6, #not a font size, but a point size
        color = "firebrick",
        face = "bold",
        fill = NA
    )+
    
    theme_classic()+
    theme(axis.text=element_text(size=16))+
    theme(axis.title=element_text(size=16))+
    theme(text = element_text(family = "sans"))

dev.off()

#Build contmap for trait X
fit<-fastAnc(tree,X,vars=TRUE,CI=TRUE)
#Print model fit to screen
print(paste(c("FastAnc ML modelfit for",x_trait)))
print(fit)
obj <- contMap(tree,X,plot=F)
tree_direction <- "rightwards"
inverse_green_colorscheme <-c('black','springgreen3','yellow','white')
obj <- setMap(obj,colors=inverse_green_colorscheme)

#Write contmap for trait X to file
pdf(paste0(output_dir,x_trait,"_asr_contmap.pdf"))
par(mai=c(12.12,1,1.1,1.1))
plot(obj,direction=tree_direction,legend=0.7*max(nodeHeights(tree)),fsize=c(0.222,0.9))
axis(1)
title(xlab="time from the root (mya)")
dev.off()

#Build contmap for trait y
fit<-fastAnc(tree,Y,vars=TRUE,CI=TRUE)
#Print model fit to screen
print(paste(c("FastAnc ML modelfit for",y_trait)))
print(fit)
obj <- contMap(tree,Y,plot=F)
tree_direction <- "leftwards"
inverse_green_colorscheme <-c('black','springgreen3','yellow','white')
obj <- setMap(obj,colors=inverse_green_colorscheme)

#Write contmap for trait X to file
pdf(paste0(output_dir,y_trait,"_asr_contmap_leftwards.pdf"))
par(mai=c(12.12,1,1.1,1.1))
plot(obj,direction=tree_direction,legend=0.7*max(nodeHeights(tree)),fsize=c(0.222,0.9))
axis(1)
title(xlab="time from the root (mya)")
dev.off()

