#Phylogenetic generalized least squares regression and phylogenetic generalized ANOVA
#load libraries
library(ape)
library(nlme)
library(geiger)
library(caper)
library(MuMIn)
library(rr2)
library(phytools)
library(stringr)
library(MASS)
library(robustbase)
library(ggplot2)
library(remotes)
library(dplyr)

#create this file only the first time then we append to it
pgl.output = data.frame("Disease_Trait", "Taxa_String", "N_Unique_Sample", "Taxa", "N_Microbes", "Compartment", "Package", "Model", "pVal", "R_Squared", "Adj_R_Squared", "x_Trait_Slope_95CI", "Intercept_95CI", "Parameters", "Estimated_Parameter_95CI", "AIC", "AICc", "Minimum_AIC", "Minimum_AICc")
pgl.output
write.table(pgl.output,file="../output/alpha_diversity_taxa/GCMP_dominant_taxa_vs_disease_susceptibility_compartments_pgls.csv",append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#input the trait table
trait.table.input <-"../input/GCMP_trait_table_with_abundances_and_adiv_and_metadata_and_growth_data_pcoa.csv"
trait.table <- read.csv(trait.table.input, header=TRUE, row.names = 1)
trait.table.df <- data.frame(trait.table)

#extract the number of microbes and package used
n.microbe <-"463791"
package <- "caper"

#Select only data from the Percent disease column and columns that have genus level microbe data.
small.data.test <- trait.table.df %>% select('perc_dis', contains("___D_5__"))

#convert the index column to 1st column and then convert it back to index column
small.data.test <-cbind(host_genus=rownames(small.data.test),small.data.test)
rownames(small.data.test) <-1:nrow(small.data.test)
rownames(small.data.test) <- small.data.test$host_genus

#import the phylogenetic tree
tree <- "../input/huang_roy_genus_tree.newick"
tree <- read.tree(tree)

#get y trait column
y.trait.column = small.data.test$perc_dis

#run each of the columns with genus information through the loop and compare it to the perc_dis column.
for (column.name in colnames(small.data.test)){
  if (column.name == "perc_dis" || column.name =="host_genus"){
    print("Skipping perc_dis or host_genus!")
    next
  }
  
  x.trait.column <- small.data.test[,column.name]
  
  #Extract the compartment from the x trait column name
  compartment <-str_extract(column.name,"^[^_]+(?=_)")
  
  #Extract only the final portion of the taxa name. In this case the genus
  filter.taxa.name <-str_extract(column.name,"[^_]+$")
  
  
  small.df <- data.frame(small.data.test$host_genus, x.trait.column, y.trait.column)
  
  #rename the column name to just host_genus
  names(small.df)[names(small.df)=='small.data.test.host_genus']<-'host_genus'
  
  #move host_genus as an index column as well
  rownames(small.df) <- small.df$host_genus
  
  #drop rows without data
  small.filtered.df <- na.omit(small.df)
  
  #make sure that data and tree use same species names
  #merge into a single tree object
  tree.obj <- name.check(tree,small.filtered.df)
  
  #take out tree tips that do not have data
  #we have data not in tree for some species. These do not have disease data
  pruned.tree<-drop.tip(tree, tree.obj$tree_not_data)
  pruned.tree.obj <-name.check(pruned.tree,small.filtered.df)
  
  
  #analyze data using the package caper which combines the data and tree into one file
  x.data <- small.filtered.df$x.trait.column
  y.data <- small.filtered.df$y.trait.column
  
  comp.data <- comparative.data(pruned.tree, small.filtered.df, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)

  #skip any data that only has 1 unique value since PGLS can not be run on this.
  #This is for unique x trait values
  unique.x.trait.values <- unique(x.data)
  n.unique.x.trait.values <- length(unique.x.trait.values)
  print(unique.x.trait.values)
  print(n.unique.x.trait.values)
  unique.val <- print(paste(n.unique.x.trait.values))
  if (n.unique.x.trait.values <= 2){
    print("This column doesn't have enough x values. Skipping")
    next
  }
  #This is for unique y trait values
  unique.y.trait.values <- unique(y.data)
  n.unique.y.trait.values <- length(unique.y.trait.values)
  print(unique.y.trait.values)
  print(n.unique.y.trait.values)
  if (n.unique.y.trait.values <= 2){
    print("This column doesn't have enough y values. Skipping")
    next
  }
  
  #run data using the basic model
  model.formula <- as.formula("y.trait.column ~ x.trait.column")
  print(model.formula)
  model.3<-pgls(model.formula, data=comp.data)
  model.3.summary <- summary(model.3)
  #prepare the formula and paramters for the output
  model.3.formula<-"pgls(x_trait~y_trait, data=comp.data)"
  model.3.param <-"lambda=1,delta=1,kappa=1"
  estimate.param <-"NA"
  #print(model.3.summary)
  
  #function for the output based on the model you are running
  summary.output <- function(model){
    aic <- AIC(model)
    aicc <- AICc(model)
    #extract only the coefficient output so that we can more easily extract the p value
    summary.model<-coef(summary(model))
    #get only the p value from the summary table for the model
    model.pval<- (summary.model)[2,4]
    #extract the intercept, its standard deviation, and 95% CI.
    model.intercept<-(summary.model)[1,1]
    model.intercept.sd<-(summary.model)[1,2]
    model.intercept.95CI.pos<-(model.intercept.sd+1.96)
    model.intercept.95CI.neg<-(model.intercept.sd-1.96)
    model.intercept.output<-paste(model.intercept, "+", model.intercept.95CI.pos, "-", model.intercept.95CI.neg)
    #extract the slope, its standard deviation, and its 95% CI.
    model.slope<-(summary.model)[2,1]
    model.slope.sd<-(summary.model)[2,2]
    model.slope.95CI.pos<-(model.slope.sd+1.96)
    model.slope.95CI.neg<-(model.slope.sd-1.96)
    model.slope.output<-paste(model.slope, "+", model.slope.95CI.pos, "-", model.slope.95CI.neg)
    #extract the R squared and adjusted r squared values for model 3
    model.rSquared <- summary(model)$r.squared
    model.rSquared.adj<-(summary(model.3)$adj.r.squared)[1,1]
    model.rSquared.adj
    
    return(list(model.pval=model.pval, model.rSquared=model.rSquared, model.rSquared.adj=model.rSquared.adj, model.slope.output= model.slope.output, model.intercept.output=model.intercept.output, aic=aic, aicc=aicc))
  }
  model.output<-summary.output(model.3)
  #print(model.output)
  print(paste ("The pvalue for model 3 is ", model.output$model.pval, "R squared is ", model.output$model.rSquared, "Adjust R squared is", model.output$model.rSquared.adj, "and the AIC/AICc is ", model.output$aic, model.output$aicc))
  
  #create a output row for the output table
  
  output.row.model.3 <- data.frame("perc_dis", column.name, unique.val,
    filter.taxa.name, n.microbe, compartment, package, model.3.formula, 
    model.output$model.pval, model.output$model.rSquared,
    model.output$model.rSquared.adj, model.output$model.slope.output, 
    model.output$model.intercept.output, model.3.param, estimate.param, 
    model.output$aic, model.output$aicc)
  
  #Write it to the csv file
  write.table(output.row.model.3, file="../output/alpha_diversity_taxa/GCMP_dominant_taxa_vs_disease_susceptibility_compartments_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
  
  #run next model that has lambda as ML
  model.3l.formula<-"pgls(x_trait~y_trait, data=comp.data, lambda=ML,kappa=1,delta=1)"
  model.3l.param <-"lambda=ML,delta=1,kappa=1"
  model.3l<-pgls(model.formula, data=comp.data,lambda="ML",kappa=1,delta=1)
  model.3l.summary <- summary(model.3l)
  #print(model.3l.summary)
  
  model.output<-summary.output(model.3l)
  #print(model.output)
  print(paste ("The pvalue for model 3l is ", model.output$model.pval, "R squared is ", model.output$model.rSquared, "Adjust R squared is", model.output$model.rSquared.adj, "and the AIC/AICc is ", model.output$aic, model.output$aicc))
  
  #extract the branch length parameter
  lambda<-model.3l$param["lambda"]
  lambda
  lambda.CI<-pgls.confint(model.3l,"lambda")$ci.val
  lambda.CI
  lambda.output<-paste(lambda, "-/+", lambda.CI)
  
  #create a output row for the output table
  print(paste("Model slope output:",model.output$model.slope.output))
  output.row.model.3l <- data.frame("perc_dis", column.name,
      n.unique.x.trait.values, filter.taxa.name, n.microbe, 
      compartment, package, model.3l.formula, model.output$model.pval, 
      model.output$model.rSquared, model.output$model.rSquared.adj, 
      model.output$model.slope.output, model.output$model.intercept.output,
      model.3l.param, lambda.output, model.output$aic, model.output$aicc)
  #print(output.row.model.3l)
  
  #Write the row to the csv file
  write.table(output.row.model.3l, file="../output/alpha_diversity_taxa/GCMP_dominant_taxa_vs_disease_susceptibility_compartments_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
  
  #run data setting delta to ML
  model.3d.formula<-"pgls(x_trait~y_trait, data=comp.data, lambda=1,kappa=1,delta=ML)"
  model.3d.param <-"lambda=l,delta=ML,kappa=1"
  model.3d<-pgls(model.formula, data=comp.data,lambda=1,kappa=1,delta="ML")
  summary(model.3d)
  
  model.output<-summary.output(model.3d)
  #print(model.output)
  print(paste ("The pvalue for model 3d is ", model.output$model.pval, "R squared is ", model.output$model.rSquared, "Adjust R squared is", model.output$model.rSquared.adj, "and the AIC/AICc is ", model.output$aic, model.output$aicc))
  
  #extract the branch length parameter
  delta<-model.3d$param["delta"]
  delta.CI<-pgls.confint(model.3d,"delta")$ci.val
  delta.CI
  delta.output<-paste(delta, "-/+", delta.CI)
  
  #create a output row for the output table
  output.row.model.3d <- data.frame("perc_dis", column.name, 
      n.unique.x.trait.values, filter.taxa.name, n.microbe, compartment, 
      package, model.3d.formula, model.output$model.pval, 
      model.output$model.rSquared, model.output$model.rSquared.adj, 
      model.output$model.slope.output, model.output$model.intercept.output, 
      model.3d.param, delta.output, model.output$aic, model.output$aicc)
  #print(output.row.model.3d)
  
  #Write the row to the csv file
  write.table(output.row.model.3d, file="../output/alpha_diversity_taxa/GCMP_dominant_taxa_vs_disease_susceptibility_compartments_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
  
  #run model setting kappa to ML
  model.3k.formula<-"pgls(x_trait~y_trait, data=comp.data, lambda=1,kappa=ML,delta=1)"
  model.3k.param <-"lambda=1,delta=1,kappa=ML"
  model.3k<-pgls(model.formula, data=comp.data,lambda=1,kappa="ML",delta=1)
  summary(model.3k)
  
  model.output<-summary.output(model.3k)
  #print(model.output)
  print(paste ("The pvalue for model 3k is ", model.output$model.pval, "R squared is ", model.output$model.rSquared, "Adjust R squared is", model.output$model.rSquared.adj, "and the AIC/AICc is ", model.output$model.aic, model.output$model.aicc))
  
  #extract the branch length parameter
  kappa<-model.3k$param["kappa"]
  kappa.CI<-pgls.confint(model.3k,"kappa")$ci.val
  kappa.CI
  kappa.output<-paste(kappa, "-/+", kappa.CI)
  
  #create a output row for the output table
  output.row.model.3k <- data.frame("perc_dis", column.name, 
      n.unique.x.trait.values, filter.taxa.name, n.microbe, compartment, 
      package, model.3k.formula, model.output$model.pval, 
      model.output$model.rSquared, model.output$model.rSquared.adj, 
      model.output$model.slope.output, model.output$model.intercept.output, 
      model.3k.param, kappa.output, model.output$aic, model.output$aicc)
  #print(output.row.model.3k)
  
  
  #Write the row to the csv file
  write.table(output.row.model.3k, file="../output/alpha_diversity_taxa/GCMP_dominant_taxa_vs_disease_susceptibility_compartments_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
  

}

