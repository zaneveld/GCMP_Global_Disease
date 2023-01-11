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

#create this file only the first time then we append to it


print("Usage: Rscript pgl_phylto_taxa_gcmp_no_loop_v2.r <path_to_trait_table> <path_to_tree> <x_trait_column> <y_trait_column> <filter_column> <filter_value> <output_dir_suffix>")

#Get the user input and assign variables
args <- commandArgs(trailingOnly = TRUE)
trait_table_fp <- args[1]
tree_fp <- args[2]
x_trait <- args[3]
y_trait <- args[4]
filter_column <- args[5]
filter_value <- args[6]
user_specified_output_dir <- args[7]

#Set output directory for all the PGLS models
#results of this analysis will be created in a subdirectory
output_dir_for_all_pgls = "../output/PGLS_results"
if (!is.na(user_specified_output_dir) & (user_specified_output_dir != 'None')){
  output_dir_for_all_pgls = paste0(args[7],"/")
}
dir.create(output_dir_for_all_pgls, showWarnings = FALSE)
#Allow users to add a special suffix to the output

suffix <- args[8]

output_dir <- paste0(output_dir_for_all_pgls,"PGLS_",x_trait,"_vs_",y_trait)
pgl.output = data.frame("Disease_Trait", "Taxa_string", "N_Unique_Samples", "Taxa", "N_Microbes", "Compartment", "Package", "Model", "pVal", "R_Squared", "Adj_R_Squared", "x_Trait_Slope_95CI", "Intercept_95CI", "Parameters", "Estimated_Parameter_95CI", "AIC", "AICc", "Minimum_AIC", "Minimum_AICc")
pgl.output
write.table(output_dir,pgl.output,file="PGLS_restults.csv",append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#Add the filter column and value to directory name, if they were provided
if (!is.na(filter_column) & (filter_column != 'None') & (filter_value != 'None')){
  output_dir <- paste0(output_dir,"_",filter_column,"_is_",gsub("\\;","_",filter_value))
}

#Add a special suffix to the output dir, if user requested one
if (!is.na(suffix)){
  output_dir <- paste0(output_dir,"_",suffix)
}

print(paste("Output results to:",output_dir))
dir.create(output_dir,showWarning = FALSE)

#Output results to a log file and screen
sink(paste0(output_dir,x_trait,"_vs_",y_trait,"_results_log.txt"))

print("PGLS Analysis Report")
print (paste("Analyzing", x_trait, "vs.", y_trait))
print(paste("Trait table filepath:",trait_table_fp))
print(paste("Tree filepath:",tree_fp))
print(paste("Filtering data based on column:",filter_column))
print(paste("Including data only if filter column value is:", filter_value))


#filter column and the value to filter by within the column
filter_column <- filter_column
filter_value <- filter_value
filter_value_last <-sub('.*\\;', '', filter_value)


#Pick the columns to analyze from the various tissue compartments (mucus, skeleton, tissue, all)
#x traits can be: will be the bacteria of interest (copied strait from trait table)
x_trait_column <-x_trait
#only extract the last part of the name
filter_taxa_name <-str_extract(x_trait_column,"[^_]+$")

#y trait disease can be: sum_dis, sum_healthy, sum_bbd, sum_wb, perc_healthy, perc_dis, perc_bbd, perc_wd
y_trait_column<-y_trait

trait_table <- read.csv(trait_table_fp,header=TRUE, row.names=1)

#filter trait table by the column and values listed (add a ! before trait_table after the comma)
#if (!is.na(filter_column)){
  #trait_table <- subset(trait_table,!(trait_table[[filter_column]] == filter_value))
#}

if (!is.na(filter_column)){
  trait_table <- subset(trait_table,trait_table[[filter_column]] == filter_value)
}



#extract the compartment to be used, number of microbes, package

compartment <-str_extract(x_trait_column,"^[^_]+(?=_)")
compartment
n_microbe <-sum(trait_table$sum_total,na.rm=TRUE)
#n_microbe <- sum(trait_table.df$sample_column,na.rm=TRUE)
n_microbe

package <- "caper"

#sample_number <-trait_table[[sample_column]]
#reduce trait table to just columns
trait_table_small <- trait_table[,c(x_trait_column,y_trait_column)]
trait_table_small <- na.omit(trait_table_small)
x_trait <-trait_table_small[[x_trait_column]]
y_trait <-trait_table_small[[y_trait_column]]

sample_column <- unique(x_trait)
sample_column 
n_unique_sample <- length(sample_column)
n_unique_sample

#load the newick created previously at the genus level.
tree <- read.tree(tree_fp)


#make sure that data and tree use same species names
coral_obj <- name.check(tree,trait_table_small)
coral_obj

#take out tree tips that do not have data
#we have data not in tree for some species. These do not have disease data
tree.prune<-drop.tip(tree, coral_obj$tree_not_data)
tree.prune_obj <-name.check(tree.prune,trait_table_small)
tree.prune_obj

#delete the extra rows from the object
trait_table.prune<-trait_table_small[-which(rownames(trait_table_small)%in% tree.prune_obj$data_not_tree),]

#commented out parts that are not needed when you don't need to prune the table
prune_tree_data <- name.check(tree.prune,trait_table.prune)
prune_tree_data
plot(tree.prune)

trait_table.prune$host_genus<-rownames(trait_table.prune)

#trait_table_small$host_genus<-rownames(trait_table_small)

#analyze data using the package caper which combines the data and tree into one file

comp.data <- comparative.data(tree.prune, trait_table.prune, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)
#comp.data <- comparative.data(tree.prune, trait_table_small, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)
model.formula <- as.formula(paste0(y_trait_column,'~',x_trait_column))

#run data using the basic model
model.3_formula<-"pgls(y_trait~x_trait, data=comp.data)"
model.3_param <-"lambda=1,delta=1,kappa=1"
model.3<-pgls(model.formula, data=comp.data,lambda=1,kappa=1,delta=1)
summary(model.3)

#extract the AIC and AICc from Model 3
model.3_aic <- AIC(model.3)
model.3_aic
model.3_aicc <- AICc(model.3)
model.3_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.3<-coef(summary(model.3))
#get only the p value from the summary table for model 3
model.3_pval<- (summary_model.3)[2,4]
#extract the intercept, its standard deviation, and 95% CI.
model.3_intercept<-(summary_model.3)[1,1]
model.3_intercept.sd<-(summary_model.3)[1,2]
model.3_intercept.95CI.pos<-(model.3_intercept.sd+1.96)
model.3_intercept.95CI.neg<-(model.3_intercept.sd-1.96)
model.3_intercept_output<-paste(model.3_intercept, "+", model.3_intercept.95CI.pos, "-", model.3_intercept.95CI.neg)
#extract the slope, its standard deviation, and its 95% CI.
model.3_slope<-(summary_model.3)[2,1]
model.3_slope.sd<-(summary_model.3)[2,2]
model.3_slope.95CI.pos<-(model.3_slope.sd+1.96)
model.3_slope.95CI.neg<-(model.3_slope.sd-1.96)
model.3_slope_output<-paste(model.3_slope, "+", model.3_slope.95CI.pos, "-", model.3_slope.95CI.neg)
#extract the R squared and adjusted r squared values for model 3
model.3_rSquared <- summary(model.3)$r.squared
model.3_rSquared.adj<-(summary(model.3)$adj.r.squared)[1,1]
model.3_rSquared.adj
#estimated parameters (none for this model)
estimated.param <- "NA"
print(paste ("The pvalue for model 3 is ", model.3_pval, "R squared is ", model.3_rSquared, "Adjust R squared is", model.3_rSquared.adj, "and the AIC/AICc is ", model.3_aic, model.3_aicc))
output_row.model.3 <- data.frame(y_trait_column, x_trait_column, filter_column, filter_value, n_microbe, compartment, package, model.3_formula, model.3_pval, model.3_rSquared, model.3_rSquared.adj,model.3_slope_output, model.3_intercept_output,model.3_param,estimated.param, model.3_aic, model.3_aicc)
print(output_row.model.3)

#Write it to the csv file
write.table(output_row.model.3, file=pgl.output, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#run data setting lambda to ML
model.3l_formula<-"pgls(y_trait~x_trait, data=comp.data, lambda=ML,kappa=1,delta=1)"
model.3l_param <-"lambda=ML,delta=1,kappa=1"
model.3l<-pgls(model.formula, data=comp.data,lambda="ML",kappa=1,delta=1)
summary(model.3l)

#extract the AIC and AICc from Model 3
model.3l_aic <- AIC(model.3l)
model.3l_aic
model.3l_aicc <- AICc(model.3l)
model.3l_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.3l<-coef(summary(model.3l))
#get only the p value from the summary table for model 3
model.3l_pval<- (summary_model.3l)[2,4]
#extract the intercept and its standard deviation.
model.3l_intercept<-(summary_model.3l)[1,1]
model.3l_intercept.sd<-(summary_model.3l)[1,2]
model.3l_intercept.95CI.pos<-(model.3l_intercept.sd+1.96)
model.3l_intercept.95CI.neg<-(model.3l_intercept.sd-1.96)
model.3l_intercept_output<-paste(model.3l_intercept, "+", model.3l_intercept.95CI.pos, "-", model.3l_intercept.95CI.neg)
#extract the slope and its standard deviation.
model.3l_slope<-(summary_model.3l)[2,1]
model.3l_slope.sd<-(summary_model.3l)[2,2]
model.3l_slope.95CI.pos<-(model.3l_slope.sd+1.96)
model.3l_slope.95CI.neg<-(model.3l_slope.sd-1.96)
model.3l_slope_output<-paste(model.3l_slope, "+", model.3l_slope.95CI.pos, "-", model.3l_slope.95CI.neg)
#extract the R squared and adjusted r squared values for model 3
model.3l_rSquared <- summary(model.3l)$r.squared
model.3l_rSquared.adj<-(summary(model.3l)$adj.r.squared)[1,1]
#estimated parameter
estimated.param <- "lambda"
#extract the branch length parameter
lambda<-model.3l$param["lambda"]
lambda
lambda.CI<-pgls.confint(model.3l,"lambda")$ci.val
lambda.CI
lambda_output<-paste(lambda, "-/+", lambda.CI)

print(paste ("The pvalue for model 3l is ", model.3l_pval, "R squared is ", model.3l_rSquared, "Adjust R squared is", model.3l_rSquared.adj, "the AIC/AICc is ", model.3l_aic, model.3_aicc, "and the estimated paramter value of lambda is", lambda))

#create to output row for the data that allows lambda to be chosen.
output_row.model.3l <- data.frame(y_trait_column, x_trait_column, filter_column, filter_value, n_microbe, compartment, package, model.3l_formula, model.3l_pval, model.3l_rSquared, model.3l_rSquared.adj, model.3l_slope_output, model.3l_intercept_output, model.3l_param, lambda_output, model.3l_aic, model.3l_aicc)
print(output_row.model.3l)

#Write it to the csv file
write.table(output_row.model.3l, file=pgl.output, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)

#run data setting delta to ML
model.3d_formula<-"pgls(y_trait~x_trait, data=comp.data, lambda=1,kappa=1,delta=ML)"
model.3d_param <-"lambda=l,delta=ML,kappa=1"
model.3d<-pgls(model.formula, data=comp.data,lambda=1,kappa=1,delta="ML")
summary(model.3d)
#extract the AIC and AICc from Model 3d
model.3d_aic <- AIC(model.3d)
model.3d_aic
model.3d_aicc <- AICc(model.3d)
model.3d_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.3d<-coef(summary(model.3d))
#get only the p value from the summary table for model 3
model.3d_pval<- (summary_model.3d)[2,4]
#extract the intercept and its standard deviation.
model.3d_intercept<-(summary_model.3d)[1,1]
model.3d_intercept.sd<-(summary_model.3d)[1,2]
model.3d_intercept.95CI.pos<-(model.3d_intercept.sd+1.96)
model.3d_intercept.95CI.neg<-(model.3d_intercept.sd-1.96)
model.3d_intercept_output<-paste(model.3d_intercept, "+", model.3d_intercept.95CI.pos, "-", model.3d_intercept.95CI.neg)
#extract the slope and its standard deviation.
model.3d_slope<-(summary_model.3d)[2,1]
model.3d_slope.sd<-(summary_model.3d)[2,2]
model.3d_slope.95CI.pos<-(model.3d_slope.sd+1.96)
model.3d_slope.95CI.neg<-(model.3d_slope.sd-1.96)
model.3d_slope_output<-paste(model.3d_slope, "+", model.3d_slope.95CI.pos, "-", model.3d_slope.95CI.neg)
#extract the R squared and adjusted r squared values for model 3
model.3d_rSquared <- summary(model.3d)$r.squared
model.3d_rSquared.adj<-(summary(model.3d)$adj.r.squared)[1,1]
#estimated parameter
esimated.param <- "delta"
delta<-model.3d$param["delta"]
delta.CI<-pgls.confint(model.3d,"delta")$ci.val
delta.CI
delta_output<-paste(delta, "-/+", delta.CI)
print(paste ("The pvalue for model 3d is ", model.3d_pval, "R squared is ", model.3d_rSquared, "Adjust R squared is", model.3d_rSquared.adj, "the AIC/AICc is ", model.3d_aic, model.3d_aicc, "and the estimated parameter of delta is", delta))

#create the output row when delta is not fixed.
output_row.model.3d <- data.frame(y_trait_column, x_trait_column, filter_column, filter_value, n_microbe, compartment, package, model.3d_formula, model.3d_pval, model.3d_rSquared, model.3d_rSquared.adj,model.3d_slope_output, model.3d_intercept_output, model.3d_param, delta_output, model.3d_aic, model.3d_aicc)
print(output_row.model.3d)

#Write it to the csv file
write.table(output_row.model.3d, file=pgl.output, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)

#run data setting kappa to ML
model.3k_formula<-"pgls(y_trait~x_trait, data=comp.data, lambda=1,kappa=ML,delta=1)"
model.3k_param <-"lambda=1,delta=1,kappa=ML"
model.3k<-pgls(model.formula, data=comp.data,lambda=1,kappa="ML",delta=1)
summary(model.3k)
#extract the AIC and AICc from Model 3k
model.3k_aic <- AIC(model.3k)
model.3k_aic
model.3k_aicc <- AICc(model.3k)
model.3k_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.3k<-coef(summary(model.3k))
#get only the p value from the summary table for model 3k
model.3k_pval<- (summary_model.3k)[2,4]
#extract the intercept and its standard deviation.
model.3k_intercept<-(summary_model.3k)[1,1]
model.3k_intercept.sd<-(summary_model.3k)[1,2]
model.3k_intercept.95CI.pos<-(model.3k_intercept.sd+1.96)
model.3k_intercept.95CI.neg<-(model.3k_intercept.sd-1.96)
model.3k_intercept_output<-paste(model.3k_intercept, "+", model.3k_intercept.95CI.pos, "-", model.3k_intercept.95CI.neg)
#extract the slope and its standard deviation.
model.3k_slope<-(summary_model.3k)[2,1]
model.3k_slope.sd<-(summary_model.3k)[2,2]
model.3k_slope.95CI.pos<-(model.3k_slope.sd+1.96)
model.3k_slope.95CI.neg<-(model.3k_slope.sd-1.96)
model.3k_slope_output<-paste(model.3k_slope, "+", model.3k_slope.95CI.pos, "-", model.3k_slope.95CI.neg)

#extract the R squared and adjusted r squared values for model 3
model.3k_rSquared <- summary(model.3k)$r.squared
model.3k_rSquared.adj<-(summary(model.3k)$adj.r.squared)[1,1]
#estimated parameter
esimated.param <- "kappa"
kappa<-model.3k$param["kappa"]
kappa.CI<-pgls.confint(model.3k,"kappa")$ci.val
kappa.CI
kappa_output<-paste(kappa, "-/+", kappa.CI)
print(paste ("The pvalue for model 3k is ", model.3k_pval, "R squared is ", model.3k_rSquared, "the AIC/AICc is ", model.3k_aic, model.3k_aicc, "and the estimated parameter for kappa is", kappa))


#write a data frame with each AIC for each combination and identify the AIC and AICc with 
#the lowest value.
aic.df <- data.frame(Model_3_AIC=model.3_aic, Model_3l_AIC=model.3l_aic, Model_3d_AIC=model.3d_aic, Model_3k_AIC=model.3k_aic)
aic.df
aic.df_columns <- colnames(aic.df)
aic.df_columns
aic.df$min <- apply(aic.df,1,which.min)
aic.df$min
aic.df$min_model <- colnames(aic.df)[aic.df$min]
aic.df$min_model

aicc.df <- data.frame(Model_3_AICc=model.3_aicc, Model_3l_AICc=model.3l_aicc, Model_3d_AICc=model.3d_aicc, Model_3k_AICc=model.3k_aicc)
aicc.df
aicc.df_columns <- colnames(aicc.df)
aicc.df_columns
aicc.df$min <- apply(aicc.df,1,which.min)
aicc.df$min
aicc.df$min_model <- colnames(aicc.df)[aicc.df$min]
aicc.df$min_model

#create the output row for the data when kappa is not fixed.
output_row.model.3k <- data.frame(y_trait_column, x_trait_column, filter_column, filter_value, n_microbe, compartment, package, model.3k_formula, model.3k_pval, model.3k_rSquared, model.3k_rSquared.adj, model.3k_slope_output, model.3k_intercept_output, model.3k_param, kappa_output, model.3k_aic, model.3k_aicc, aic.df$min_model, aicc.df$min_model)
print(output_row.model.3k)

#Write it to the csv file
write.table(output_row.model.3k, file=pgl.output, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


