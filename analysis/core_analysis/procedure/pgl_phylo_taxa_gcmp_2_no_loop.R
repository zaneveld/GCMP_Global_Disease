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
pgl.output = data.frame("Disease_Trait", "Taxa_string", "N_Unique_Samples", "Taxa", "N_Microbes", "Compartment", "Package", "Model", "pVal", "R_Squared", "Adj_R_Squared", "x_Trait_Slope_95CI", "Intercept_95CI", "Parameters", "Estimated_Parameter_95CI", "AIC", "AICc", "Minimum_AIC", "Minimum_AICc")
pgl.output
write.table(pgl.output,file="../output/alpha_diversity_taxa/GCMP_taxa_specific_alpha_compartments_pgls.csv",append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#filter column and the value to filter by within the column
filter_column <- "functional_group_sensu_darling"
filter_value <- "Generalist"
filter_value_last <-sub('.*\\;', '', filter_value)


#formula

trait_table_input <-"../input/GCMP_trait_table_with_abundances_and_adiv_and_metadata_and_growth_data_pcoa.csv"
#Pick the columns to analyze from the various tissue compartments (mucus, skeleton, tissue, all)
#x traits can be: will be the bacteria of interest (copied strait from trait table)
x_trait_column <-"all_D_0__Bacteria___D_1__Proteobacteria___D_2__Gammaproteobacteria___D_3__Oceanospirillales___D_4__Endozoicomonadaceae___D_5__Endozoicomonas"
#only extract the last part of the name
filter_taxa_name <-str_extract(x_trait_column,"[^_]+$")

#y trait disease can be: sum_dis, sum_healthy, sum_bbd, sum_wb, perc_healthy, perc_dis, perc_bbd, perc_wd
y_trait_column<-"growth_rate_mm_per_year"

trait_table <- read.csv(trait_table_input,header=TRUE, row.names=1)

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
#need the new genus level tree!
tree <- "../input/huang_roy_genus_tree.newick"
tree <- read.tree(tree)


#make sure that data and tree use same species names
coral_obj <- name.check(tree,trait_table_small)
coral_obj

#take out tree tips that do not have data
#we have data not in tree for some species. These do not have disease data
tree.prune<-drop.tip(tree, coral_obj$tree_not_data)
tree.prune_obj <-name.check(tree.prune,trait_table_small)
tree.prune_obj

#delete the extra rows from the object
#trait_table.prune<-trait_table_small[-which(rownames(trait_table_small)%in% tree.prune_obj$data_not_tree),]

#commented out parts that are not needed when you don't need to prune the table
#prune_tree_data <- name.check(tree.prune,trait_table.prune)
#prune_tree_data
#plot(tree.prune)

#trait_table.prune$host_genus<-rownames(trait_table.prune)

trait_table_small$host_genus<-rownames(trait_table_small)

#analyze data using the package caper which combines the data and tree into one file

#comp.data <- comparative.data(tree.prune, trait_table.prune, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)
comp.data <- comparative.data(tree.prune, trait_table_small, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)
model.formula <- as.formula(paste0(x_trait_column,'~',y_trait_column))

#run data using the basic model
model.3_formula<-"pgls(x_trait~y_trait, data=comp.data)"
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
write.table(output_row.model.3, file="../output/alpha_diversity_taxa/GCMP_taxa_alpha_diversity_all_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#run data setting lambda to ML
model.3l_formula<-"pgls(x_trait~y_trait, data=comp.data, lambda=ML,kappa=1,delta=1)"
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
write.table(output_row.model.3l, file="../output/alpha_diversity_taxa/GCMP_taxa_alpha_diversity_all_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)

#run data setting delta to ML
model.3d_formula<-"pgls(x_trait~y_trait, data=comp.data, lambda=1,kappa=1,delta=ML)"
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
write.table(output_row.model.3d, file="../output/alpha_diversity_taxa/GCMP_taxa_alpha_diversity_all_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)

#run data setting kappa to ML
model.3k_formula<-"pgls(x_trait~y_trait, data=comp.data, lambda=1,kappa=ML,delta=1)"
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
write.table(output_row.model.3k, file="../output/alpha_diversity_taxa/GCMP_taxa_alpha_diversity_all_pgls.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)

#I don't think these will require these outputs so stop here!
#run phylomorphspace and PICs on 
#create this file only the first time then we append to it
pic.output = data.frame("Disease_Trait", "Diversity_Trait", "Filtered_Column", "Filtered_Value", "N_Microbes", "Compartment", "Model", "pVal", "R_Squared")
write.table(pic.output,file="../output/beta_diversity_taxa/GCMP_taxa_beta_diversity_all_pic.csv",append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#need to redefine the traits for some reason.
#comment out section if don't need to prune the trait table
trait_table.prune<-trait_table_small[-which(rownames(trait_table_small)%in% tree.prune_obj$data_not_tree),]
#trait_table.prune<-trait_table_small
#Pick the columns to analyze
#y trait disease can be: sum_dis, sum_healthy, sum_bbd, sum_wb, perc_healthy, perc_dis, perc_bbd, perc_wd
y_trait_column <-"perc_wd"
#x traits can be: observed_features_compartment, dominance_compartment, gini_index_compartment, simpson_e_compartment, faith_pd_compartment
x_trait_column<-"skeleton_weighted_unifrac_ordination_PC3"
#trait_table <- read.csv(trait_table_input,header=TRUE,row.names=1)
x_trait <-trait_table.prune[[x_trait_column]]
y_trait <-trait_table.prune[[y_trait_column]]

#make plots of the traits and save as a pdf 
phylomorphospace(tree.prune,trait_table.prune[,c(x_trait_column,y_trait_column)],xlab=x_trait_column,ylab=y_trait_column)
pdf(paste("../output/beta_diversity_taxa/", filter_value_last,"_",x_trait_column,"_",y_trait_column,"_beta_diversity_phylomorphospace.pdf",sep=""))
phylomorphospace(tree.prune,trait_table.prune[,c(x_trait_column,y_trait_column)],xlab=x_trait_column,ylab=y_trait_column)
points(x_trait,y_trait,pch=21,bg="grey",cex=1.4)
dev.off()


# Give them names
names(x_trait) <- rownames(trait_table.prune)
names(y_trait) <- rownames(trait_table.prune)

#for contrasts, you should positivize them, since the order doesn't matter. This is NOT taking absolute value.
PositivizeContrasts <- function(x_trait, y_trait) {
  #Cache the sign of x so it doesn't
  #change as we reflect points about x-axis!
  sign_of_x <- sign(x_trait)
  x_trait.positivized <- x_trait * sign_of_x
  y_trait.positivized <- y_trait * sign_of_x
  return(cbind(x_trait.positivized, y_trait.positivized))
}




# Calculate PICs
x_trait_pic <- pic(x_trait, tree.prune)
#print(x_trait_pic)
y_trait_pic <- pic(y_trait, tree.prune)
#print(y_trait_pic)

positivized.results <- PositivizeContrasts(x_trait_pic, y_trait_pic)
x_trait_positive <-positivized.results[,1]
y_trait_positive <-positivized.results[,2]

# Make a pic model
pic_model <- lm(x_trait ~ y_trait - 1)
# plot pic results
plot(x_trait ~ y_trait,xlab=x_trait_column,ylab=y_trait_column,bg='gray',pch=16)
abline(a = 0, b = coef(pic_model))

# Make a pic model
pic_model <- lm(x_trait_positive ~ y_trait_positive - 1)
# plot pic results
plot(x_trait_positive ~ y_trait_positive,xlab=x_trait_column,ylab=y_trait_column,bg='gray',pch=16)
abline(a = 0, b = coef(pic_model))

#Add a 95% CI to the PIC
pic_df <- data.frame(x_trait_positive,y_trait_positive)
ggplot(pic_df, aes(x_trait_positive,y_trait_positive)) +
  geom_smooth(method = "lm", se = TRUE, col = "black", formula = y~x -1)+
  geom_point(size = 3, col = "firebrick")+
  labs(x = paste("Contrast in ", x_trait_column), y = paste("Contrast in ", y_trait_column))+
  theme_classic()

#Save raw PIC contrasts as a pdf
pdf(paste("../output/beta_diversity_taxa/",filter_value_last,"_",x_trait_column,"_",y_trait_column,"_beta_diversity_pic_scatter_YX.pdf",sep=""))
plot(x_trait_positive ~ y_trait_positive,xlab=x_trait_column,ylab=y_trait_column,bg='gray',pch=16)
abline(a = 0, b = coef(pic_model))
dev.off()

#save the PIC contrasts with the 95% confidence interval as a pdf
pdf(paste("../output/beta_diversity_taxa/",filter_value_last,"_",x_trait_column,"_",y_trait_column,"_beta_diversity_pic_scatter_YX_confidence.pdf"))
ggplot(pic_df, aes(x_trait_positive,y_trait_positive)) +
  geom_smooth(method = "lm", se = TRUE, col = "black", formula = y~x -1)+
  geom_point(size = 3, col = "firebrick")+
  labs(x = paste("Contrast in ", x_trait_column), y = paste("Contrast in ", y_trait_column))+
  theme_classic()
dev.off()

#summarize the results
summary(pic_model)
rSquared_PIC <-summary(pic_model)$r.squared
pval_PIC <- anova(pic_model)$'Pr(>F)'[1]
#pval_PIC.rlm <- anova(pic_rlm)$'Pr(>F)'[1]
model_PIC <- "lm(y_trait_positive ~ x_trait_positive -1)"

#save results to a table

output_row.pic <- data.frame(y_trait_column, x_trait_column, filter_column, filter_value, n_microbe, compartment, model_PIC, pval_PIC, rSquared_PIC)
print(output_row.pic)

#Write it to the csv file
write.table(output_row.pic, file="../output/beta_diversity_taxa/GCMP_taxa_beta_diversity_all_pic.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


# plot contMap
x_trait_reconstruction<-contMap(tree.prune,x_trait)
plot(x_trait_reconstruction,direction="rightwards")
#save plot as pdf
pdf(paste("../output/beta_diversity_taxa/",filter_value_last,"_",x_trait_column,"_beta_diversity_contmap.pdf",sep=""))
plot(x_trait_reconstruction,direction="rightwards")
dev.off()

y_trait_reconstruction<-contMap(tree.prune,y_trait)
plot(y_trait_reconstruction,direction="leftwards")
#save plot as pdf
pdf(paste("../output/beta_diversity_taxa/",filter_value_last,"_",y_trait_column,"_beta_diversity_contmap.pdf",sep=""))
plot(y_trait_reconstruction,direction="leftwards")
dev.off()


#plot at likelyhood surface of the lambda parameter
lm.lk<-pgls.profile(model.5, which="lambda")
plot(lm.lk)

#Here we want to run the correction on only the best fit models.
#Is there a way to group by best model?
#correct p-values from the csv file for multiple comparisons
#extract columns we want
pgl_input_table <- "./symbiodinium_pics/PGL_Sym_output_its2&clade.csv"
pgl_table <- read.csv(pgl_input_table,header=TRUE)
pgl_table.df <- data.frame(pgl_table)
p_value_columns <- pgl_table.df[,c("Immune_Trait", "Diversity_Trait","pVal_Model_1", "pVal_Model_2", "pVal_Model_3", "pVal_Model_4")]
#p_value_columns
#Select Immune Trait Values
#Immune traits can be: TIR_total, TIR_total_unique, LRR_total, LRR_total_unique, Lectin_total, Lectin_unique
Immune_trait_variable <- p_value_columns[p_value_columns$Immune_Trait == "Lectin_total",]
#Immune_trait <- p_value_columns$Immune_trait_variable
Immune_trait_variable
#Select by Diversity Trait Value
#Diversity traits can be: obs_asvs, ln_asvs, dominance, gini_index, simpson_e, faith_pd
Diversity_trait_variable <- p_value_columns[p_value_columns$Diversity_Trait =="simpson_e",]
Diversity_trait_variable
#Select the p value column you will be using
#p value columns can be "p_val_Model_1","p_val_Model_2","p_val_Model_3","p_val_model_4"
p_value_select <- Diversity_trait_variable$pVal_Model_1
p_value_select
#run p adjust for bonferroni
adjusted_bonferroni <- p.adjust(p_value_select,method="bonferroni")
adjusted_bonferroni
#Might not need this
#run p adjust for hochberg
adjusted_hochberg <- p.adjust(p_value_select,method="hochberg")
adjusted_hochberg
#run p adjust for BH (fdr)
adjusted_bh <- p.adjust(p_value_select,method="BH")
adjusted_bh

