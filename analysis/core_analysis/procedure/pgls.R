#Phylogenetic generalized least squares regression and phylogenetic generalized ANOVA
#load libraries
#library(ape)
#library(nlme)
library(geiger)
library(caper)
#library(MuMIn)
#library(rr2)
#library(phytools)
library(stringr)
library(MASS)
#library(robustbase)
library(ggplot2)
#library(remotes)



print("Usage: Rscript pgls.r <path_to_trait_table> <path_to_tree> <x_trait_column> <y_trait_column> <filter_column> <filter_value> <output_dir>")

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
output_dir_for_all_pgls = "../output/PIC_results"
if (!is.na(user_specified_output_dir) & (user_specified_output_dir != 'None')){
  output_dir_for_all_pgls = paste0(args[7],"/")
}
print("Creating output directory for all PGLS results")
dir.create(output_dir_for_all_pgls, showWarnings = FALSE)


#For this project we co-organize PIC and PGLS results together
#TODO: for other projects you may want to make 'PIC' into 'PGLS' here
output_dir <- paste0(output_dir_for_all_pgls,"PIC_",x_trait,"_vs_",y_trait)

#Add the filter column and value to directory name, if they were provided
if (!is.na(filter_column) & (filter_column != 'None') & (filter_value != 'None')){
  output_dir <- paste0(output_dir,"_",filter_column,"_is_",gsub("\\;","_",filter_value))
}



#Add the trailing slash
output_dir <- paste0(output_dir,"/")
print(paste("Output results to:",output_dir))
dir.create(output_dir,showWarning = FALSE)

#print("Starting log file...")
#Output results to a log file and screen
#NOTE: currently the python parser relies on stdout so logging breaks it
#sink(paste0(output_dir,x_trait,"_vs_",y_trait,"_PGLS_results_log.txt"))

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

#y trait disease can be: sum_dis, sum_healthy, sum_bbd, sum_wb, perc_healthy, perc_dis, perc_bbd, perc_wd
y_trait_column<-y_trait

trait_table <- read.csv(trait_table_fp,header=TRUE, row.names=1,sep="\t")

#filter trait table by the column and values listed (add a ! before trait_table after the comma to invert)
if (!is.na(filter_column) & filter_column != "None" & filter_value != "None"){
  print("Filtering trait table based on user-specified filter_column and filter_value")
  trait_table <- subset(trait_table,trait_table[[filter_column]] == filter_value)
}

#extract the compartment to be used, number of microbes, package

#compartment <-str_extract(x_trait_column,"^[^_]+(?=_)")
compartment <- "-"

n_microbe <-sum(trait_table$sum_total,na.rm=TRUE)

package <- "caper"

#reduce trait table to just columns used in the analysis

print("About to create a reduced trait table")
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

print(paste0("Coral object:",coral_obj))

#take out tree tips that do not have data in the trait table
#(we have trait data for some species that are not in the tree)

#Skip pruning the tree if there are no non-overlapping taxa
if (!is.atomic(coral_obj)){
    pruned_tree<-drop.tip(tree, coral_obj$tree_not_data)
}else{
    pruned_tree <- tree
}

missing_taxa <-name.check(pruned_tree,trait_table_small)

if (!is.atomic(missing_taxa) > 0){
#delete the extra rows from the trait table that are not in the tree
filtered_trait_table<-trait_table_small[-which(rownames(trait_table_small)%in% missing_taxa$data_not_tree),]
}else{
filtered_trait_table <- trait_table_small
}
#commented out parts that are not needed when you don't need to prune the table
pruned_tree_data <- name.check(pruned_tree,filtered_trait_table)

filtered_trait_table$host_genus<-rownames(filtered_trait_table)


#analyze data using the package caper which combines the data and tree into one file
comp.data <- comparative.data(pruned_tree, filtered_trait_table, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)
model.formula.as.text <- paste0(y_trait_column,'~',x_trait_column)
model.formula <- as.formula(model.formula.as.text)

#Define a function to run PGLS and (mostly) organize the results into a dataframe

run_pgls <- function(formula,data,...,model_name="model",lambda=1,kappa=1,delta=1,x_trait="x_trait",y_trait="y_trait",output_dir=NA){
    print(paste0("Running PGLS on model:",model_name))
    bounds = list(lambda = c(1e-6,1.0),kappa=c(1e-6,1),delta=c(1e-6,1))
    model_fit <- pgls(formula,data,lambda=lambda,kappa=kappa,delta=delta,bounds=bounds)
    print("Model fit:")
    print(model_fit)
    
    #Save diagnostic plots
    if (!is.na(output_dir)){
       print(paste0("Outputting PGLS plots to:",output_dir))
       par(mfrow=c(2,2))
       pdf(paste0(output_dir,x_trait,"_vs_",y_trait,"_",model_name,"_pgls_diagnostics.pdf")) 
       plot(model_fit)
       dev.off()
       par(mfrow=c(1,1))
    }
    parameter_description <-paste0("lambda=",lambda," delta=",delta,"kappa=",kappa)
    print(paste0("Summary(",model_name,")"))
    
    print(summary(model_fit))
    aicc <- model_fit$aicc
    print(paste0("AICc:",aicc))
    aic <- model_fit$aic
    print(paste0("AIC:",aic))
    r_squared <- summary(model_fit)$r.squared
    print(paste0("R2:",r_squared))    

    #Extract model coefficients, standard deviations, and 95% CIs
    model_coef <- coef(summary(model_fit))
    #Model coefficients are organized into columns
    #To avoid magic numbers, these are recorded below
    means_col_idx <- 1
    sd_col_idx <-2
    p_value_col_idx <-4

    intercept_row_idx <- 1
    x_trait_row_idx <-2 

        
    intercept <- model_coef[intercept_row_idx,means_col_idx]
    x_trait_slope <- model_coef[x_trait_row_idx,means_col_idx]
    
    intercept_stdev <- model_coef[intercept_row_idx,sd_col_idx]
    x_trait_slope_stdev <- model_coef[x_trait_row_idx,sd_col_idx]

    intercept_lower_CI95 <- intercept - intercept_stdev*1.96
    intercept_upper_CI95 <- intercept + intercept_stdev*1.96

    x_trait_slope_lower_CI95 <- x_trait_slope - x_trait_slope_stdev*1.96
    x_trait_slope_upper_CI95 <- x_trait_slope + x_trait_slope_stdev*1.96
    print(paste0("intercept = ",intercept,"(95% CI ",intercept_lower_CI95,"-",intercept_upper_CI95,")"))
    x_trait_slope_95CI = paste0(x_trait_slope_lower_CI95," - ",x_trait_slope_upper_CI95)

    p <- model_coef[x_trait_row_idx,p_value_col_idx]    
    #extract the branch length parameter
    #TODO: this should really be a loop
    #NOTE: does not accomodate multiple branch length transformations

    estimated_parameter <- "None"
    if (lambda == "ML"){
        estimated_parameter <- "lambda"
    }
    if (delta == "ML"){
        estimated_parameter <- "delta"
    }
    if (kappa == "ML"){
        estimated_parameter <- "kappa"
    }
  
    if (estimated_parameter == "None"){
        parameter_value = 1
        parameter_CI95 = "Fixed"
        parameter_estimate_text = "All parameters fixed" 
    }
    if (estimated_parameter != "None"){
        parameter_value = model_fit$param[estimated_parameter]
        parameter_CI95 = pgls.confint(model_fit,estimated_parameter)$ci.val
        parameter_estimate_text = paste(estimated_parameter,":",parameter_value,"(95% CI ",parameter_CI95[1]," - ",parameter_CI95[2],")")
    }


   
    
    results_dataframe <- data.frame(model_name = model_name,compartment=compartment,branch_length_transformation = parameter_description,AIC = aic, AICc = aicc,estimated_parameter = estimated_parameter,parameter_value = parameter_estimate_text,R2 = r_squared,slope=x_trait_slope,intercept=intercept,x_trait_slope_95CI=x_trait_slope_95CI,x_trait_slope_stdev = x_trait_slope_stdev,p = p)
    results_dataframe$x_trait <- x_trait
    results_dataframe$y_trait <- y_trait
    return(results_dataframe)  
}

model_label = "BM"
bm_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda=1,kappa=1,delta=1,x_trait=x_trait_column,y_trait=y_trait_column)
print("BM results:")
print(bm_results)

model_label = "BM_Lambda"
bm_lambda_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda="ML",kappa=1,delta=1,x_trait=x_trait_column,y_trait=y_trait_column,output_dir=output_dir)
print("BM + Lambda Results")
print(bm_lambda_results)

model_label = "BM_Kappa"
bm_kappa_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda=1,kappa="ML",delta=1,x_trait=x_trait_column,y_trait=y_trait_column,output_dir=output_dir)
print("BM + Kappa results")
print(bm_kappa_results)

model_label = "BM_Delta"
bm_delta_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda=1,kappa=1,delta="ML",x_trait=x_trait_column,y_trait=y_trait_column,output_dir = output_dir)
print("BM Delta Results")
print(bm_delta_results)

table_template = data.frame(x_trait=c(),y_trait = c(),compartment=c(),model_name = c(), branch_length_transformation = c(), AIC = c(), AICc = c(),p= c(),R2=c(),estimated_parameter = c(),parameter_value = c(),R2 = c())

combined_results_df = do.call(rbind, list(bm_results, bm_lambda_results, bm_kappa_results,bm_delta_results ))

min_aicc_index <- which.min(combined_results_df$AICc)
print(paste("Index of lowest AICc:",min_aicc_index))

best_model <- combined_results_df$model_name[min_aicc_index]
print(paste("Min AICc Model:",best_model))

is_best_model <- combined_results_df$model_name == best_model

combined_results_df$best_model <- is_best_model


pgls_results_file = paste0(output_dir,"PGLS_results.tsv")

print(paste0("Writing PGLS results file:",pgls_results_file))
write.table(combined_results_df,file=pgls_results_file,append=FALSE, sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)


quit()

