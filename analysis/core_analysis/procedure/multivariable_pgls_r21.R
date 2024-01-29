#Phylogenetic generalized least squares regression and phylogenetic generalized ANOVA
#load libraries

library(geiger)
library(caper)

library(stringr)
library(MASS)
library(ggplot2)




print("Usage: Rscript pgls.r <path_to_trait_table> <path_to_tree> <x_trait_columns> <y_trait_column> <filter_column> <filter_value> <output_dir>")
print("Note: x_trait_columns should be comma delimited, like trait1,trait2,trait3, etc")

#Get the user input and assign variables
args <- commandArgs(trailingOnly = TRUE)
trait_table_fp <- args[1]
tree_fp <- args[2]
x_traits <- args[3]
y_trait <- args[4]
filter_column <- args[5]
filter_value <- args[6]
user_specified_output_dir <- args[7]



#Extract a vector of trait column names from
#comma separated list
x_trait_columns <- unlist(strsplit(x_traits, ","))
print(paste0("Found x trait columns:",x_trait_columns))

#Redundantly define a Y trait column 
#for consistency with the X trait variable
y_trait_column <- y_trait

output_dir_for_all_pgls = "./PGLS_results/"

#Create a function to make the overall directory holding all PIC/PGLS results
create_overall_output_dir <- function(output_dir_for_all_pgls,user_specified_output_dir){
	#Set output directory for all the PGLS models
	#results of this analysis will be created in a subdirectory

	if (!is.na(user_specified_output_dir) & (user_specified_output_dir != 'None')){
  		output_dir_for_all_pgls = paste0(user_specified_output_dir,"/")
	}
	print("Creating output directory for all PGLS results")
	dir.create(output_dir_for_all_pgls, showWarnings = FALSE)		
}

create_overall_output_dir(output_dir_for_all_pgls,user_specified_output_dir)



#create function to automatically name the output directory given 
#command line arguments
set_output_dir_name <- function(output_dir_for_all_pics,x_traits,y_trait,filter_column,
  filter_value){
    output_dir <- paste0(output_dir_for_all_pics,"PIC_",x_traits,"_vs_",y_trait)

    #Add filter column and value to directory name, if they were provided
    if (!is.na(filter_column) & (filter_column != 'None') & (filter_value != 'None')){
        output_dir <- paste0(output_dir,"_",filter_column,"_is_",gsub("\\;","_",filter_value))
    }
    
    output_dir <- paste0(output_dir,"/")

    print(paste("Outputting results to:",output_dir))   
    return(output_dir)
}

#Create output subdirectory for this analysis specifically
output_dir <- set_output_dir_name(output_dir_for_all_pgls,x_traits,y_trait,filter_column,
  filter_value)
dir.create(output_dir,showWarning = FALSE)
print(paste0("Created output dir for this analysis:",output_dir))

#Set formula for the regression analysis
formula_str <- paste(y_trait, "~",  paste(x_trait_columns, collapse = " * "))
print(paste0("Setting formula string to:",formula_str))


#print("Starting log file...")
#Output results to a log file and screen
#NOTE: currently the python parser relies on stdout so logging breaks it
#sink(paste0(output_dir,x_trait,"_vs_",y_trait,"_PGLS_results_log.txt"))

print("PGLS Analysis Report")
print("--------------------------------------------------")
print (paste("Analyzing", x_traits, "vs.", y_trait))
print(paste("Trait table filepath:",trait_table_fp))
print(paste("Tree filepath:",tree_fp))
print(paste("Filtering data based on column:",filter_column))
print(paste("Including data only if filter column value is:", filter_value))
print(paste("Output dir:",output_dir))
print(paste("Regression formula:",formula_str))
print("--------------------------------------------------")

#Load trait table
trait_table <- read.table(trait_table_fp,comment.char="",header=T,row.names=1,as.is=T,sep="\t")

#Load tree
tree <- read.tree(tree_fp)

#Generate a vector with all relevant column names
trait_columns = c(y_trait,x_trait_columns)

# Filter the trait table to include only the specified columns
filtered_trait_table <- trait_table[, trait_columns]

#Optionally filter the trait table based on a specific column value
#prior to harmonizing the tree and trait table

#Remove trailing \; characters in filter value
#filter_value <-sub('.*\\;', '', filter_value)

if (!is.na(filter_column) & filter_column != 'None' & filter_value != 'None'){
    trait_table <- subset(trait_table,trait_table[[filter_column]] == filter_value)
}

# Define a function for removing rows with missing data
remove_taxa_with_missing_data <- function(column, tree, trait_table) {
  cat("Analyzing Trait Column:", column, "\n")

  #Coerce trait column to numeric, possibly inducing NAs
  trait_table[,column]<-as.numeric(trait_table[,column])
  #Drop rows with an NA for that trait
  trait_table<- trait_table[!is.na(trait_table[,column]),]
  return(trait_table)
}

#Filter specified trait columns for missing data
for (column in trait_columns) {
  filtered_trait_table <- remove_taxa_with_missing_data(column, tree, filtered_trait_table)
  print(filtered_trait_table)
}

#Define a function to drop tree tips or trait table rows that don't match
harmonize_tree_and_trait_table <-function(trait_table,tree){
	#Filter table to match tree tips
	filtered_trait_table <- trait_table[rownames(trait_table) %in% tree$tip.label,]
	#Filter tree to match table
	filtered_tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% rownames(trait_table)])
	#Dichotomize tree
	tree <- multi2di(filtered_tree)
	return(list(tree = filtered_tree, trait_table = filtered_trait_table))
}



#Drop mismatched data or tips between trait table and tree
harmonized_data <- harmonize_tree_and_trait_table(filtered_trait_table,tree)
trait_table <- harmonized_data$trait_table
tree <- harmonized_data$tree
print(paste0("tree class:",class(tree)))

# Add a host_name column based on rownames
host_name <- rownames(trait_table)
trait_table <- cbind(host_name = host_name, trait_table)

#Record the filtered tree and trait table (as .tsv)
write.tree(tree,paste0(output_dir,x_traits,"_vs_",y_trait,"_filtered_tree.newick"))
write.table(trait_table,paste0(output_dir,x_traits,"_vs_",y_trait,"_filtered_table.tsv"),sep="\t", row.names=FALSE)


package <- "caper"

#reduce trait table to just columns used in the analysis
#n_unique_sample <- length(unique(sample_column))
n_unique_sample <- nrow(trait_table)

#make sure that data and tree use same species names
coral_obj <- name.check(tree,trait_table)

print(paste0("Coral object:",coral_obj))

#Take out tree tips that do not have data in the trait table
#(we have trait data for some species that are not in the tree)
#Skip pruning the tree if there are no non-overlapping taxa
if (!is.atomic(coral_obj)){
    pruned_tree<-drop.tip(tree, coral_obj$tree_not_data)
}else{
    pruned_tree <- tree
}

coral_obj <-name.check(pruned_tree,trait_table)
print(paste0("Missing taxa?",coral_obj))

#name.check returns the text "OK" if all taxa match
#or a list of taxa if there are mismatches.
#So we have to do a little work to figure out which has happened

if (!is.atomic(coral_obj) > 0){
	#delete the extra rows from the trait table that are not in the tree
	filtered_trait_table<-trait_table[-which(rownames(trait_table)%in% coral_obj$data_not_tree),]
}else{
	filtered_trait_table <- trait_table
}

#commented out parts that are not needed when you don't need to prune the table
coral_obj <- name.check(pruned_tree,filtered_trait_table)

filtered_trait_table$host_genus<-rownames(filtered_trait_table)


#analyze data using the package caper which combines the data and tree into one file
comp.data <- comparative.data(pruned_tree, filtered_trait_table, names.col="host_genus", vcv.dim=2, warn.dropped = TRUE)

#Write a formula for the x traits
#model.formula.as.text <- paste(y_trait_column, "~",  paste(x_trait_columns, collapse = " * "))
model.formula.as.text <- paste(y_trait_column, "~",  paste(x_trait_columns, collapse = " + "))
print(paste0("Using formula:",model.formula.as.text))

model.formula <- as.formula(model.formula.as.text)

#Define a function to run PGLS and (mostly) organize the results into a dataframe
print("Defining a PGLS function")

#Define a function to parse PGLS model coefficients (automatically called by run_pgls)
parse_model_coefficients <-function(model_fit){
	
	# Create a list to store coefficients
	coefficients_list <- list()
	
	
	#Model coefficients are organized into columns
    #To avoid magic numbers, these are recorded below
    model_summary <- summary(model_fit)
    model_summary_rownames <- rownames(model_summary)
    
    

	#Model fit is NOT a normal model summary,
	#but instead a special pgls object
	#So we can't easily convert to standard
	#DataFrame
	
	#Get TOTAL number of parameters fit, 
	#including interaction effects and the intercept    
    n_slopes <- model_fit$k
	
	#Extract model coefficients, standard deviations, and 95% CIs
    model_coef <- coef(summary(model_fit))
	model_coef_table <- as.data.frame(model_coef)
	
	statistic_names <- colnames(model_coef_table)
	
	model_coef_table["CI upper"] <- model_coef_table["Estimate"]+model_coef_table["Std. Error"] * 1.96
	model_coef_table["CI lower"] <- model_coef_table["Estimate"]-model_coef_table["Std. Error"] * 1.96
	
	# Loop over coefficients and store them in the list
	
	fit_parameter_names <- rownames(model_coef_table)
	fit_parameter_names <- fit_parameter_names[1:length(fit_parameter_names)]
	model_coef_table$x_trait <-fit_parameter_names
	
	#Add data on R2, AIC, AICc
    aicc <- model_fit$aicc
    aic <- model_fit$aic
    r_squared <- summary(model_fit)$r.squared  
    
    model_coef_table$R2<-r_squared
    model_coef_table$AIC<-as.numeric(aic)
    model_coef_table$AICc<-as.numeric(aicc)
    model_coef_table$n_parameters <- n_slopes
	#Return the model coefficient table as a DataFrame
	return(model_coef_table)
	
	
}



#Define a function to label compartment 
get_compartment <-function(x_trait_columns){
	
	print(x_trait_columns)
	x_traits<-(paste(x_trait_columns,collapse=","))
	print(paste0("x traits:",x_traits))
	compartment <- ""
	for (curr_compartment in c("mucus","tissue","skeleton"))
	{
		if (grepl(curr_compartment,x_traits))
			{
				if (nchar(compartment) >= 1){
					compartment <- paste(compartment,curr_compartment,sep=",")
				}
				
				if (nchar(compartment) < 1){
					compartment <- curr_compartment
				}
				
			}
	}
	return(compartment)
}

#Define a function to run PGLS (and parse model coefficients into a DataFrame)
run_pgls <- function(formula,data,model_name="model",lambda=1,kappa=1,delta=1,x_traits="x_traits",y_trait="y_trait",output_dir=NA){
    
    print(x_traits)
    compartment <- get_compartment(x_traits)
    print(paste0("Compartment (detected from x_traits):",compartment))
    
    print(paste0("Running PGLS on model:",model_name))
    bounds = list(lambda = c(1e-6,1.0),kappa=c(1e-6,1),delta=c(1e-6,1))
    model_fit <- pgls(formula,data,lambda=lambda,kappa=kappa,delta=delta,bounds=bounds)
    
    x_traits_label<-(paste(x_traits,collapse=","))
    #Save diagnostic plots
    if (!is.na(output_dir)){
       print(paste0("Outputting PGLS plots to:",output_dir))
       par(mfrow=c(2,2))
       pdf(paste0(output_dir,x_traits_label,"_vs_",y_trait,"_",model_name,"_pgls_diagnostics.pdf")) 
       plot(model_fit)
       dev.off()
       #Reset subpanel arrangement for next plot
       par(mfrow=c(1,1))
    }
    parameter_description <-paste0("lambda=",lambda," delta=",delta,"kappa=",kappa)
    print(paste0("Summary(",model_name,")"))
    	
	coefficients_df = parse_model_coefficients(model_fit)
	
    #Summarize evolutionary model for results dataframe
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
    coefficients_df$compartment <- compartment
    coefficients_df$x_traits <-x_traits_label
    coefficients_df$model_name <- model_name
    coefficients_df$parameter_value <- parameter_value
    coefficients_df$estimated_parameter <- parameter_estimate_text
    return(coefficients_df)
}
    


print("About to run PGLS")
model_label = "BM"

bm_results <-run_pgls(model.formula,comp.data,model_name=model_label,lambda=1,
  kappa=1,delta=1,x_traits=x_trait_columns,y_trait=y_trait_column,output_dir=output_dir)


model_label = "BM_Lambda"
bm_lambda_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda="ML",kappa=1,delta=1,x_trait=x_trait_columns,y_trait=y_trait_column,output_dir=output_dir)


model_label = "BM_Kappa"
bm_kappa_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda=1,kappa="ML",delta=1,x_trait=x_trait_columns,y_trait=y_trait_column,output_dir=output_dir)

model_label = "BM_Delta"
bm_delta_results <- run_pgls(model.formula,comp.data,model_name=model_label,lambda=1,kappa=1,delta="ML",x_trait=x_trait_columns,y_trait=y_trait_column,output_dir = output_dir)


table_template = data.frame(x_traits=c(),y_trait = c(),compartment=c(),model_name = c(), branch_length_transformation = c(), AIC = c(), AICc = c(),p= c(),R2=c(),estimated_parameter = c(),parameter_value = c(),R2 = c())

combined_results_df = do.call(rbind, list(bm_results, bm_lambda_results, bm_kappa_results,bm_delta_results ))
combined_results_df$formula <- model.formula.as.text
min_aicc_index <- which.min(combined_results_df$AICc)
best_model <- combined_results_df$model_name[min_aicc_index]
is_best_model <- combined_results_df$model_name == best_model

combined_results_df$n <- n_unique_sample
combined_results_df$y_trait <- y_trait
combined_results_df$delta_AICc <- combined_results_df$AICc - min(combined_results_df$AICc)
combined_results_df$best_model <- is_best_model
combined_results_df$filter_column <- filter_column
combined_results_df$filter_value <- filter_value
combined_results_df$trait_table <- trait_table_fp
combined_results_df$tree <- tree_fp

rownames(combined_results_df) <- NULL

#Sort table into specified order
order <- c("model_name","compartment","formula","n_parameters","x_traits","y_trait", "x_trait", "R2", "Pr(>|t|)","Estimate", "Std. Error","t value", "CI lower","CI upper","estimated_parameter","parameter_value","n","AIC","AICc","delta_AICc","best_model","filter_column","filter_value","trait_table","tree")
combined_results_df <- combined_results_df[, order]


print("****************FINAL TABLE**********************")
print(combined_results_df)
print("*************************************************")

pgls_results_file = paste0(output_dir,"PGLS_results.tsv")

print(paste0("Writing PGLS results file:",pgls_results_file))
write.table(combined_results_df,file=pgls_results_file,append=FALSE, sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)


quit()

