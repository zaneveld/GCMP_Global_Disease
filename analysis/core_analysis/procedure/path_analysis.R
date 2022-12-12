##R Code for phylogenetic causality analysis##

library(ggplot2)
library(phylopath)
library(ape)

#Make sure you are in the analysis/core_analysis/procedure directory

#import data and tree
coral_disease_data <- read.table("../output/GCMP_trait_table_with_abundances_and_adiv_and_metadata_and_growth_data.tsv", header=TRUE, sep="\t")
coral_tree <- read.tree('../output/huang_roy_genus_tree.newick')
rownames(coral_disease_data) <- coral_disease_data$X
filtered_table <- subset(coral_disease_data, select = c("growth_rate_mm_per_year","tissue_D_0__Bacteria___D_1__Proteobacteria___D_2__Gammaproteobacteria___D_3__Oceanospirillales___D_4__Endozoicomonadaceae___D_5__Endozoicomonas","perc_dis"))

# new column names:
# GR = growth_rate_mm_per_year
# EN = tissue...Endozoicomonas
# DS = perc_dis

colnames(filtered_table) <- c("GR","EN","DS")
filtered_table

#Define the possible (biologically logical) models:
# __ notates causation ( B <- A, aka B is caused by A )
models <- define_model_set(
  EN__DS_plus_GR    = c(EN ~ DS + GR),
  EN__GR__DS  = c(EN ~ GR, GR ~ DS),
  GR__DS__EN   = c(GR ~ DS, DS ~ EN),
  GR__DS_plus_EN  = c(GR ~ DS + EN),
  DS__EN__GR  = c(DS ~ EN, EN ~ GR),
  DS__GR__EN  = c(DS ~ GR, GR ~ EN),
  EN__GR    = c(EN ~ GR),
  GR__EN = c(GR ~ EN),
  GR__DS = c(GR ~ DS),
  EN__DS  = c(EN ~ DS),
  DS__EN = c(DS ~ EN),
  DS__GR = c(DS ~ GR),
  DS_and_EN__GR = c(EN ~ GR, DS ~ GR),
  DS_and_GR__EN = c(DS ~ EN, GR ~ EN)
)

#Plot models & save
plot_model_set(models)
ggsave("../output/phylopath/model_set.pdf", plot = last_plot())

#Evaluation of hypotheses using both Brownian (BM) and Pagel's Lamba (PL) evolutionary models
p_BM <- phylo_path(models, data = filtered_table, tree = coral_tree, model = "BM")
p_PL <- phylo_path(models, data = filtered_table, tree = coral_tree, model = "lambda")

#Summary results & graph of support for each model 
#Brownian
s_BM <- summary(p_BM)
s_BM
write.csv(s_BM, "../output/phylopath/phylopath_summary_BM.csv")
plot(s_BM)
ggsave("../output/phylopath/phylopath_model_summary_BM.pdf", plot = last_plot())

#Pagel's Lambda
s_PL <- summary(p_PL)
s_PL
write.csv(s_PL, "../output/phylopath/phylopath_summary_PL.csv")
plot(s_PL)
ggsave("../output/phylopath/phylopath_model_summary_PL.pdf", plot = last_plot())


#Selecting and fitting a final model
#First look at the "best" model fit (although be aware of CICc criteria)
#Brownian Motion
best_model_BM <- best(p_BM)
best_model_BM
plot(best_model_BM)
ggsave("../output/phylopath/phylopath_best_model_BM.pdf", plot = last_plot())

#Pagel's Lambda
best_model_PL <- best(p_PL)
best_model_PL
plot(best_model_PL)
ggsave("../output/phylopath/phylopath_best_model_PL.pdf", plot = last_plot())

#Next check the "average" models, which likely better represent the causality
#due to the multiple model fits as determined by CICc

#Average models:
#Brownian
average_model_BM <- average(p_BM, avg_method = "full")
plot(average_model_BM, algorithm = "mds", curvature = 0.1)
ggsave("../output/phylopath/phylopath_average_model_BM.pdf", plot = last_plot())

#Pagel's Lambda
average_model_PL <- average(p_PL, avg_method = "full")
plot(average_model_PL, algorithm = "mds", curvature = 0.1)
ggsave("../output/phylopath/phylopath_average_model_PL.pdf", plot = last_plot())











