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

#Define the models:
# __ notates causation ( B <- A, aka B is caused by A )
models <- define_model_set(
  EN__DS_plus_GR    = c(EN ~ DS + GR),
  #two    = c(EN ~ DS + GR, GR ~ DS),
  EN__GR__DS  = c(EN ~ GR, GR ~ DS),
  GR__DS__EN   = c(GR ~ DS, DS ~ EN),
  GR__DS_plus_EN  = c(GR ~ DS + EN),
  #six   = c(GR ~ DS + EN, EN ~ DS),
  DS__EN__GR  = c(DS ~ EN, EN ~ GR),
  DS__GR__EN  = c(DS ~ GR, GR ~ EN),
  #nine   = c(DS ~ GR + EN, EN ~ GR)
  EN__GR    = c(EN ~ GR),
  GR__EN = c(GR ~ EN),
  GR__DS = c(GR ~ DS),
  EN__DS  = c(EN ~ DS),
  DS__EN = c(DS ~ EN),
  DS__GR = c(DS ~ GR),
  DS_and_EN__GR = c(EN ~ GR, DS ~ GR),
  DS_and_GR__EN = c(DS ~ EN, GR ~ EN)
)

plot_model_set(models)
ggsave("../output/phylopath/model_set.pdf", plot = last_plot())

#Evaluation of hypotheses
p_BM <- phylo_path(models, data = filtered_table, tree = coral_tree, model = "BM")
p_PL <- phylo_path(models, data = filtered_table, tree = coral_tree, model = "lambda")

s_BM <- summary(p_BM)
s_BM
s_PL <- summary(p_PL)
s_PL

#Summary of support for each model

plot(s_BM)
ggsave("../output/phylopath/phylopath_model_summary_BM.pdf", plot = last_plot())
plot(s_PL)
ggsave("../output/phylopath/phylopath_model_summary_PL.pdf", plot = last_plot())

#Selecting and fitting a final model
#Brownian Motion
best_model_BM <- best(p_BM)
best_model_BM
plot(best_model_BM)
ggsave("../output/phylopath/phylopath_best_model_BM.pdf", plot = last_plot())

#Pagal's lambda
best_model_PL <- best(p_PL)
best_model_PL
plot(best_model_PL)
ggsave("../output/phylopath/phylopath_best_model_PL.pdf", plot = last_plot())


#Average models:
average_model_BM <- average(p_BM, avg_method = "full")
plot(average_model_BM, algorithm = "mds", curvature = 0.1)
ggsave("../output/phylopath/phylopath_average_model_BM.pdf", plot = last_plot())
average_model_PL <- average(p_PL, avg_method = "full")
plot(average_model_PL, algorithm = "mds", curvature = 0.1)
ggsave("../output/phylopath/phylopath_average_model_PL.pdf", plot = last_plot())











