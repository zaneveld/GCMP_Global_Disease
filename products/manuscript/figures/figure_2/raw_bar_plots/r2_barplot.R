### Code for making an R2 bar plot for publication ###

#R2 of disease vs. variable on y axis and alpha/beta diversity metrics on x axis, grouped by compartment type

#load libraries
library(dplyr)
library(ggplot2)

#read in table (data extracted from PGLS_gcmp_output, R2 from best fitting models for alpha div)
r2.table <- read.csv("../input/gcmp_r2_table.csv", header = TRUE)
#Make sure the order is right
r2.table$diversity_metric <- factor(x = r2.table$diversity_metric, levels = c("richness", "evenness", "dominance"))

ggplot(r2.table, aes(fill=compartment, y=R2, x=diversity_metric)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("#098F94","#e98A5C", "#8D3F95")) +
  theme_classic()

ggsave("output/r2_barplot_alphadiv.pdf")
