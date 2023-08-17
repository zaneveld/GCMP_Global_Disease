### Code for making an R2 bar plot for publication ###

#R2 of disease vs. variable on y axis and alpha/beta diversity metrics on x axis, grouped by compartment type

#load libraries
library(dplyr)
library(ggplot2)

#read in table (data extracted from PGLS_gcmp_output, R2 from best fitting models for alpha div)
r2.table <- read.csv("input/gcmp_r2_table_revised.csv", header = TRUE)
#Make sure the order is right

r2.table$diversity_metric <- factor(x = r2.table$diversity_metric, levels = c("richness", "evenness", "dominance", "WPC1", "WPC2", "WPC3"))

#Group by compartment & fill by diversity metric

ggplot(r2.table, aes(fill=factor(diversity_metric), y=R2, x=diversity_category)) + 
  geom_bar(position="dodge", stat="identity") +
  #scale_fill_manual(values = "#098F94") +
  guides(fill=guide_legend(title="Compartment")) +
  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) +
  theme_classic() +
  facet_wrap(~compartment)

ggsave("output/r2_barplot_color_by_compartment_new.pdf")

#Group by compartment and fill by metric
ggplot(r2.table, aes(fill=factor(compartment, levels = c("mucus", "tissue", "skeleton")), y=R2, x=diversity_metric)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("#098F94", "#e98A5C", "#8D3F95")) +
  guides(fill=guide_legend(title="Compartment")) +
  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) +
  theme_classic() +
  facet_wrap(~compartment)

ggsave("output/r2_barplot_color_by_compartment_grouped_by_metric_new.pdf")

