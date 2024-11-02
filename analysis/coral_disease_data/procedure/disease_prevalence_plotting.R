##Plots for disease by genus - GCMP ###
library(ggplot2)
library(tidyverse)

#Bring in disease data by genus & cleanup
dis.all <- read.csv("../output/disease_by_genus.csv", header = T)
dis.all <- dis.all[,-1] #remove the first indexing column (X)
dis.all.nona <- na.omit(dis.all) #remove NAs (there are 2)

#What does it look like:
ggplot(dis.all.nona, aes(reorder(host_genus, percent_diseased), percent_diseased)) +
  geom_col(fill = "gray", width = 0.7)+
  coord_flip() +
  labs(x = "Coral Genus", y = "Disease Percentage (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, vjust = 1, hjust = 1), 
        panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))



#What about a stacked column plot where half is diseased and half is healthy?

#extract a sorted order of host genus names by disease prevalence to use to re-order bars at the end
dis.all.nona.sort <- arrange(dis.all.nona, desc(percent_diseased))
factor_names <- as.list(dis.all.nona.sort$host_genus)

#remove everything except the %healthy and %diseased:
dis.health.percent <- dis.all.nona.sort[,c(1,8,9)]
#Melt to a long-form
dis.health.percent.long <- melt(dis.health.percent, id = c("host_genus"))
#plot
stack_bar <-
  dis.health.percent.long%>%
  ggplot(aes(fill = variable, y = value, x = host_genus))+
  geom_bar(position = "fill", stat =  "identity")+ #, colour = "white"
  #scale_fill_brewer(palette = "Spectral")+
  scale_fill_manual(values = c("#A54657", "#F1A66A"), labels = c("Healthy", "Diseased")) +
  xlab("Coral Genus") +
  ylab("Disease Prevalence (%)")+
  labs(fill = '') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(hjust = 0, size = 6))

stack_bar$data$host_genus <- factor(x = stack_bar$data$host_genus, levels = factor_names) #reorders the plot data by names extracted above
stack_bar

ggsave(stack_bar, file = "../output/disease_prevalence_by_genus.jpg", width = 12, height = 5) #needs to be long to see the names okay
ggsave(stack_bar, file = "../output/disease_prevalence_by_genus.pdf", width = 12, height = 5) #needs to be long to see the names okay


#How about just disease
#Remove healthy observations
dis.percent <- dis.all.nona.sort[,c(1,9)]


bar <-
  dis.percent%>%
  ggplot(aes(y = percent_diseased, x = host_genus))+
  geom_col(position = "stack")+ #, colour = "white"
  #scale_fill_brewer(palette = "Spectral")+
 # scale_fill_manual(values = host_genus) +
  xlab("Coral Genus") +
  ylab("Disease Susceptibility (%)")+
  labs(fill = '') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(hjust = 0, size = 8))

bar$data$host_genus <- factor(x = bar$data$host_genus, levels = factor_names) #reorders the plot data by names extracted above
bar

ggsave(bar, file = "../output/disease_susc_by_genus.jpg", width = 12, height = 5) #needs to be long to see the names okay
ggsave(bar, file = "../output/disease_susc_by_genus.pdf", width = 12, height = 5) #needs to be long to see the names okay




#What about one with a tree??

library(phytools)

ctree <- read.tree('../../core_analysis/output/huang_roy_genus_tree.newick')
disease_table <- read.csv("../../core_analysis/input/GCMP_trait_table_with_abundances_and_adiv_and_metadata_zeros.csv", header = T, row.names = 1)
disease <-setNames(disease_table$perc_dis,rownames(disease_table))
disease_color <- c("lightcoral")

pdf(file="../output/total_disease_barplot_with_phylogram.pdf")
plotTree.barplot(ctree, disease, xlab = "Total Disease")
dev.off()


phytools::plotTree.wBars(ctree, disease, type = "phylogram", fsize = 0.2, xlab = "Total Disease")
#, col=disease_color,xlab="Total Disease")
#dev.off()

