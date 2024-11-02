#Plots for growth rate by genus
library(ggplot2)
library(tidyverse)
#Bring in the combined trait table with growth rate

gcmp.all <- read.csv("../input/GCMP_trait_table_with_abundances_and_adiv_and_metadata_and_growth_data_pcoa.csv")

#Identify the indices for the two columns we need. 
colnames(gcmp.all)
#host_genus_id [393]
#growth_rate_mm_per_year [404]
growth <- gcmp.all[,c(393,404)]


growth.sort <- arrange(growth, desc(growth_rate_mm_per_year))
factor_names <- as.list(growth.sort$host_genus_id)

bar_plot <-
  growth.sort %>%
  ggplot(aes(y = growth_rate_mm_per_year, x = host_genus_id)) +
  geom_bar(stat =  "identity", fill = "#A54657") + 
  #scale_fill_manual(values = c("Acropora" = '#F3C300', "Dipsastraea" = '#008856', "Gardineroseris" = '#875692', 
  #                 "Goniastrea"= '#F38400', "Hydnophora" = '#A1CAF1', "Isopora" = '#BE0032', 
  #                  "Merulina" = '#C2B280',  "Montastraea" = '#222222',"Montipora" = '#848482',  "Orbicella" = '#E68FAC', '#0067A5', 
  #                  "Pavona" = '#F99379', "Platygyra" = '#604E97', "Pocillopora" = '#F6A600', "Porites" = '#B3446C', '#DCD300', 
  #                  "Psammacora" = '#882D17', "Steriatopora" = '#8DB600', "siderastrea" = '#654522', 
  #                  "Stylophora" = '#E25822', "Turbinaria" = '#2B3D26')+
  xlab("Coral Genus") +
  ylab("Growth Rate (mm/year)")+
  labs(fill = '') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(hjust = 0, size = 6))

bar_plot$data$host_genus_id <- factor(x = bar_plot$data$host_genus_id, levels = factor_names) #reorders the plot data by names extracted above
bar_plot

ggsave(bar_plot, file = "../output/growth_by_genus.pdf", width = 7, height = 5) #needs to be long to see the names okay

#perc_dis [9]

#Can we do a combined scatter plot?

#we have 19 taxa with growth rate 

kelly_colors = c('#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032', 
                 '#C2B280',  '#222222','#848482',  '#E68FAC', '#0067A5', 
                 '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', 
                 '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')

scatter_plot <- ggplot(gcmp.all, aes(y = perc_dis, x = growth_rate_mm_per_year, group = host_genus_id)) + 
  geom_point(aes(color = host_genus_id), size = 4) +
  scale_color_manual(values = kelly_colors) +
  xlab("Disease Prevalence (%)") +
  ylab("Growth Rate (mm/year)")+
  labs(color = 'Coral Genus') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(hjust = 0, size = 12))

scatter_plot

library(cowplot)

growth_plots<-cowplot::plot_grid(bar_plot, scatter_plot, labels = "auto", label_fontface = "plain")
growth_plots
ggsave(growth_plots, file = "../output/growth_disease_plots/growth_by_disease_plots.pdf", width = 12, height = 5) #needs to be long to see the names okay
ggsave(growth_plots,file = "../output/growth_disease_plots/growth_by_disease_plots.jpg", width = 12, height = 5 )



#Disease bars
gcmp.dis <- read.csv("../input/GCMP_trait_table_with_abundances_and_adiv_and_metadata_zeros.csv")
dis <- gcmp.dis[,c(9, 393)] 
dis.sort <- arrange(dis, desc(perc_dis))
dis.sort <- na.omit(dis.sort)
dis_factor_names <- as.list(dis.sort$host_genus_id)
dis.sort$host_genus_id <- as.factor(dis.sort$host_genus_id)


dis.bar <-
  ggplot(dis.sort, aes(y = perc_dis, x = host_genus_id))+
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

dis.bar$data$host_genus_id <- factor(x = dis.bar$data$host_genus_id, levels = dis_factor_names) #reorders the plot data by names extracted above
dis.bar

print(dis_factor_names)

ggsave(bar, file = "../output/disease_susc_by_genus.jpg", width = 12, height = 5) #needs to be long to see the names okay
ggsave(bar, file = "../output/disease_susc_by_genus.pdf", width = 12, height = 5) #needs to be long to see the names okay



