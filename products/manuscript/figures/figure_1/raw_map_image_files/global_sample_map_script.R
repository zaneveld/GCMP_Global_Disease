##Script for making a global map with sampling locations##
##Set working directory
#setwd("GCMP_Global_Disease/products/manuscript/figures/")

#Load libraries
library(tidyverse)
library(rvest)
library(magrittr)
library(ggmap)
library(stringr)

##For final mapping code go to line 161##

#Bring in world map data
world <- map_data("world")
world

#Plot the map to make sure it works 
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region)
  ) 

#bring in GCMP data in sample summary

gcmp <- read.csv("../supplementary_materials/gcmp_sample_summary.csv", header = TRUE)

#Check to see if it is easy to map the gcmp data:
ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id = region),
    color = "white", fill = "gray", size = 0.1) +
  geom_point(data = gcmp, aes(Longitude, Latitude, color = Coral.Genus,size=N),
    alpha = 0.5, position="dodge") +
  #labs(x = NULL, y = NULL, color = NULL)+
  theme_void() +
  theme(legend.position = "none")

#Position dodge doesn't seem to work on map - might be another way to offset the overlapping dots
#but not sure right now

ggsave("/figure_1/worldmap_withlocation.pdf", plot = last_plot())


##Can we get the number of genera and the sample size for each dot?
sum.gcmp <- ddply(gcmp, c("Reef.Name", "Longitude", "Latitude"), summarise,
                  NumGen = length(Coral.Genus),
                  sum.tot = sum(N)
)

#then plot dots sized by number of genera and/or sample size:
ggplot() + geom_map(data = world, map = world, aes(long, lat, map_id = region),
                   color = "white", fill = "gray", size = 0.1) +
  geom_point(data = sum.gcmp, aes(Longitude, Latitude, color = Reef.Name, size=NumGen),
    alpha = 0.5, show.legend = FALSE) +
  geom_text() +
  labs(x = NULL, y = NULL, color = NULL)+
  theme_void() #+
  #theme(legend.position = "none")


#Let's attempt making pie charts to show genera by sample size:
#First we need to make our data in long format
long.gcmp <- gcmp %>% spread(Coral.Genus, N) #make a column for each coral genus
long.gcmp[is.na(long.gcmp)] <- 0 #set NAs to 0
long.gcmp$Sample.Size <- rowSums(long.gcmp[,7:76]) #Add back a sample size column since we split this up by genus

View(long.gcmp) #check it all worked

#Trial setting some colours - we need 69! 
mycolors <- c('#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032',
              '#C2B280',  '#222222','#848482',  '#E68FAC', '#0067A5', '#BF40BF',
              '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#E37383',
              '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26', 'F89880',
              '#FFDEAD', '#FFBF00', '#A42A04', '#722F37', '#5D3FD3', '#301934',
              '#EC5800', '#DAA520', '#8A9A5B', '#DFFF00', '#097969', '#818589',
              '#8A9A5B', '#6082B6', '#483C32', '#A0522D', '#CC7722', '#40B5AD', 
              '#93C572', '#B4C424', '#FFFAA0', '#ECFFDC', '#FBCEB1', '#FF2400',
              '#E0115F', '#FF3131', '#702963', '#C3B1E1', '#A95C68', '#FF00FF', 
              '#C9A9A6', '#DC143C', '#811331', '#FFAC1C', '#00FF7F', '#808000',
              '#00FFFF', '#C0C0C0', "#D2B48C", '#8B4513', '#4A0404', '#D27D2D',
              '#2AAA8A', '#0BDA51', '#FF69B4', '#F89880', '#E30B5C')

ggplot() + geom_map(data = world, map = world, aes(long, lat, map_id = region),
                    color = "white", fill = "gray", size = 0.1) +
  geom_scatterpie(data = long.gcmp, aes(x=Longitude, y = Latitude, group=Reef.Name, r = (Sample.Size/10), color = mycolors),
                  cols=colnames(long.gcmp[,7:76]), color= "Black", alpha=0.6, show.legend = FALSE) +
  labs(x = NULL, y = NULL, color = NULL)+
  theme_void()

#the colors didn't work - it just stays at default
#perhaps there is another way to set colors for the pies?

##What about combining by political area to make the map more readable

#First average the lat and long by political area (is this okay to do??)
polarea.gcmp <- ddply(gcmp, c("Ocean", "Political.Area"), summarise,
                      Pol.lat = mean(Latitude), 
                      Pol.long = mean(Longitude))
#Get sample sizes
genus.gcmp <- ddply(gcmp, c("Coral.Genus", "Political.Area"), summarise,
                    Sample.Size = sum(N))
#Combine the two
polarea.gcmp <- left_join(polarea.gcmp, genus.gcmp, by = "Political.Area")

#Make into long form as above
polarea.gcmp.long <- polarea.gcmp %>% spread(Coral.Genus, Sample.Size)
polarea.gcmp.long[is.na(polarea.gcmp.long)] <- 0 #set NAs to 0
polarea.gcmp.long$Sample.Size <- rowSums(polarea.gcmp.long[,5:74]) #Add back in a sample size since we split it up
View(polarea.gcmp.long)


ggplot() + geom_map(data = world, map = world, aes(long, lat, map_id = region),
                    color = "white", fill = "gray", size = 0.1) +
  geom_scatterpie(data = polarea.gcmp.long, aes(x=Pol.long, y = Pol.lat, group=Political.Area, color = mycolors, r = 10),
                  cols=colnames(polarea.gcmp.long[,5:74]), color= NA, alpha=0.8, show.legend = FALSE) +
  labs(x = NULL, y = NULL, color = NULL)+
  theme_void()

#Still didn't do the colour thing, so set to default. 
ggsave("worldmap_with_polarea_pies.pdf", plot = last_plot())
# run again with "show.legend = TRUE" to get the legend. I don't usually export
# with the legend and the graph as the legend is too big, it squishes up the map
ggsave("worldmap_polarea_pies_legend.pdf", plot = last_plot())


#The two locations in the Caribbean/Central American Pacific are almost on top of each other
#according to lat and long (and plotted on a small map)
#Let's see if we can shuffle these a little, and make the pies bigger so it's more visible.

#Change the lat long of Panama to be a little more into the Pacific
#And the lat long of Colombia to the middle of the Caribbean

polarea.gcmp.long.newcoord <- polarea.gcmp.long
polarea.gcmp.long.newcoord$Pol.lat[polarea.gcmp.long.newcoord$Political.Area == "Panama"] <- 4.6
polarea.gcmp.long.newcoord$Pol.long[polarea.gcmp.long.newcoord$Political.Area == "Panama"] <- -89.0

polarea.gcmp.long.newcoord$Pol.lat[polarea.gcmp.long.newcoord$Political.Area == "Colombia"] <- 17.1
polarea.gcmp.long.newcoord$Pol.long[polarea.gcmp.long.newcoord$Political.Area == "Colombia"] <- -78.0

#Plot
ggplot() + geom_map(data = world, map = world, aes(long, lat, map_id = region),
                    color = "white", fill = "gray", size = 0.1) +
  geom_scatterpie(data = polarea.gcmp.long.newcoord, aes(x=Pol.long, y = Pol.lat, group=Political.Area, color = mycolors, r = 10),
                  cols=colnames(polarea.gcmp.long.newcoord[,5:74]), color= NA, alpha=0.8, show.legend = FALSE) +
  labs(x = NULL, y = NULL, color = NULL)+
  theme_void()

ggsave("worldmap_with_polarea_pies_newcoord.pdf", plot = last_plot())


##Check the stats
#Sample size for each location + # of genera to put on the map (in illustrator)

gcmp.stats <- ddply(gcmp, c("Ocean", "Political.Area"), summarise,
                    Sample.size = sum(N), 
                    NGen = length(Coral.Genus))
View(gcmp.stats)


######################
######Final code!#####
######################
##What about making the pies families and ordered by the tree (how do we do this?)
#Can we change the colours to make warm colours robust and cool colours complex?

#First bring in the full dataset

all.gcmp <- read.csv("~/Desktop/Desktop Misc/GCMP_decontaminated_1000_sample-metadata.csv", header = TRUE)
View(all.gcmp)

all.gcmp$reef_name <- as.factor(all.gcmp$reef_name)
levels(all.gcmp$reef_name)

all.gcmp <- all.gcmp %>% filter(host_genus_id != "Not applicable")
all.gcmp <- all.gcmp %>% filter(host_genus_id != "Missing: Not collected")

all.gcmp.group <- ddply(all.gcmp, c("family", "host_genus_id", "complex_robust", "ocean", "political_area", "site_name", "reef_name", "latitude", "longitude"), summarise,
                      N = length(sample_name_mf)
)
View(all.gcmp.group)
#family
#complex_robust

#remove outgroups
#remove Missing data
all.gcmp.group <- all.gcmp.group %>% filter(complex_robust != "outgroup")
all.gcmp.group <- all.gcmp.group %>% filter(family != "Missing: Not collected")


all.gcmp.group$family <- factor(all.gcmp.group$family, levels = c("Acroporidae", "Agariciidae", "Astrocoeniidae",
                                                                      "Dendrophylliidae", "Discosomatidae", "Euphylliidae", "Poritidae",
                                                                      "Siderastreidae",
                                                                      "Coscinaraeidae", "Diploastreidae", "Fungiidae", "incertae sedis", "Lobophylliidae", 
                                                                      "Merulinidae", "Montastreidae", "Mussidae", "Pocilloporidae", "Psammocoridae"))

#First average the lat and long by political area (is this okay to do??)
all.gcmp.group$latitude <- as.numeric(all.gcmp.group$latitude)
all.gcmp.group$longitude <- as.numeric(all.gcmp.group$longitude)
polarea.gcmp.all <- ddply(all.gcmp.group, c("ocean", "political_area"), summarise,
                      Pol.lat = mean(latitude), 
                      Pol.long = mean(longitude))
#Get sample sizes
family.gcmp.all <- ddply(all.gcmp.group, c("family", "political_area"), summarise,
                    Sample.Size = sum(N))
#Combine the two
polarea.gcmp.all <- left_join(polarea.gcmp.all, family.gcmp.all, by = "political_area")

#Make into long form as above
polarea.gcmp.all.long <- polarea.gcmp.all %>% spread(family, Sample.Size)
polarea.gcmp.all.long[is.na(polarea.gcmp.all.long)] <- 0 #set NAs to 0
polarea.gcmp.all.long$Sample.Size <- rowSums(polarea.gcmp.all.long[,5:21]) #Add back in a sample size since we split it up
View(polarea.gcmp.all.long)

#Make a list of colours, using cool colours for complex and warm colours for robust
#Since the families should already be in order of complex-robust, we should be able to just order the 16
#different colours to match
#8 cool colours
#10 warm colours

#clade.colors <- c("#7209B7", "#560BAD", "#480CA8", "#3A0CA3","#3F37C9", "#4361EE", "#4895EF", "#4CC9F0",
#                  "#FFB950", "#FFAD33", "#FF931F", "#FF7E33", "#FA5E1F", "#EC3F13", "#B81702", "#A50104", "#8E0103", "#7A0103")

##Okay need to manually apply the colours to the families using the scall_fill_manual 
#Need to make warm colours opposite!

#The two locations in the Caribbean/Central American Pacific are almost on top of each other
#according to lat and long (and plotted on a small map)
#Let's see if we can shuffle these a little, and make the pies bigger so it's more visible.

#Change the lat long of Panama to be a little more into the Pacific
#And the lat long of Colombia to the middle of the Caribbean

polarea.gcmp.all.long.newcoord <- polarea.gcmp.all.long
polarea.gcmp.all.long.newcoord$Pol.lat[polarea.gcmp.all.long.newcoord$political_area == "Panama"] <- 4.6
polarea.gcmp.all.long.newcoord$Pol.long[polarea.gcmp.all.long.newcoord$political_area == "Panama"] <- -89.0

polarea.gcmp.all.long.newcoord$Pol.lat[polarea.gcmp.all.long.newcoord$political_area == "Colombia"] <- 17.1
polarea.gcmp.all.long.newcoord$Pol.long[polarea.gcmp.all.long.newcoord$political_area == "Colombia"] <- -78.0


ggplot() + geom_map(data = world, map = world, aes(long, lat, map_id = region),
                    color = "white", fill = "gray", size = 0.1) +
  geom_scatterpie(data = polarea.gcmp.all.long.newcoord, aes(x=Pol.long, y = Pol.lat, group=political_area, color = mycolors, r = 8),
                  cols=colnames(polarea.gcmp.all.long.newcoord[,5:21]), color= NA, alpha=0.8, show.legend = TRUE) +
  #scale_fill_manual(values=sample(clade.colors)) +
  scale_fill_manual(values = c("Acroporidae" = "#7209B7", "Agariciidae" = "#560BAD", "Astrocoeniidae" = "#480CA8", "Dendrophylliidae" = "#3A0CA3",
                               "Discosomatidae" = "#3F37C9", "Euphylliidae" = "#4361EE", "Poritidae" = "#4895EF", "Siderastreidae" = "#4CC9F0",
                               "Coscinaraeidae" = "#FFB950", "Diploastreidae" = "#FFAD33", "Fungiidae" = "#FF931F", "incertae sedis" = "#FF7E33", 
                               "Lobophylliidae" = "#FA5E1F", "Merulinidae" = "#EC3F13", "Montastreidae" = "#B81702", "Mussidae" = "#A50104",
                               "Pocilloporidae" = "#8E0103", "Psammocoridae" = "#7A0103")) +
  labs(x = NULL, y = NULL, color = NULL)+
  theme_void()


ggsave("figure_1/worldmap_with_polarea_pies_newcoord_cladecolor.pdf", plot = last_plot())







