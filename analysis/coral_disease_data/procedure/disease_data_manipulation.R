##Data manipulation for complementary disease data##
#Combining FRRP, HICORDIS and Lamb disease data
#Trial with combing FRRP & HICORDIS on just disease vs. not diseased, summarised by genus

#Make sure you're in the working folder "coral_disease_data"
setwd("~/Documents/OSUDocs/Projects/Disease_LHS/GCMP_Global_Disease/analysis/coral_disease_data")
library(plyr)
library(dplyr)

##We start with the FRRP and the HICRODIS since their data format is the same

#First fix up the frrp datasheet and get the host genus & species names
frrp <- read.csv("input/FRRP_DRM/FRRP_DRM_2005_2020_downloaded_10082021.csv", fileEncoding="UTF-8-BOM"
)
frrp_name_map <- read.csv("input/FRRP_DRM/FRRP_DRM_name_mappings.csv", fileEncoding="UTF-8-BOM"
)
frrp_disease_map <- read.csv("input/FRRP_DRM/frrp_disease_name_mappings.csv", fileEncoding="UTF-8-BOM"
)

#make a corresponding column named species_code
#frrp <-rename(frrp, Species = species_code)
frrp <- rename(frrp, species_code = Species)
#frrp_name_map_read <- colnames (frrp_name_map)[1]
#frrp_name_map_read
frrp_name_map <- rename(frrp_name_map, species_code = species_code)
#frrp_name_map <- rename(frrp_name_map, species_code = ï..species_code)

#merge with mapping data to incorporate the host_genus and host_species names 
#according to their species code
frrp <- merge(x = frrp, y = frrp_name_map, by.x = "species_code", by.y= "species_code")
#merge with mapping data to incorporate centralized disease names
frrp$Disease.Conditions <- na_if(frrp$Disease.Conditions, "" ) #changes all empty cells to NAs
frrp <- merge(x=frrp, y=frrp_disease_map, by = "Disease.Conditions")
#frrp <- merge(x=frrp, y=frrp_disease_map, by.x = "Disease.Conditions", by.y = "ï..Disease.Conditions")

##let's change the diseases to a binary - yes or no/diseased vs healthy
frrp$disease_binary <- ifelse(frrp$disease_name != "healthy", "diseased", "healthy")

#let's make binary columns for each disease name to get counts (trial with common disease)
#trial healthy, black band, and white disease
frrp$bbd <- ifelse(frrp$disease_name == "black band disease", 1, 0)
frrp$wd <- ifelse(frrp$disease_category == "white disease", 1, 0)


#Now bring in the HICORDIS data
hicordis <- read.csv("input/HICORDIS/inputs/Caldwell_et_al_HICORDIS_data.csv", fileEncoding="UTF-8-BOM")
hicordis_disease_map <- read.csv("input/HICORDIS/inputs/hicordis_disease_name_mappings.csv", fileEncoding="UTF-8-BOM")

#add centralized disease names 
hicordis <- merge(hicordis, hicordis_disease_map, by = "Disease_Type")

#let's try to get the data in the same format as frrp
hicordis$disease_binary <- ifelse(hicordis$Disease_Type != "No_Disease", "diseased", "healthy")
hicordis <- rename(hicordis, c("host_genus" = "Genus", "host_name" = "Species"))

#There is a wrong genus name in HICORDIS so here we fix it
hicordis$host_genus <- as.factor(hicordis$host_genus)
levels(hicordis$host_genus)[levels(hicordis$host_genus) == "Psammacora"] <- "Psammocora"
levels(hicordis$host_genus)

#let's make binary columns for each disease name to get counts (trial with diseases present across all datasets)
#trial healthy, black band and white diseases
hicordis$bbd <- ifelse(hicordis$disease_name == "black band disease", 1, 0)
hicordis$wd <- ifelse(hicordis$disease_category == "white disease", 1, 0)

#let's extract the disease binary data with the species data for each dataset
frrp_ex <- select(frrp, "host_genus", "host_name", "disease_binary", "bbd", "wd" )
hicordis_ex <- select(hicordis, "host_genus", "host_name", "disease_binary", "bbd", "wd")

#Now that each dataset has exactly the same variables, we can do an rbind
#this merges the data
disease_data <- rbind(frrp_ex, hicordis_ex)
#Since each row is an observation, make a "count" variable so that we can summarise
disease_data$count_total <- 1 #all observations get a 1
disease_data$count_disease <- ifelse(disease_data$disease_binary != "healthy", 1, 0) #only diseased gets a 1
disease_data$count_healthy <- ifelse(disease_data$disease_binary == "healthy", 1, 0)

#Let's look at summarizing based on genus both the N of diseased & healthy
sum_frrp_hicordis_genus <- ddply(disease_data, "host_genus", summarise,
             sum_total = sum(count_total),
             sum_dis = sum(count_disease),
             sum_healthy = sum(count_healthy),
             sum_bbd = sum(bbd),
             sum_wd = sum(wd)) #can add here specific diseases??

View(sum_frrp_hicordis_genus)

##let's bring in Lamb data (non-public dataset)
#This is in a different format than the others, each line is a reef and each species has the
#total number + total diseased number so we can add it in here
lamb <- read.csv("../../../joleahs_private_disease_data/Australia_Disease_Data.csv") #this is hidden & not on github
lamb_ex <- select(lamb, "Genus", "TotalCorals", "TotalDiseased", "TotalHealthy", "TotalBBD", "TotalWS")
str(lamb_ex)
#Two species names are miss spelled or in the case of the Pacific Montastraea the genus name has changed
lamb_ex <- lamb_ex %>%
  mutate(Genus = case_when(Genus == "Montastrea" ~ "Paramontastraea", 
                               Genus == "Tubastrea" ~ "Tubastraea",
                               TRUE~Genus))
#It is not easy to change the names of the diseases for Lamb data bc they are column names
#For reference (or to change in the future):
#TotalBBD = black band disease
#TotalCompromisedORBleached = compromised
#TotalGA = growth anomaly
#TotalWS = white syndrome
#TotalBrB = brown band
#TotalSEB = skeletal eroding band
#TotalAtN = atramentous necrosis

#rename the columns to match with the other two data sets, both with total counts & with individual diseases
lamb_ex <- rename(lamb_ex, c("host_genus" = "Genus", "count_total" = "TotalCorals", "count_disease" = "TotalDiseased",
                              "count_healthy" = "TotalHealthy", "bbd" = "TotalBBD", "wd" = "TotalWS"))

sum_lamb_genus <- ddply(lamb_ex, "host_genus", summarise,
                      sum_total = sum(count_total),
                      sum_dis = sum(count_disease),
                      sum_healthy = sum(count_healthy),
                      sum_bbd = sum(bbd),
                      sum_wd = sum(wd))


#Combine the FRRP/HICORDIS with the Lamb
sum_dis_genus <- rbind(sum_frrp_hicordis_genus, sum_lamb_genus)


#Sum the columns
sum_dis_genus <- ddply(sum_dis_genus, "host_genus", summarise,
                       sum_total = sum(sum_total),
                       sum_dis = sum(sum_dis),
                       sum_healthy = sum(sum_healthy),
                       sum_bbd = sum(sum_bbd),
                       sum_wd = sum(sum_wd))

#Check the percentage of diseased observations
sum_dis_genus$perc_healthy <- (sum_dis_genus$sum_healthy/sum_dis_genus$sum_total)*100
sum_dis_genus$perc_dis <- (sum_dis_genus$sum_dis/sum_dis_genus$sum_total)*100
sum_dis_genus$perc_bbd <- (sum_dis_genus$sum_bbd/sum_dis_genus$sum_total)*100
sum_dis_genus$perc_wd <- (sum_dis_genus$sum_wd/sum_dis_genus$sum_total)*100
View(sum_dis_genus)

#Write your CSV by genus - it will include counts for total healthy, total diseased
#plus total black band and total white disease (includes white plague, white syndrome, white band, white pox)
write.csv(sum_dis_genus, "output/disease_by_genus.csv")

