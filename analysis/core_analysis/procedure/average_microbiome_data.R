
##Task is to create a new data table for all numeric microbial data averaged by genus

library(dplyr)

##Make sure you are in the "core_analysis/procedure" folder, then read in the mapping file
map <- read_tsv("../input/GCMP_EMP_map_r28_no_empty_samples.txt")

#The column classes in some instances are incorrect - true numeric values are classified as character 
#due to the presence of string characters in numeric columns (e.g., "Missing:Not collected" and "Unknown")
#So we need to change these columns to numeric

#Let's find the column numbers of all the columns so that we can pass it to our function
col.numbers <- colnames(map)
View(col.numbers) #Unfortunately the columns need to be chosen manually and passed to an object "i"

i <- c(8,18,22,23,24,25,26,27,32,48,49,50,54,57,59,62,65,68,
       71,74,86,87,110,112,114,116,118,120,127,
       128,130,132,134,137,140,146,148,149,150,151,152,153)

#The following will re-classify the columns above into numeric values
#FYI: it will coerce some of the string data (e.g., "Missing:Not collected) into NAs and produce warnings
#This is okay for what we are doing downstream
map[ , i] <- apply(map[ , i], 2,        
                    function(x) as.numeric(as.character(x)))
#This will coerce some of the string data into NAs and produce warnings - this is okay

#Check that the code worked
sapply(map, class)

#We can use summarise_if to get the mean of the numeric columns only
#WARNING: This will subset to just the numeric columns.
map.mean <- map %>%
  group_by(map$host_genus) %>%
  summarise_if(is.numeric, mean)

write_tsv(map.mean, "../output/GCMP_EMP_map_r28_no_empty_samples_numeric_only_averaged_by_genus.txt")
