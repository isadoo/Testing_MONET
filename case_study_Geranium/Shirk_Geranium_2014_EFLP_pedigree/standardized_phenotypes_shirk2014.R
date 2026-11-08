library(tidyr)
library(dplyr)

setwd("MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree/doi_10_5061_dryad_nm41q__v20140605")
#This sudy includes two experiment in a common garden in two separated files

data1 <- read.csv("Gercar qst data.csv", header = T) #experiment 1
IDd <- data.frame(ID=1:length(data1$pop))
data1 <- cbind(IDd, data1)
column_names <- colnames(data1)[6:length(data1)] 

reshaped_data <- pivot_longer( 
  data1,
  cols = all_of(column_names), 
  names_to = "trait_id", 
  values_to = "trait"
)

final_data <- reshaped_data %>%
  select(individual = ID, block = blk, trait, trait_id) %>%
  arrange(trait_id, individual)


head(final_data)
unique(final_data$trait_id)

setwd("/Users/samuelmichel/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree")
write.csv(final_data, "shirk2014_exp1_phenotypes_standardized.csv", row.names = FALSE)

setwd("/Users/samuelmichel/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree/doi_10_5061_dryad_nm41q__v20140605")

data2 <- read.table("Gercar qst SI data.txt", header = T) #experiment 2
IDd <- data.frame(ID=1:length(data2$pop))
data2 <- cbind(IDd, data2)
column_names <- colnames(data2)[9:length(data2)] 

reshaped_data <- pivot_longer( 
  data2,
  cols = all_of(column_names), 
  names_to = "trait_id", 
  values_to = "trait"
)

final_data <- reshaped_data %>%
  select(individual = ID, block = blk, trait, trait_id) %>%
  arrange(trait_id, individual)


head(final_data)
unique(final_data$trait_id)

setwd("/Users/samuelmichel/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree")
write.csv(final_data, "shirk2014_exp2_phenotypes_standardized.csv", row.names = FALSE)
