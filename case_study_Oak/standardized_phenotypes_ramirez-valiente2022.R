library(tidyr)
library(dplyr)

setwd("~/Documents/GitHub/MetaAnalysis_LocalAdaptation/Ramirez-Valiente_Quercus_2022_SNP_pedigree")

## FAGINEA GH
data <- read.csv("Ramirez-Valiente_Quercus_Faginea_GH_2022_Phenotypes.csv", sep = ";", header=T, dec = ",")
#There seems to be a mistake in the naming of the two treatments: there should only be WW and D treatment
#I'll replace the W treatment by WW
data$Treatment[data$Treatment=="W"] <- "WW"
column_names <- colnames(data)[6:ncol(data)]

reshaped_data <- pivot_longer( 
  data,
  cols = all_of(column_names), 
  names_to = "trait_id", 
  values_to = "trait"
)

reshaped_data <- reshaped_data %>%
  unite("trait_id", trait_id, Treatment, sep=".Treat")

final_data <- reshaped_data %>%
  select(individual = ID, Block, trait, trait_id) %>%
  arrange(trait_id, individual)


head(final_data)
unique(final_data$trait_id)

write.csv(final_data, "ramirez-valiente2022_faginea_GH_phenotypes_standardized.csv", row.names = FALSE)

## FAGINEA OC
data <- read.csv("Ramirez-Valiente_Quercus_Faginea_OC_2022_Phenotypes.csv", sep = ";", header=T, dec = ",")

column_names <- colnames(data)[6:ncol(data)]

reshaped_data <- pivot_longer( 
  data,
  cols = all_of(column_names), 
  names_to = "trait_id", 
  values_to = "trait"
)

final_data <- reshaped_data %>%
  select(individual = ID, Row, Column, trait, trait_id) %>%
  arrange(trait_id, individual)


head(final_data)
unique(final_data$trait_id)

write.csv(final_data, "ramirez-valiente2022_faginea_OC_phenotypes_standardized.csv", row.names = FALSE)

##LUSITANICA GH
data <- read.csv("Ramirez-Valiente_Quercus_Lusitanica_GH_2022_Phenotypes.csv", sep = ";", header=T, dec = ",")

column_names <- colnames(data)[6:ncol(data)]

reshaped_data <- pivot_longer( 
  data,
  cols = all_of(column_names), 
  names_to = "trait_id", 
  values_to = "trait"
)

reshaped_data <- reshaped_data %>%
  unite("trait_id", trait_id, Treatment, sep=".Treat")

final_data <- reshaped_data %>%
  select(individual = ID, Block, trait, trait_id) %>%
  arrange(trait_id, individual)


head(final_data)
unique(final_data$trait_id)

write.csv(final_data, "ramirez-valiente2022_lusitanica_GH_phenotypes_standardized.csv", row.names = FALSE)
