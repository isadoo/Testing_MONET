library(tidyr)
library(dplyr)

source("~/Documents/GitHub/MetaAnalysis_LocalAdaptation/calculate_coancestries.r")
source("~/Documents/GitHub/MetaAnalysis_LocalAdaptation/kinship_from_pedigree.r")

setwd("~/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree")

### EXP 1
F1_pedigree <- read.csv("shirk2014_exp1_pedigree_standardized.csv")
parent_genotypes <- read.csv("shirk2014_genotypes_dosage.csv")

## Parent
parent_genotypes <- parent_genotypes %>% arrange(pop)
parent_id_pop <- parent_genotypes %>% select(id=ID, population_id=pop)
#Number of pop
n_pop_parent <- length(unique(parent_id_pop$population_id))
#parent dosage
parent_dosage <- parent_genotypes[,-c(1:3)]

## F1
if (identical(F1_pedigree$sire_pop, F1_pedigree$dam_pop)) {
  F1_pedigree <- F1_pedigree %>% arrange(dam_pop)
  F1_pedigree <- F1_pedigree %>% 
    group_by(dam_pop) %>% 
    arrange(dam, .by_group = TRUE) %>% 
    ungroup()
  F1_id_pop <- F1_pedigree %>% select(id, population_id = dam_pop) 
}
#Number of pop
n_pop_F1 <- length(unique(F1_id_pop$population_id))

## population_individual_id
if (n_pop_parent == n_pop_F1 & n_pop_parent==length(unique(c(parent_id_pop$population_id, F1_id_pop$population_id)))) {
  print(paste("There is", n_pop_F1, "populations"))
  id_pop_parent_F1 <- rbind(parent_id_pop, F1_id_pop)
} else {
  print("Number of population is not equivalent between dataframes!")
}


coancestries <- calculate_coancestries(genetic_data_parents = parent_dosage,
                                       data_type_parents = "dosage",
                                       genotyped_parent_populations = parent_id_pop$population_id,
                                       genetic_data_F1 = NA, 
                                       data_type_F1 = NA,
                                       number_of_populations = NA, 
                                       population_individual_id = id_pop_parent_F1, 
                                       column_individual = "id", 
                                       column_population = "population_id", 
                                       pedigree = F1_pedigree, 
                                       BiAllelic = FALSE,
                                       founders_phenotyped = FALSE,
                                       all_parents_genotyped = FALSE)

##Visual inspection
png("ThetaP_exp1_shirk2014.png", width = 1000, height = 1000, type = "cairo")
image(coancestries$Theta.P, axes=F)
dev.off()

png("TheM_exp1_shirk2014.png", width = 1000, height = 1000, type = "cairo")
image(coancestries$M, axes=F)
dev.off()

image(coancestries$M[1:200,1:200], axes=F)

##Saving data
saveRDS(coancestries$Theta.P, file = "shirk2014_exp1_ThetaP.rds")
saveRDS(coancestries$M, file = "shirk2014_exp1_TheM.rds")

### EXP 2
F1_pedigree <- read.csv("shirk2014_exp2_pedigree_standardized.csv")

## F1
if (identical(F1_pedigree$sire_pop, F1_pedigree$dam_pop)) {
  F1_pedigree <- F1_pedigree %>% arrange(dam_pop)
  F1_pedigree <- F1_pedigree %>% 
    group_by(dam_pop) %>% 
    arrange(dam, .by_group = TRUE) %>% 
    ungroup()
  F1_id_pop <- F1_pedigree %>% select(id, population_id = dam_pop) 
}
#Number of pop
n_pop_F1 <- length(unique(F1_id_pop$population_id))

## population_individual_id
if (n_pop_parent == n_pop_F1 & n_pop_parent==length(unique(c(parent_id_pop$population_id, F1_id_pop$population_id)))) {
  print(paste("There is", n_pop_F1, "populations"))
  id_pop_parent_F1 <- rbind(parent_id_pop, F1_id_pop)
} else {
  print("Number of population is not equivalent between dataframes!")
}


coancestries <- calculate_coancestries(genetic_data_parents = parent_dosage,
                                       data_type_parents = "dosage",
                                       genotyped_parent_populations = parent_id_pop$population_id,
                                       genetic_data_F1 = NA, 
                                       data_type_F1 = NA,
                                       number_of_populations = NA, 
                                       population_individual_id = id_pop_parent_F1, 
                                       column_individual = "id", 
                                       column_population = "population_id", 
                                       pedigree = F1_pedigree, 
                                       BiAllelic = FALSE,
                                       founders_phenotyped = FALSE,
                                       all_parents_genotyped = FALSE)

##Visual inspection
png("ThetaP_exp2_shirk2014.png", width = 1000, height = 1000, type = "cairo")
image(coancestries$Theta.P, axes=F)
dev.off()

png("TheM_exp2_shirk2014.png", width = 1000, height = 1000, type = "cairo")
image(coancestries$M, axes=F)
dev.off()

image(coancestries$M[1:200,1:200], axes=F)

##Saving data
saveRDS(coancestries$Theta.P, file = "shirk2014_exp2_ThetaP.rds")
saveRDS(coancestries$M, file = "shirk2014_exp2_TheM.rds")
