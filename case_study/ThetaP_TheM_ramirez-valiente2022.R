library(tidyr)
library(dplyr)

options(bitmapType = "cairo")
devtools::load_all("/users/ijeronim/mywork/LAVA")

setwd("/users/ijeronim/mywork/MetaAnalysis_LocalAdaptation/Ramirez-Valiente_Quercus_2022_SNP_pedigree")

### Quercus faginea Greenhouse experiment ---------------------------------------------------
F1_pedigree <- read.csv("ramirez-valiente2022_faginea_GH_pedigree_standardized.csv")
parent_genotypes <-read.csv("ramirez-valiente2022_faginea_GH_genotype_dosage.csv")

## F1
if (identical(F1_pedigree$sire_pop, F1_pedigree$dam_pop)) {
  
  selected_pop <- intersect(unique(parent_genotypes$pop), unique(F1_pedigree$dam_pop))
  F1_pedigree <- F1_pedigree %>% 
    filter(dam_pop %in% selected_pop) %>%
    mutate(dam_pop=as.character(dam_pop), sire_pop=as.character(sire_pop)) %>%
    arrange(dam_pop)
  F1_pedigree <- F1_pedigree %>% 
    group_by(dam_pop) %>% 
    arrange(dam, .by_group = TRUE) %>% 
    ungroup()
  F1_id_pop <- F1_pedigree %>% select(id, population_id = dam_pop) 
  
}
#Number of pop
n_pop_F1 <- length(unique(F1_id_pop$population_id))

## Parent
parent_genotypes <- parent_genotypes %>% 
  filter(pop %in% selected_pop) %>%
  mutate(pop=as.character(pop)) %>%
  arrange(pop)
parent_id_pop <- parent_genotypes %>% select(id=id, population_id=pop)
#Number of pop
n_pop_parent <- length(unique(parent_id_pop$population_id))
#parent dosage
parent_dosage <- parent_genotypes[,-c(1:2)]

## population_individual_id
if (n_pop_parent == n_pop_F1 & n_pop_parent==length(unique(c(parent_id_pop$population_id, F1_id_pop$population_id)))) {
  print(paste("There is", n_pop_F1, "populations"))
  id_pop_parent_F1 <- rbind(parent_id_pop, F1_id_pop)
} else {
  print("Number of population is not equivalent between dataframes!")
}

coancestries <- calculate_coancestries(genetic_data_parents=parent_dosage,
                                       genotyped_parent_populations=parent_id_pop$population_id, 
                                       genetic_data_F1 = NA, 
                                       population_individual_id = id_pop_parent_F1, 
                                       column_individual = "id", 
                                       column_population = "population_id", 
                                       pedigree = F1_pedigree,
                                       all_parents_genotyped = FALSE)

##Visual inspection
png("ThetaP_ramirez-valiente2022_faginea_GH.png", width = 1000, height = 1000, type = "cairo")
image(coancestries$Theta.P, axes=F)
dev.off()

png("TheM_ramirez-valiente2022_faginea_GH.png", width = nrow(coancestries$M), height = nrow(coancestries$M), type = "cairo")
image(coancestries$M, axes=F)
dev.off()

image(coancestries$M[1:200,1:200], axes=F)

##Saving data
saveRDS(coancestries$Theta.P, file = "ramirez-valiente2022_faginea_GH_ThetaP.rds")
saveRDS(coancestries$M, file = "ramirez-valiente2022_faginea_GH_TheM.rds")

### Quercus faginea Outdoor common garden experiment -----------------------------------------
# Individuals come from the same wild populations as the GH experiment, so the parental
# coancestry matrix (Theta.P) can be reused. We only need to compute the kinship matrix
# (TheM) from the OC pedigree.
F1_pedigree_OC <- read.csv("ramirez-valiente2022_faginea_OC_pedigree_standardized.csv")

## Order pedigree by population (and by dam within population), matching the GH workflow
if (identical(F1_pedigree_OC$sire_pop, F1_pedigree_OC$dam_pop)) {
  F1_pedigree_OC <- F1_pedigree_OC %>%
    mutate(dam_pop = as.character(dam_pop), sire_pop = as.character(sire_pop)) %>%
    arrange(dam_pop) %>%
    group_by(dam_pop) %>%
    arrange(dam, .by_group = TRUE) %>%
    ungroup()
}

## Compute kinship matrix from pedigree
TheM_OC <- kinship_from_pedigree(F1_pedigree_OC[, c("id", "sire", "dam")])

## Visual inspection
png("TheM_ramirez-valiente2022_faginea_OC.png", width = nrow(TheM_OC), height = nrow(TheM_OC), type = "cairo")
image(TheM_OC, axes = F)
dev.off()

## Saving data
saveRDS(TheM_OC, file = "ramirez-valiente2022_faginea_OC_TheM.rds")

### Quercus lusitanica Greenhouse experiment ---------------------------------------------------
F1_pedigree <- read.csv("ramirez-valiente2022_lusitanica_GH_pedigree_standardized.csv")
parent_genotypes <-read.csv("ramirez-valiente2022_lusitanica_GH_genotype_dosage.csv")

## F1
if (identical(F1_pedigree$sire_pop, F1_pedigree$dam_pop)) {
  
  selected_pop <- intersect(unique(parent_genotypes$pop), unique(F1_pedigree$dam_pop))
  F1_pedigree <- F1_pedigree %>% 
    filter(dam_pop %in% selected_pop) %>%
    mutate(dam_pop=as.character(dam_pop), sire_pop=as.character(sire_pop)) %>%
    arrange(dam_pop)
  F1_pedigree <- F1_pedigree %>% 
    group_by(dam_pop) %>% 
    arrange(dam, .by_group = TRUE) %>% 
    ungroup()
  F1_id_pop <- F1_pedigree %>% select(id, population_id = dam_pop) 
  
}
#Number of pop
n_pop_F1 <- length(unique(F1_id_pop$population_id))

## Parent
parent_genotypes <- parent_genotypes %>% 
  filter(pop %in% selected_pop) %>%
  mutate(pop=as.character(pop)) %>%
  arrange(pop)
parent_id_pop <- parent_genotypes %>% select(id=id, population_id=pop)
#Number of pop
n_pop_parent <- length(unique(parent_id_pop$population_id))
#parent dosage
parent_dosage <- parent_genotypes[,-c(1:2)]

## population_individual_id
if (n_pop_parent == n_pop_F1 & n_pop_parent==length(unique(c(parent_id_pop$population_id, F1_id_pop$population_id)))) {
  print(paste("There is", n_pop_F1, "populations"))
  id_pop_parent_F1 <- rbind(parent_id_pop, F1_id_pop)
} else {
  print("Number of population is not equivalent between dataframes!")
}

coancestries <- calculate_coancestries(genetic_data_parents=parent_dosage,
                                       genotyped_parent_populations=parent_id_pop$population_id, 
                                       genetic_data_F1 = NA, 
                                       population_individual_id = id_pop_parent_F1, 
                                       column_individual = "id", 
                                       column_population = "population_id", 
                                       pedigree = F1_pedigree,
                                       all_parents_genotyped = FALSE)

##Visual inspection
png("ThetaP_ramirez-valiente2022_lusitanica_GH.png", width = 1000, height = 1000, type = "cairo")
image(coancestries$Theta.P, axes=F)
dev.off()

png("TheM_ramirez-valiente2022_lusitanica_GH.png", width = nrow(coancestries$M), height = nrow(coancestries$M), type = "cairo")
image(coancestries$M, axes=F)
dev.off()

image(coancestries$M[1:200,1:200], axes=F)

##Saving data
saveRDS(coancestries$Theta.P, file = "ramirez-valiente2022_lusitanica_GH_ThetaP.rds")
saveRDS(coancestries$M, file = "ramirez-valiente2022_lusitanica_GH_TheM.rds")
