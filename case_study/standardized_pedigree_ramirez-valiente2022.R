library(tidyr)
library(dplyr)

setwd("~/Documents/GitHub/MetaAnalysis_LocalAdaptation/Ramirez-Valiente_Quercus_2022_SNP_pedigree")

data <- read.csv("Ramirez-Valiente_Quercus_Faginea_GH_2022_Phenotypes.csv", sep = ";", header=T, dec = ",")
#There seems to be a mistake in the naming of the two treatments: there should only be WW and D treatment
#I'll replace the W treatment by WW
data$Treatment[data$Treatment=="W"] <- "WW"
for (i in levels(as.factor(data$Family))) {
  if (length(unique(data$Population[data$Family==i]))!=1){
    print("the family identifiers are not unique for each family")
    print(i)
  }
}
pop <- data$Population
data <- data %>% unite("Family", Population, Family)

p <- data.frame(id = data$ID, pop=pop, family=data$Family)
p <- distinct(p)

create.pedigree <- function(data, family.type=c("half.sib.dam", "half.sib.sire", "full.sib"), same.pop=T) {
  ## data is a data frame with the following columns:
  ## ID, population of origin, family
  
  pedigree <- matrix(nrow = nrow(data), ncol = 5)
  colnames(pedigree) <- c("id", "sire", "dam", "sire_pop", "dam_pop")
  
  pedigree[,1] <- data[,1] #set the ID into the pedigree data frame
  if (family.type=="half.sib.dam" & same.pop) { 
    #if the half sib families are maternal families
    pedigree[,2] <- 1:nrow(pedigree) #sire
    pedigree[,3] <- data[,3] #dam
    pedigree[,4] <- data[,2] #sire_pop
    pedigree[,5] <- data[,2] #dam_pop
  } else if (family.type == "half.sib.sire" & same.pop) {
    #if the half sib families are paternal families
    pedigree[,2] <- data[,3] #sire
    pedigree[,3] <- 1:nrow(pedigree) #dam
    pedigree[,4] <- data[,2] #sire_pop
    pedigree[,5] <- data[,2] #dam_pop
  } else if (family.type=="full.sib" & same.pop) { 
    #if the families are full sib
    pedigree[,2] <- data[,3] #sire
    pedigree[,3] <- data[,3] #dam
    pedigree[,4] <- data[,2] #sire_pop
    pedigree[,5] <- data[,2] #dam_pop
  } 
  return(as.data.frame(pedigree))
}

pedigree <- create.pedigree(p, family.type = "half.sib.dam")

write.csv(pedigree, "ramirez-valiente2022_faginea_GH_pedigree_standardized.csv", row.names = FALSE)

data <- read.csv("Ramirez-Valiente_Quercus_Faginea_OC_2022_Phenotypes.csv", sep = ";", header=T, dec = ",")
for (i in levels(as.factor(data$Family))) {
  if (length(unique(data$Population[data$Family==i]))!=1){
    print("the family identifiers are not unique for each family")
    print(i)
  }
}

p <- data.frame(id = data$ID, pop=data$Population, family=data$Family)
p <- distinct(p)

pedigree <- create.pedigree(p, family.type = "half.sib.dam")

write.csv(pedigree, "ramirez-valiente2022_faginea_OC_pedigree_standardized.csv", row.names = FALSE)

data <- read.csv("Ramirez-Valiente_Quercus_Lusitanica_GH_2022_Phenotypes.csv", sep = ";", header=T, dec = ",")
for (i in levels(as.factor(data$Family))) {
  if (length(unique(data$Population[data$Family==i]))!=1){
    print("the family identifiers are not unique for each family")
    print(i)
  }
}
pop <- data$Population
data <- data %>% unite("Family", Population, Family)

p <- data.frame(id = data$ID, pop=pop, family=data$Family)
p <- distinct(p)

pedigree <- create.pedigree(p, family.type = "half.sib.dam")

write.csv(pedigree, "ramirez-valiente2022_lusitanica_GH_pedigree_standardized.csv", row.names = FALSE)
