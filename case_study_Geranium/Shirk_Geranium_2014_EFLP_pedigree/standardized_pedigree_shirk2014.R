library(tidyr)
library(dplyr)

setwd("MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree/doi_10_5061_dryad_nm41q__v20140605")
#This sudy includes two experiment in a common garden in two separated files

data1 <- read.csv("Gercar qst data.csv", header = T) #experiment 1
IDd <- data.frame(ID=1:length(data1$pop))
data1 <- cbind(IDd, data1)

for (i in levels(as.factor(data1$fam))) {
  if (length(unique(data1$pop[data1$fam==i]))!=1){
    print("the family identifiers are not unique for each family")
    print(i)
  }
}

population <- data1$pop
data1 <- data1 %>% unite("family", pop, fam, sep = "_")

p <- data.frame(id = data1$ID, pop=population, family=data1$family)

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

setwd("/Users/samuelmichel/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree")
write.csv(pedigree, "shirk2014_exp1_pedigree_standardized.csv", row.names = FALSE)

setwd("/Users/samuelmichel/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree/doi_10_5061_dryad_nm41q__v20140605")

data2 <- read.table("Gercar qst SI data.txt", header = T) #experiment 2
IDd <- data.frame(ID=1:length(data2$pop))
data2 <- cbind(IDd, data2)

for (i in levels(as.factor(data2$fam))) {
  if (length(unique(data2$pop[data2$fam==i]))!=1){
    print("the family identifiers are not unique for each family")
    print(i)
  }
}

population <- data2$pop
data2 <- data2 %>% unite("family", pop, fam, sep = "_")

p <- data.frame(id = data2$ID, pop=population, family=data2$family)

pedigree <- create.pedigree(p, family.type = "half.sib.dam")

setwd("/Users/samuelmichel/Documents/GitHub/MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree")
write.csv(pedigree, "shirk2014_exp2_pedigree_standardized.csv", row.names = FALSE)
