setwd("MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree/doi_10_5061_dryad_nm41q__v20140605")
library(hierfstat)
library(dplyr)

dat <- read.csv("Gercar allozyme genotypes.csv", header = T, sep=",")
#here I assumed that each loci is encoded by one column
#some loci have very little polymorphism (i included them all but one might need to removed some of them for specific analysis)
dosage <- fstat2dos(dat[,-c(1,2)])
ID <- 1:nrow(dosage)
dosage <- cbind(ID, dat[,1:2], dosage)

setwd("MetaAnalysis_LocalAdaptation/Shirk_Geranium_2014_EFLP_pedigree")
write.csv(dosage, "shirk2014_genotypes_dosage.csv", row.names = F)
