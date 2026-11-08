# This script is part of the breeding tests with simulated data we will run
# The current script is the one related to the breeding design and number of families.
# We will do 3 dams and 3 sires per population. (2 offspring per mating; total 18 offspring per population)
# For this breeding design we need to sample 6 individuals per population (now we have 10) and then say half is sires and half dams.
# We run for all population structures
# Used in the paper for wdiff 4.6, wvar 10 (both correlations)

#setwd("Chapter2")
source("../../Package/counting_blocks_matrix.r")
source("../../Package/kinship_from_pedigree.r")
source("../../Package/create_F1.r")
source("../../Package/calculate_coancestries.r")

library(hierfstat)
library(gaston)
library(JGTeach)

args <- commandArgs(trailingOnly = TRUE)

# Set parameters
replicate_number <- as.integer(args[1])
Population_structure <- args[2] # "IM"
number_of_pop <- as.integer(args[3]) # 18
generations <- as.integer(args[4]) # 1000 for IM
correlation <- ifelse(args[5] == "NA", "", paste0("_", args[5]))
wdiff <- args[6]
wvar <- args[7]
experiment_name <- args[8]

# Build file paths
folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
Simulation_directory <- paste0("../quantinemo_", Population_structure, "_", number_of_pop, "pop/", folder_name, "/", folder_name, "_rep", replicate_number, "/")

# File names
Neutral_file_name <- paste0("neutral_data_g", generations, ".dat")
Quanti_file_name <- paste0("quanti_trait_g", generations, ".dat")

# Full file paths
neutral_file <- paste0(Simulation_directory, Neutral_file_name)
quanti_file <- paste0(Simulation_directory, Quanti_file_name)

# Output directory for intermediate files
intermediate_dir <- paste0("intermediate/", experiment_name, "/", folder_name, "/")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# Debug: Print paths
cat("Simulation directory:", Simulation_directory, "\n")
cat("Neutral file:", neutral_file, "\n")
cat("Quanti file:", quanti_file, "\n")
cat("Intermediate dir:", intermediate_dir, "\n")

# Read data
sim <- read.fstat(fname = neutral_file)
sim_quanti <- read.fstat(fname = quanti_file)

dos <- biall2dos(sim[,-1])
dos_quanti <- biall2dos(sim_quanti[,-1], diploid = TRUE)

# ===== SAMPLE 6 INDIVIDUALS PER POPULATION =====
pop_original <- sim[,1]
pop_sizes_original <- table(pop_original)
np <- length(pop_sizes_original)

# Set parameters for breeding design
n_individuals_per_pop <- 6
n_sire <- 3
n_dam <- 3

# Sample individuals from each population
selected_indices <- c()
for (i in 1:np) {
  pop_indices <- which(pop_original == i)
  sampled_indices <- sample(pop_indices, n_individuals_per_pop, replace = FALSE)
  selected_indices <- c(selected_indices, sampled_indices)
  
}

# Apply sampling to all data
sim <- sim[selected_indices, ]
sim_quanti <- sim_quanti[selected_indices, ]
dos <- dos[selected_indices, ]
dos_quanti <- dos_quanti[selected_indices, ]
pop <- sim[,1]

cat("Total individuals after sampling:", nrow(sim), "\n")
cat("Individuals per population:", table(pop), "\n")
# ===== END SAMPLING =====

# BREEDING PARAMETERS
offsrping_per_mating <- no <- 2
n_loci_qtl <- 100
maplength_quanti <- 50

# Split populations equally into sires and dams (3 sires, 3 dams per population)
pop_sizes <- table(sim[,1])

sire <- c()
dam <- c()

for (i in 1:np) {
  n_ind <- pop_sizes[i]
  # First half are sires, second half are dams
  n_sire_pop <- n_sire
  n_dam_pop <- n_dam

  sire <- c(sire, rep((i - 1) * n_ind + 1:(n_sire_pop), each = n_dam_pop * offsrping_per_mating))
  dam <- c(dam, rep((i - 1) * n_ind + (n_sire_pop + 1):n_ind, each = offsrping_per_mating, times = n_sire_pop))
}

# Founders as NA
nft <- sum(pop_sizes)
sire <- c(rep(NA, nft), sire)
dam <- c(rep(NA, nft), dam)
nt <- length(sire)

# Non-founders
nf <- -c(1:nft)

# Create pedigree
ped_forBreeding <- data.frame(ind = 1:nt, sire = sire, dam = dam)

# Create F1 data
F1_full_data <- create_F1(sim, sim_quanti, ped_forBreeding)

# Building population vector for F1
sire_per_pop <- as.vector(pop_sizes)/2
dam_per_pop <- as.vector(pop_sizes)/2
offspring_per_pop_F1 <- sire_per_pop * dam_per_pop * 2
pop_F1 <- rep(1:np, offspring_per_pop_F1)
pop_P_and_F1 <- c(pop, pop_F1)

# Population individual ID dataframe
population_individual_id_df <- data.frame(population = pop_F1, individual = length(pop):(length(pop_P_and_F1)-1))

# F1 data excluding founders
dos_F1only_neutral <- F1_full_data$dos_F1_neutral[nf,]

# CREATE TRAITS
trait <- c()
dosage_quanti <- F1_full_data$dos_F1_quanti[nf,] 
Y <- rowSums((dosage_quanti-1)*0.2)
err <- rnorm(length(Y), mean = 0, sd = 1)
Y <- Y + err
mean_Y <- mean(Y)
centered_Y <- Y - mean_Y
var_Y <- var(centered_Y)
trait <- c(trait, (centered_Y / sqrt(var_Y)))

# Neutral trait
trait_Neutral <- c()
dosage_quanti_Neutral <- dos_F1only_neutral[,sample(1:ncol(dos_F1only_neutral), n_loci_qtl)]
Y_neutral <- rowSums((dosage_quanti_Neutral-1)*0.2)
err_neutral <- rnorm(length(Y_neutral), mean = 0, sd = 1)
Y_neutral <- Y_neutral + err_neutral
mean_Y_neutral <- mean(Y_neutral)
centered_Y_neutral <- Y_neutral - mean_Y_neutral
var_Y_neutral <- var(centered_Y_neutral)
trait_Neutral <- c(trait_Neutral, (centered_Y_neutral / sqrt(var_Y_neutral)))

# Create trait dataframes
NTraits = 1
individual <- rep(1:(length(trait)/NTraits), NTraits)
trait_id <- rep(1:NTraits, each = (length(trait)/NTraits))
population <- rep(pop_F1, NTraits)

trait_df_pop <- data.frame(individual, trait, trait_id, population)
trait_df_pop_Neutral <- data.frame(individual, trait_Neutral, trait_id, population)

# Calculate coancestries
coancestries_dosage <- calculate_coancestries(genetic_data_parents = dos,
                                              genotyped_parent_populations = pop,
                                              genetic_data_F1 = dos_F1only_neutral, 
                                              population_individual_id = population_individual_id_df,
                                              column_individual = "individual", 
                                              column_population = "population",
                                              all_parents_genotyped = TRUE)

Theta.P <- coancestries_dosage$Theta.P
The.M <- coancestries_dosage$The.M

# Calculate FST for founders
fst.founders <- fs.dosage(dos, pop)

# Create dataf for QST-FST analysis
n_ind <- n_sire + n_dam
nft <- np * (n_sire + n_dam)
individuals <- nft:(n_ind+nft-1)
populations <- rep(1:np, each = n_sire * n_dam * no)
sires <- rep(rep(1:n_sire, each = n_dam * no), np)
dams <- rep(rep(1:n_dam, each = no), n_sire * np)

dataf <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = trait)
dataf_neutral <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = trait_Neutral)

# Prepare data for driftsel
bed <- as.bed.matrix(dos)
n_snps <- length(bed@p)
block <- 1:100
num_blocks <- floor(n_snps/100)
list_blocks <- lapply(1:num_blocks, function(i) {
  current_block <- block + 100 * (i - 1)
})

# Pedigree for Driftsel
sire_F1 <- rep(1:n_sire, each = n_dam * no)
dam_F1 <- rep(n_sire + 1:n_dam, each = no, n_sire)
sire_F1 <- rep(0:(np - 1) * (n_dam + n_sire), each = (n_dam * n_sire * no)) + sire_F1
dam_F1 <- rep(0:(np - 1) * (n_dam + n_sire), each = (n_dam * n_sire * no)) + dam_F1

ped_driftsel <- data.frame(id = (max(dam_F1) +1):((n_sire*n_dam*no*np)+max(dam_F1)), 
                          sire = sire_F1, 
                          dam = dam_F1, 
                          sire.pop = pop_F1, 
                          dam.pop = pop_F1)

covars <- data.frame(id = 1:(n_sire*n_dam*no*np), dummy = rep(1,(n_sire*n_dam*no*np)))
trait_driftsel <- data.frame(individual, trait)
trait_driftsel_neutral <- data.frame(individual, trait_Neutral)

# SAVE ALL INTERMEDIATE DATA
intermediate_file <- paste0(intermediate_dir, "prepared_data_rep", replicate_number, ".RData")

save(
  # Basic parameters
  replicate_number, Population_structure, number_of_pop, generations, folder_name,
  n_sire, n_dam, no, np, nft,
  
  # Trait data
  trait_df_pop, trait_df_pop_Neutral,
  
  # Coancestry matrices
  Theta.P, The.M,
  
  # FST data
  fst.founders,
  
  # QST-FST data
  dataf, dataf_neutral,
  
  # Driftsel data
  bed, list_blocks, ped_driftsel, covars, trait_driftsel, trait_driftsel_neutral,
  pop, pop_F1,
  
  file = intermediate_file
)

cat("Data preparation completed. Saved to:", intermediate_file, "\n")