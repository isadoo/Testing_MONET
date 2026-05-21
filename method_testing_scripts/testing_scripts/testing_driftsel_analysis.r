source("tools/driftsel.r")

library(gaston)
library(hierfstat)
library(JGTeach)

args <- commandArgs(trailingOnly = TRUE)

# Set parameters
replicate_number <- as.integer(args[1])
Population_structure <- args[2]
number_of_pop <- as.integer(args[3])
generations <- as.integer(args[4])
correlation <- ifelse(args[5] == "NA", "", paste0("_", args[5]))
wdiff <- args[6]
wvar <- args[7]
experiment_name <- args[8]

# Load prepared data
folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
intermediate_dir <- paste0("intermediate/", experiment_name, "/", folder_name, "/")
intermediate_file <- paste0(intermediate_dir, "prepared_data_rep", replicate_number, ".RData")

cat("Loading data from:", intermediate_file, "\n")
load(intermediate_file)

# Bootstrap coancestry function
boot_coan <- function(list_blocks, bed, np, indperpop){
  sampled_blocks <- sample(list_blocks, 100, replace = TRUE)
  samp_list <- unlist(sampled_blocks)
  fst_mat <- fs.dosage(bed[,samp_list], pop=rep(1:np, each=indperpop))$FsM
  off <- fst_mat[row(fst_mat)!=col(fst_mat)]
  coan <- (fst_mat - min(off))/(1 - min(off))
}

cat("Running Driftsel analysis...\n")

# Create tensor for bootstrapped coancestries
indperpop <- n_sire + n_dam
tensor.Theta.P <- replicate(25, boot_coan(list_blocks, bed, np, indperpop = indperpop))

# Run driftsel for adaptive trait
cat("Analyzing adaptive trait...\n")
samp <- MH(tensor.Theta.P, ped_driftsel, covars, trait_driftsel, 3000, 1000, 2, alt=T)
s <- neut.test(samp$pop.ef, samp$G, samp$theta,, silent=T)
driftsel_result <- list(S_Value = mean(s))

# Run driftsel for neutral trait
cat("Analyzing neutral trait...\n")
samp_neutral <- MH(tensor.Theta.P, ped_driftsel, covars, trait_driftsel_neutral, 3000, 1000, 2, alt=T)
s_neutral <- neut.test(samp_neutral$pop.ef, samp_neutral$G, samp_neutral$theta,, silent=T)
driftsel_result_Neutral <- list(S_Value = mean(s_neutral))

# Save Driftsel results
driftsel_results <- data.frame(
  replicate_number = replicate_number,
  folder_name = folder_name,
  S_value_Driftsel = driftsel_result$S_Value,
  S_value_Driftsel_Neutral = driftsel_result_Neutral$S_Value
)

# Save results
results_dir <- paste0("results/",experiment_name, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

driftsel_file <- paste0(results_dir, "driftsel_results_", folder_name, ".csv")

write.table(
  driftsel_results, file = driftsel_file, append = TRUE, 
  sep = ",", col.names = !file.exists(driftsel_file), row.names = FALSE
)

cat("Driftsel analysis completed. Results saved to:", driftsel_file, "\n")