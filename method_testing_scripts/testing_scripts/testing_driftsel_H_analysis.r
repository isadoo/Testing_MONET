source("tools/driftsel.r")
source("tools/habitat_optima.r")

library(gaston)
library(hierfstat)
library(JGTeach)

args <- commandArgs(trailingOnly = TRUE)

# Set parameters
replicate_number <- as.integer(args[1])
Population_structure <- args[2]
number_of_pop <- as.integer(args[3])
generations <- as.integer(args[4])
corr <- gsub("^_+|_+$", "", args[5])
wdiff <- args[6]
wvar <- args[7]
experiment_name <- args[8]

# Load prepared data
folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
intermediate_dir <- paste0("intermediate/", experiment_name, "/", folder_name, "/")
intermediate_file <- paste0(intermediate_dir, "prepared_data_rep", replicate_number, ".RData")

cat("Loading data from:", intermediate_file, "\n")
load(intermediate_file)

# Habitat = simulated optima, same as MONET w/ habitat; H.test drops the 1st column (pop id)
cat("Getting habitat optima for:", Population_structure, "correlation =", corr, "wdiff =", wdiff, "\n")
habitat_optima <- get_habitat_optima(Population_structure, corr, wdiff)
if (length(habitat_optima) != np) {
  stop("Number of habitat values (", length(habitat_optima), ") != number of populations (", np, ")")
}
habitat <- data.frame(pop = 1:np, habitat = habitat_optima)

# Bootstrap coancestry function
boot_coan <- function(list_blocks, bed, np, indperpop){
  sampled_blocks <- sample(list_blocks, 100, replace = TRUE)
  samp_list <- unlist(sampled_blocks)
  fst_mat <- fs.dosage(bed[,samp_list], pop=rep(1:np, each=indperpop))$FsM
  off <- fst_mat[row(fst_mat)!=col(fst_mat)]
  coan <- (fst_mat - min(off))/(1 - min(off))
}

cat("Running Driftsel H-test analysis...\n")

# Create tensor for bootstrapped coancestries
indperpop <- n_sire + n_dam
tensor.Theta.P <- replicate(25, boot_coan(list_blocks, bed, np, indperpop = indperpop))

# Run driftsel for adaptive trait
cat("Analyzing adaptive trait...\n")
samp <- MH(tensor.Theta.P, ped_driftsel, covars, trait_driftsel, 3000, 1000, 2, alt=T)
s <- neut.test(samp$pop.ef, samp$G, samp$theta,, silent=T)
h <- H.test(samp$pop.ef, samp$G, samp$theta, habitat, silent=T)

# Run driftsel for neutral trait
cat("Analyzing neutral trait...\n")
samp_neutral <- MH(tensor.Theta.P, ped_driftsel, covars, trait_driftsel_neutral, 3000, 1000, 2, alt=T)
s_neutral <- neut.test(samp_neutral$pop.ef, samp_neutral$G, samp_neutral$theta,, silent=T)
h_neutral <- H.test(samp_neutral$pop.ef, samp_neutral$G, samp_neutral$theta, habitat, silent=T)

# Save Driftsel H-test results (S-test on the same MCMC chain as _Hrun)
driftsel_results <- data.frame(
  replicate_number = replicate_number,
  folder_name = folder_name,
  H_value_Driftsel = h,
  H_value_Driftsel_Neutral = h_neutral,
  S_value_Driftsel_Hrun = mean(s),
  S_value_Driftsel_Hrun_Neutral = mean(s_neutral),
  n_nearPD_H = attr(h, "n.nearPD"),
  n_nearPD_H_Neutral = attr(h_neutral, "n.nearPD")
)

results_dir <- paste0("results/",experiment_name, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

driftsel_file <- paste0(results_dir, "driftsel_H_results_", folder_name, ".csv")

write.table(
  driftsel_results, file = driftsel_file, append = TRUE,
  sep = ",", col.names = !file.exists(driftsel_file), row.names = FALSE
)

cat("Driftsel H-test analysis completed. Results saved to:", driftsel_file, "\n")
