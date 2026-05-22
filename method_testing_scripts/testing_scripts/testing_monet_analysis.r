source("tools/monet.r")
source("tools/counting_blocks_matrix.r")

library(gaston)
library(hierfstat)
library(JGTeach)
library(brms)

#args <- commandArgs(trailingOnly = TRUE)

Set parameters
replicate_number <- as.integer(args[1])
Population_structure <- args[2]
number_of_pop <- as.integer(args[3])
generations <- as.integer(args[4])
correlation <- ifelse(args[5] == "NA", "", paste0("_", args[5]))
wdiff <- args[6]
wvar <- args[7]
experiment_name <- args[8]

# replicate_number <- 1
# Population_structure <- "IM"
# number_of_pop <- 18
# generations <- 1000
# correlation <- "NA"
# wdiff <- 10
# wvar <- 10
# experiment_name <- "quantinemo_IM_18pop"

# Load prepared data
folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
intermediate_dir <- paste0("intermediate/", experiment_name, "/", folder_name, "/")
intermediate_file <- paste0(intermediate_dir, "prepared_data_rep", replicate_number, ".RData")

cat("Loading data from:", intermediate_file, "\n")
load(intermediate_file)

# Run MONET analysis
cat("Running MONET analysis...\n")

monet_result <- monet(Theta.P, 
                      The.M, 
                      trait_dataframe = trait_df_pop, 
                    column_individual = "individual", 
                    column_trait = "trait")

monet_result_Neutral <- monet(Theta.P,
                              The.M, 
                              trait_dataframe = trait_df_pop_Neutral, 
                              column_individual = "individual", 
                              column_trait = "trait_Neutral")

# Save MONET results
monet_results <- data.frame(
  replicate_number = replicate_number,
  folder_name = folder_name,
  p_value_MONET = monet_result$log_ratio$p_value,
  log_ratio_MONET = monet_result$log_ratio$mean_log_ratio,
  p_value_MONET_Neutral = monet_result_Neutral$log_ratio$p_value,
  log_ratio_MONET_Neutral = monet_result_Neutral$log_ratio$mean_log_ratio
)

# Save results
results_dir <- paste0("results/",experiment_name, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

monet_file <- paste0(results_dir, "monet_results_", folder_name, ".csv")

write.table(
  monet_results, file = monet_file, append = TRUE, 
  sep = ",", col.names = !file.exists(monet_file), row.names = FALSE
)

cat("MONET analysis completed. Results saved to:", monet_file, "\n")