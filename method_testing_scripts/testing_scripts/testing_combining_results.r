
# Set parameters  
replicate_number <- as.integer(args[1])
Population_structure <- args[2]
number_of_pop <- as.integer(args[3])
generations <- as.integer(args[4])
correlation <- ifelse(args[5] == "NA", "", paste0("_", args[5]))
wdiff <- args[6]
wvar <- args[7]
experiment_name <- args[8]

folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
results_dir <- paste0("results/",experiment_name, "/")

# Read individual results files
monet_file <- paste0(results_dir, "monet_results_", folder_name, ".csv")
qstfst_file <- paste0(results_dir, "qstfst_results_", folder_name, ".csv")
driftsel_file <- paste0(results_dir, "driftsel_results_", folder_name, ".csv")

# Check if all files exist
if (!file.exists(monet_file) || !file.exists(qstfst_file) || !file.exists(driftsel_file)) {
  cat("Warning: Not all result files exist for replicate", replicate_number, "folder", folder_name, "\n")
  cat("MONET file exists:", file.exists(monet_file), "\n")
  cat("QSTFST file exists:", file.exists(qstfst_file), "\n")
  cat("Driftsel file exists:", file.exists(driftsel_file), "\n")
  quit(status = 1)
}

# Read results
monet_results <- read.csv(monet_file)
qstfst_results <- read.csv(qstfst_file)
driftsel_results <- read.csv(driftsel_file)

# Filter for current replicate
monet_rep <- monet_results[monet_results$replicate_number == replicate_number, ]
qstfst_rep <- qstfst_results[qstfst_results$replicate_number == replicate_number, ]
driftsel_rep <- driftsel_results[driftsel_results$replicate_number == replicate_number, ]

# Check if all methods have results for this replicate
if (nrow(monet_rep) == 0 || nrow(qstfst_rep) == 0 || nrow(driftsel_rep) == 0) {
  cat("Warning: Missing results for replicate", replicate_number, "in folder", folder_name, "\n")
  quit(status = 1)
}

# Combine results
combined_results <- data.frame(
  replicate_number = replicate_number,
  p_value_QSTFST = qstfst_rep$p_value_QSTFST,
  p_value_MONET = monet_rep$p_value_MONET,
  log_ratio_MONET = monet_rep$log_ratio_MONET,
  S_value_Driftsel = driftsel_rep$S_value_Driftsel,
  p_value_QSTFST_Neutral = qstfst_rep$p_value_QSTFST_Neutral,
  p_value_MONET_Neutral = monet_rep$p_value_MONET_Neutral,
  log_ratio_MONET_Neutral = monet_rep$log_ratio_MONET_Neutral,
  S_value_Driftsel_Neutral = driftsel_rep$S_value_Driftsel_Neutral
)

# Save to final results file
final_results_file <- paste0(results_dir, "ALLresults_", Population_structure, "_", number_of_pop, "pop_", folder_name, ".csv")

write.table(
  combined_results, file = final_results_file, append = TRUE, 
  sep = ",", col.names = !file.exists(final_results_file), row.names = FALSE
)

cat("Combined results saved to:", final_results_file, "\n")