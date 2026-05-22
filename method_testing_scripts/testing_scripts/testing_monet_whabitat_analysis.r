source("tools/monet_whabitat.r")
source("tools/counting_blocks_matrix.r")

library(gaston)
library(hierfstat)
library(JGTeach)
library(brms)

args <- commandArgs(trailingOnly = TRUE)

# Set parameters
replicate_number <- as.integer(args[1])
Population_structure <- args[2]
number_of_pop <- as.integer(args[3])
generations <- as.integer(args[4])
correlation <- ifelse(args[5] == "NA", NA, args[5])
wdiff <- as.numeric(args[6])
wvar <- as.numeric(args[7])
experiment_name <- args[8]

#This bit here is for when I want to test the script.
# For testing:
# replicate_number <- 1
# Population_structure <- "hierarchical"
# number_of_pop <- 18
# generations <- 1000
# correlation <- "swapped"
# wdiff <- 4.6
# wvar <- 22
# experiment_name <- "quantinemo_hierarchical_18pop"

#Correct correlation parameter: remove leading/trailing underscores and potential typos
if (!is.na(correlation)) {
  #taking off underscore
  correlation <- gsub("^_+|_+$", "", correlation)
  
  #i think maybe there's this typo somewhere 
  if (correlation == "groupped") {
    correlation <- "grouped"
  
  }

}



#data loading
folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
intermediate_dir <- paste0("intermediate/", experiment_name, "/", folder_name, "/")
intermediate_file <- paste0(intermediate_dir, "prepared_data_rep", replicate_number, ".RData")

cat("Loading data from:", intermediate_file, "\n")
load(intermediate_file)

#Going from parameters to habitat info. This habitat info is what we use as a covariate in monet.
#This "habitat" is actually the optima we used in the simulations.
get_habitat_optima <- function(Population_structure, correlation) {
  
  if (Population_structure == "hierarchical") {
    if (correlation == "grouped") {
      return(rep(c(-2.3, -1.38, -0.46, 0.46, 1.38, 2.3), each = 3))
    } else if (correlation == "swapped") {
      return(rep(c(-2.3, -1.15, 0, -1.15, 0, 1.15, 0, 1.15, 2.3), 2))
    }
  } else if (Population_structure == "SS") {
    if (correlation == "cline") {
      return(c(-2.3000000, -2.0578947, -1.8157895, -1.5736842, -1.3315789, 
               -1.0894737, -0.8473684, -0.6052632, -0.3631579, -0.1210526,  
               0.1210526,  0.3631579,  0.6052632,  0.8473684, 1.0894737,  
               1.3315789,  1.5736842,  1.8157895,  2.0578947,  2.3000000))
    } else if (correlation == "sine") {
      return(c(0.1210526, 1.5736842, 2.3000000, 1.5736842, 0.1210526, 
               -1.3315789, -2.3000000, -1.3315789, 0.1210526, 1.5736842, 
               2.3000000, 1.5736842, 0.1210526, -1.3315789, -2.3000000, 
               -1.3315789, 0.1210526, 1.5736842, 2.0578947, 1.3315789))
    }
  } else if (Population_structure == "IM") {
    # For IM, correlation is NA, so we don't check it
    return(c(-2.30, -2.03, -1.76, -1.49, -1.22, -0.95, -0.68, -0.41, -0.14,  
             0.14,  0.41,  0.68,  0.95,  1.22, 1.49,  1.76,  2.03,  2.30))
  }
  
  stop(paste0("Invalid parameter combination: Population_structure=", Population_structure,
              ", correlation=", correlation))
}

#the optima is what we use as habitat information. 
cat("Getting habitat optima for:", Population_structure, "correlation =", correlation, "\n")
habitat_optima <- get_habitat_optima(Population_structure, correlation)
cat("Habitat optima:", paste(round(habitat_optima, 2), collapse = ", "), "\n")

cat("Starting MONET with habitat\n")
monet_result_whabitat <- monet_whabitat(
  Theta.P, 
  The.M, 
  trait_dataframe = trait_df_pop,
  habitat_optima = habitat_optima,
  column_individual = "individual", 
  column_trait = "trait"
)

#Similarly to before, we will run MONET for neutral trait as well - here of course, habitat should not matter.
monet_result_whabitat_Neutral <- monet_whabitat(
  Theta.P,
  The.M, 
  trait_dataframe = trait_df_pop_Neutral,
  habitat_optima = habitat_optima,
  column_individual = "individual", 
  column_trait = "trait_Neutral"
)


monet_results <- data.frame(
  replicate_number = replicate_number,
  folder_name = folder_name,
  population_structure = Population_structure,
  wdiff = wdiff,
  wvar = wvar,
  correlation = ifelse(is.na(correlation), "NA", as.character(correlation)),
  
  #Results for selected trait
  p_value_MONET = monet_result_whabitat$log_ratio$p_value,
  log_ratio_MONET = monet_result_whabitat$log_ratio$mean_log_ratio,
  log_ratio_ci_lower_MONET = monet_result_whabitat$log_ratio$log_ratio_ci_lower,
  log_ratio_ci_upper_MONET = monet_result_whabitat$log_ratio$log_ratio_ci_upper,
  
  habitat_coefficient = monet_result_whabitat$habitat_analysis$optima_coefficient_mean,
  habitat_coef_ci_lower = monet_result_whabitat$habitat_analysis$optima_coefficient_ci_lower,
  habitat_coef_ci_upper = monet_result_whabitat$habitat_analysis$optima_coefficient_ci_upper,
  habitat_p_value = monet_result_whabitat$habitat_analysis$optima_p_value,
  
  #Results for neutral trait
  p_value_MONET_Neutral = monet_result_whabitat_Neutral$log_ratio$p_value,
  log_ratio_MONET_Neutral = monet_result_whabitat_Neutral$log_ratio$mean_log_ratio,
  log_ratio_ci_lower_MONET_Neutral = monet_result_whabitat_Neutral$log_ratio$log_ratio_ci_lower,
  log_ratio_ci_upper_MONET_Neutral = monet_result_whabitat_Neutral$log_ratio$log_ratio_ci_upper,
  
  habitat_coefficient_Neutral = monet_result_whabitat_Neutral$habitat_analysis$optima_coefficient_mean,
  habitat_coef_ci_lower_Neutral = monet_result_whabitat_Neutral$habitat_analysis$optima_coefficient_ci_lower,
  habitat_coef_ci_upper_Neutral = monet_result_whabitat_Neutral$habitat_analysis$optima_coefficient_ci_upper,
  habitat_p_value_Neutral = monet_result_whabitat_Neutral$habitat_analysis$optima_p_value
)


results_dir <- paste0("results/", experiment_name, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

monet_file <- paste0(results_dir, "monet_whabitat_results_", folder_name, ".csv")

write.table(
  monet_results, 
  file = monet_file, 
  append = TRUE, 
  sep = ",", 
  col.names = !file.exists(monet_file), 
  row.names = FALSE
)

cat("MONET with habitat analysis completed. Results saved to:", monet_file, "\n")

cat("\n- Summary for Selected Trait -\n")
print(monet_result_whabitat)

cat("\n- Summary for Neutral Trait -\n")
print(monet_result_whabitat_Neutral)