## ------------------ SETUP ------------------

library(dplyr)
library(stringr)

# Read data
df <- read.csv("../raw_data/3methods_full_breeding_combined_november.tsv", header = TRUE, sep = "\t")

## ------------------ THRESHOLDS ------------------
# Fine-grained p-thresholds for all methods to maximize chance of hitting FPR = 0.05
p_thresholds <- c(
  seq(0.0001, 0.001, by = 0.00005),   # Very fine for low p-values
  seq(0.001, 0.01, by = 0.0005),       # Fine for medium-low
  seq(0.01, 0.1, by = 0.001),          # Medium granularity
  seq(0.1, 1, by = 0.005)              # Coarser for high values
)
p_thresholds <- sort(unique(p_thresholds))

# For Driftsel: map p-values to folded S thresholds via Normal distribution
s_sd <- 0.10
S_fold_from_p <- function(p, s_sd = 0.10) {
  k  <- qnorm(1 - p/2, mean = 0, sd = s_sd)   # distance in |S-0.5|
  th <- 0.5 - k                               # folded S threshold
  pmin(pmax(th, 0), 0.5)
}
s_fold_thresholds <- S_fold_from_p(p_thresholds, s_sd = s_sd)

## ------------------ RATE CALCULATORS ------------------
# Calculate FPR and TPR for p-value methods (QSTFST, MONET)
calculate_rates_pvalue <- function(data, p_col, p_neutral_col, threshold) {
  TP <- sum(data[[p_col]] < threshold, na.rm = TRUE)
  FP <- sum(data[[p_neutral_col]] < threshold, na.rm = TRUE)
  total_selective <- sum(!is.na(data[[p_col]]))
  total_neutral   <- sum(!is.na(data[[p_neutral_col]]))
  
  if (total_selective == 0 || total_neutral == 0) {
    return(data.frame(TPR = NA, FPR = NA, threshold = threshold))
  }
  
  TPR <- TP / total_selective
  FPR <- FP / total_neutral
  data.frame(TPR = TPR, FPR = FPR, threshold = threshold)
}

# Calculate FPR and TPR for Driftsel (folded S-value)
calculate_rates_svalue_folded <- function(data, s_col, s_neutral_col, s_fold_threshold) {
  S_sel_fold <- pmin(data[[s_col]], 1 - data[[s_col]])
  S_neu_fold <- pmin(data[[s_neutral_col]], 1 - data[[s_neutral_col]])
  
  TP <- sum(S_sel_fold < s_fold_threshold, na.rm = TRUE)
  FP <- sum(S_neu_fold < s_fold_threshold, na.rm = TRUE)
  total_selective <- sum(!is.na(S_sel_fold))
  total_neutral   <- sum(!is.na(S_neu_fold))
  
  if (total_selective == 0 || total_neutral == 0) {
    return(data.frame(TPR = NA, FPR = NA, threshold = s_fold_threshold))
  }
  
  TPR <- TP / total_selective
  FPR <- FP / total_neutral
  data.frame(TPR = TPR, FPR = FPR, threshold = 2* s_fold_threshold)
}

## ------------------ INTERPOLATION FUNCTION ------------------
# Interpolate TPR at FPR = 0.05
# If exact FPR = 0.05 exists, use it; otherwise interpolate between surrounding points
get_tpr_at_fpr_005 <- function(roc_df, target_fpr = 0.05) {
  # Remove NA values and sort by FPR
  roc_df <- roc_df[!is.na(roc_df$FPR) & !is.na(roc_df$TPR), ]
  if (nrow(roc_df) == 0) return(NA)
  
  roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]
  
  # Check if exact match exists
  exact_match <- roc_df[abs(roc_df$FPR - target_fpr) < 1e-10, ]
  if (nrow(exact_match) > 0) {
    return(mean(exact_match$TPR))  # Average if multiple exact matches
  }
  
  # Find bracketing points
  below <- roc_df[roc_df$FPR <= target_fpr, ]
  above <- roc_df[roc_df$FPR >= target_fpr, ]
  
  if (nrow(below) == 0 || nrow(above) == 0) {
    # Cannot interpolate - target_fpr is outside the range
    return(NA)
  }
  
  # Get closest points for interpolation
  fpr1 <- max(below$FPR)
  fpr2 <- min(above$FPR)
  tpr1 <- below$TPR[below$FPR == fpr1][1]  # Take first if multiple
  tpr2 <- above$TPR[above$FPR == fpr2][1]
  
  # Linear interpolation
  if (fpr2 == fpr1) {
    return(mean(c(tpr1, tpr2)))
  }
  
  tpr_at_005 <- tpr1 + (tpr2 - tpr1) * (target_fpr - fpr1) / (fpr2 - fpr1)
  return(tpr_at_005)
}

## ------------------ PROCESS ALL SCENARIOS ------------------
process_all_scenarios <- function(data) {
  
  # Filter to only num_pops == 18 and exclude IM with 9 populations
  # Keep only: IM_18, SS, hierarchical
  data <- data %>%
    filter(num_pops == 18, model %in% c("IM", "SS", "hierarchical"))
  
  #Normalize the data: treat NA and "" as equivalent
  data <- data %>%
    mutate(
      fst_clean = ifelse(is.na(fst) | fst == "", "none", as.character(fst)),
      correlation_clean = ifelse(is.na(correlation) | correlation == "", "none", as.character(correlation)),
      wdiff_clean = as.character(wdiff),
      wvar_clean = as.character(wvar)
    )
  
  #Identify all unique scenarios by model, wdiff, wvar, correlation
  scenarios <- data %>%
    filter(!is.na(wdiff) & !is.na(wvar)) %>%
    select(model, wdiff_clean, wvar_clean, correlation_clean) %>%
    distinct() %>%
    arrange(model, wdiff_clean, wvar_clean, correlation_clean)
  
  cat("Found", nrow(scenarios), "unique scenarios (num_pops == 18 only)\n")
  
  results_list <- list()
  
  for (i in 1:nrow(scenarios)) {
    scenario <- scenarios[i, ]
    
    cat("Processing scenario", i, "of", nrow(scenarios), ":", 
        scenario$model, "wdiff=", scenario$wdiff_clean, 
        "wvar=", scenario$wvar_clean, "corr=", scenario$correlation_clean, "\n")
    
    
    scenario_data <- data %>%
      filter(
        model == scenario$model,
        wdiff_clean == scenario$wdiff_clean,
        wvar_clean == scenario$wvar_clean,
        correlation_clean == scenario$correlation_clean
      )
    
    if (nrow(scenario_data) == 0) {
      cat("  Warning: No data found for this scenario\n")
      next
    }
    
    cat("  Found", nrow(scenario_data), "replicates\n")
    
    ## --- QSTFST ---
    roc_qstfst <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_QSTFST", "p_value_QSTFST_Neutral", thresh)
    }))
    tpr_qstfst <- get_tpr_at_fpr_005(roc_qstfst)
    
    ## --- MONET ---
    roc_monet <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_MONET", "p_value_MONET_Neutral", thresh)
    }))
    tpr_monet <- get_tpr_at_fpr_005(roc_monet)
    
    ## --- Driftsel ---
    roc_driftsel <- do.call(rbind, lapply(s_fold_thresholds, function(th_sfold) {
      calculate_rates_svalue_folded(
        scenario_data,
        "S_value_Driftsel", "S_value_Driftsel_Neutral",
        th_sfold
      )
    }))
    tpr_driftsel <- get_tpr_at_fpr_005(roc_driftsel)
    
    ## --- Store results ---
    #Create a unique key for the list
    list_key <- paste(scenario$model, scenario$wdiff_clean, 
                      scenario$wvar_clean, scenario$correlation_clean, sep = "_")
    
    results_list[[list_key]] <- data.frame(
      model = scenario$model,
      wdiff = scenario$wdiff_clean,
      wvar = scenario$wvar_clean,
      correlation = ifelse(scenario$correlation_clean == "none", NA, scenario$correlation_clean),
      TPR_at_FPR005_QSTFST = tpr_qstfst,
      TPR_at_FPR005_MONET = tpr_monet,
      TPR_at_FPR005_Driftsel = tpr_driftsel,
      n_replicates = nrow(scenario_data),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  return(results_df)
}

## ------------------ RUN 

results <- process_all_scenarios(df)

cat("\nResults summary:\n")
print(results)

# Save results
output_file <- "../results/TPR_at_FPR_005_all_scenarios_november2025.csv"
write.csv(results, output_file, row.names = FALSE)

## ------------------ CREATE LATEX TABLE ------------------
library(xtable)

cat("\n=== CREATING LATEX TABLE ===\n")


results_latex <- results %>%
  mutate(
    `Population structure` = case_when(
      model == "IM" ~ "Island Model",
      model == "SS" ~ "Stepping Stones",
      model == "hierarchical" ~ "Hierarchical",
      TRUE ~ model
    ),
    `$\\Delta\\theta$` = case_when(
      wdiff == "4p6" ~ "4.6",
      TRUE ~ wdiff
    ),
    `$\\omega_o$` = wvar,
    Correlation = ifelse(is.na(correlation), "", gsub("_", "", correlation)),
    QSTFST = round(TPR_at_FPR005_QSTFST, 3),
    MONET = round(TPR_at_FPR005_MONET, 3),
    Driftsel = round(TPR_at_FPR005_Driftsel, 3)
  ) %>%
  select(`Population structure`, `$\\Delta\\theta$`, `$\\omega_o$`, 
         Correlation, QSTFST, MONET, Driftsel)


latex_table <- xtable(results_latex,
                      caption = "True Positive Rate (TPR) at False Positive Rate (FPR) = 0.05 for different population structures and parameter combinations. Results are based on 500 replicates per scenario.",
                      label = "tab:tpr_fpr005",
                      align = c("l", "l", "r", "r", "l", "r", "r", "r"))


output_latex <- "../results/TPR_at_FPR_005_table_november2025.tex"
print(latex_table,
      file = output_latex,
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity)

cat("LaTeX table saved to '", output_latex, "'\n", sep = "")


cat("\nLaTeX code:\n")
print(latex_table,
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      sanitize.colnames.function = identity)

