## ------------------ SETUP ------------------
setwd("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/graphing_results")

library(pROC)
library(ggplot2)
library(dplyr)
library(stringr)
library(ragg)

#non-interactive graphics device for headless systems
options(bitmapType = 'cairo')

## ------------------ LOAD AND COMBINE DATA ------------------
# Read run1 data (TSV) - UPDATE THIS PATH TO YOUR newSS FILE
# Example: run1_wdiff4p6_newSS_combined.tsv or similar
run1 <- read.delim("../results_tables/run1_combined_newSS.tsv", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#read 3 methods
all_runs <- read.delim("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/results_tables/3methods_full_breeding_combined_november.tsv", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check if we should use different filter values or no filter
if (sum((all_runs$wdiff == 10 | all_runs$wdiff == "10.0") & (all_runs$wvar == 10 | all_runs$wvar == "10.0"), na.rm = TRUE) == 0) {
  cat("\nWARNING: No data matches wdiff=10 & wvar=10. Using all data instead.\n")
  all_runs_filtered <- all_runs
} else {
  # Filter and explicitly remove NA rows - accept both 10 and 10.0 formats
  all_runs_filtered <- all_runs[!is.na(all_runs$wdiff) & !is.na(all_runs$wvar) & 
                                 (all_runs$wdiff == 10 | all_runs$wdiff == "10.0") & 
                                 (all_runs$wvar == 10 | all_runs$wvar == "10.0"), ]
}

# Add pops_analyzed column
run1$pops_analyzed <- 9
all_runs_filtered$pops_analyzed <- all_runs_filtered$num_pops
run1$experiment_name <- NULL
# Combine datasets
df <- rbind(run1, all_runs_filtered)

#standardize wdiff values: convert 0.0 to 0, 10.0 to 10, keep 4p6 as is
#on the new results it's in float (SS results)
df$wdiff <- as.character(df$wdiff) 
df$wdiff[df$wdiff == "0.0"] <- "0"
df$wdiff[df$wdiff == "10.0"] <- "10"

#Replace correlation values "cline" and "sine" with "_cline" and "_sine" for SS model
# Also replace "grouped" and "swapped" with "_grouped" and "_swapped" for hierarchical model
# This ensures consistency between datasets
df$correlation[df$model == "SS"] <- gsub("^cline$", "_cline", df$correlation[df$model == "SS"])
df$correlation[df$model == "SS"] <- gsub("^sine$", "_sine", df$correlation[df$model == "SS"])
df$correlation[df$model == "hierarchical"] <- gsub("^grouped$", "_grouped", df$correlation[df$model == "hierarchical"])
df$correlation[df$model == "hierarchical"] <- gsub("^swapped$", "_swapped", df$correlation[df$model == "hierarchical"])

#To see what we got (sanity check)
for (model_type in c("IM", "SS", "hierarchical")) {
  cat("\nModel:", model_type, "\n")
  subset_df <- df[!is.na(df$model) & df$model == model_type, ]
  if (nrow(subset_df) > 0) {
    print(table(subset_df$correlation, subset_df$pops_analyzed, useNA = "always"))
  }
}

Delta_theta <- "\u0394\u03B8"   # Δθ
omega_chr <- "ωₒ"           # ω

# All lines the same thickness
line_width <- 5
# Text sizes (multiplier on device pointsize)
axis_cex   <- 1.3
label_cex  <- 1.5
main_cex   <- 2.35

## ------------------ THRESHOLDS ------------------
# Keep the original p-threshold grid - more granular as in Results_ROC_Figures_november.r
p_thresholds <- seq(0.0001, 1, by = 0.001)

# Map two-tailed p-values to folded-S thresholds via Normal on |S-0.5|
s_sd <- 0.10
S_fold_from_p <- function(p, s_sd = 0.10) {
  k  <- qnorm(1 - p/2, mean = 0, sd = s_sd)   # distance in |S-0.5|
  th <- 0.5 - k                               # folded S threshold
  pmin(pmax(th, 0), 0.5)
}
s_fold_thresholds <- S_fold_from_p(p_thresholds, s_sd = s_sd)

## ------------------ RATE CALCULATORS ------------------
calculate_rates_pvalue <- function(data, p_col, p_neutral_col, threshold) {
  TP <- sum(data[[p_col]] < threshold, na.rm = TRUE)
  FP <- sum(data[[p_neutral_col]] < threshold, na.rm = TRUE)
  total_selective <- sum(!is.na(data[[p_col]]))
  total_neutral   <- sum(!is.na(data[[p_neutral_col]]))
  TPR <- TP / total_selective
  FPR <- FP / total_neutral
  data.frame(TPR = TPR, FPR = FPR, threshold = threshold)
}

# folded S = min(S, 1 - S); select if S_fold < threshold
calculate_rates_svalue_folded <- function(data, s_col, s_neutral_col, s_fold_threshold) {
  S_sel_fold <- pmin(data[[s_col]], 1 - data[[s_col]])
  S_neu_fold <- pmin(data[[s_neutral_col]], 1 - data[[s_neutral_col]])
  TP <- sum(S_sel_fold < s_fold_threshold, na.rm = TRUE)
  FP <- sum(S_neu_fold < s_fold_threshold, na.rm = TRUE)
  total_selective <- sum(!is.na(S_sel_fold))
  total_neutral   <- sum(!is.na(S_neu_fold))
  TPR <- TP / total_selective
  FPR <- FP / total_neutral
  data.frame(TPR = TPR, FPR = FPR, threshold = s_fold_threshold)
}

# AUC (trapezoid)
calculate_auc <- function(fpr, tpr) {
  ord <- order(fpr); fpr <- fpr[ord]; tpr <- tpr[ord]
  sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)])) / 2
}

# Prepare for log-x plotting (keep FPR for AUC; add FPR_plot for display)
prepare_for_log_plot <- function(roc_df, eps = 1e-6) {
  if (!nrow(roc_df)) return(roc_df)
  roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]
  roc_df$FPR_plot <- pmax(roc_df$FPR, eps)
  roc_df
}

## ------------------ PROCESS SCENARIOS ------------------
process_scenario <- function(scenario_data, scenario_name) {
  
  if (nrow(scenario_data) == 0) {
    cat("No data for:", scenario_name, "\n")
    return(NULL)
  }
  
  cat("Processing:", scenario_name, "\n")
  
  # QSTFST
  roc_qstfst <- do.call(rbind, lapply(p_thresholds, function(thresh) {
    calculate_rates_pvalue(scenario_data, "p_value_QSTFST", "p_value_QSTFST_Neutral", thresh)
  }))
  roc_qstfst$method <- "QSTFST"
  roc_qstfst$scenario <- scenario_name
  roc_qstfst <- prepare_for_log_plot(roc_qstfst)
  
  # MONET
  roc_monet <- do.call(rbind, lapply(p_thresholds, function(thresh) {
    calculate_rates_pvalue(scenario_data, "p_value_MONET", "p_value_MONET_Neutral", thresh)
  }))
  roc_monet$method <- "MONET"
  roc_monet$scenario <- scenario_name
  roc_monet <- prepare_for_log_plot(roc_monet)
  
  # Driftsel on folded S using thresholds mapped from p via Normal
  roc_driftsel <- do.call(rbind, lapply(s_fold_thresholds, function(th_sfold) {
    calculate_rates_svalue_folded(
      scenario_data,
      "S_value_Driftsel", "S_value_Driftsel_Neutral",
      th_sfold
    )
  }))
  roc_driftsel$method <- "Driftsel"
  roc_driftsel$scenario <- scenario_name
  roc_driftsel <- prepare_for_log_plot(roc_driftsel)
  
  # AUC (uses original FPR, not FPR_plot)
  auc_qstfst   <- calculate_auc(roc_qstfst$FPR, roc_qstfst$TPR)
  auc_monet     <- calculate_auc(roc_monet$FPR, roc_monet$TPR)
  auc_driftsel <- calculate_auc(roc_driftsel$FPR, roc_driftsel$TPR)
  
  auc_data <- data.frame(
    scenario = scenario_name,
    QSTFST_AUC = auc_qstfst,
    MONET_AUC = auc_monet,
    Driftsel_AUC = auc_driftsel
  )
  
  roc_data <- rbind(roc_qstfst, roc_monet, roc_driftsel)
  
  return(list(roc_data = roc_data, auc_data = auc_data))
}

## ------------------ DEFINE SCENARIOS ------------------
scenarios <- list(
  list(name = "Island Model", model = "IM", correlation = NA, pops_analyzed = c(9, 18)),
  list(name = "SS - Cline", model = "SS", correlation = "_cline", pops_analyzed = c(9, 18)),
  list(name = "SS - Sine", model = "SS", correlation = "_sine", pops_analyzed = c(9, 18)),
  list(name = "Hierarchical - Grouped", model = "hierarchical", correlation = "_grouped", pops_analyzed = c(9, 18)),
  list(name = "Hierarchical - Swapped", model = "hierarchical", correlation = "_swapped", pops_analyzed = c(9, 18))
)

## ------------------ PROCESS ALL SCENARIOS ------------------
all_results <- list()
all_auc <- list()

for (i in 1:length(scenarios)) {
  scenario <- scenarios[[i]]
  
  for (np in scenario$pops_analyzed) {
    if (is.na(scenario$correlation)) {
      # Island Model - correlation is NA
      scenario_data <- df[df$model == scenario$model & 
                          df$pops_analyzed == np &
                          is.na(df$correlation), ]
      scenario_label <- paste0(scenario$name, " (", np, " pops)")
    } else {
      # Other models - correlation values should now be standardized with underscore prefix
      scenario_data <- df[df$model == scenario$model & 
                          df$pops_analyzed == np &
                          !is.na(df$correlation) &
                          df$correlation == scenario$correlation, ]
      scenario_label <- paste0(scenario$name, " (", np, " pops)")
    }
    
    result <- process_scenario(scenario_data, scenario_label)
    
    if (!is.null(result)) {
      # Add metadata - this identifies which dataset (9 or 18 pops) this ROC curve came from
      result$roc_data$graph_id <- i
      result$roc_data$pops_analyzed <- np  
      result$roc_data$graph_name <- scenario$name
      result$auc_data$graph_id <- i
      result$auc_data$pops_analyzed <- np
      result$auc_data$graph_name <- scenario$name
      
      all_results[[scenario_label]] <- result$roc_data
      all_auc[[scenario_label]] <- result$auc_data
    }
  }
}

# Combine all results
all_roc_data <- do.call(rbind, all_results)
all_auc_data <- do.call(rbind, all_auc)


for (gid in unique(all_roc_data$graph_id)) {
  cat("\nGraph ID", gid, "-", unique(all_roc_data$graph_name[all_roc_data$graph_id == gid]), ":\n")
  subset_data <- all_roc_data[all_roc_data$graph_id == gid, ]
  print(table(subset_data$method, subset_data$pops_analyzed))
}
cat("\n")

## ------------------ PLOTTING ------------------
agg_png("results_December2025/ROC_PopNumber_Comparison_newSS.png",
    height = 3000, width = 3600, res = 300, pointsize = 10)

# Colors and line types
colours   <- c("Driftsel" = "#6E0D25" , "QSTFST" = "#62929E" , "MONET" =  "#F49D37")
linetypes <- c("9" = 2, "18" = 1)
# Line widths: decreasing thickness so overlapping lines are visible
linewidths <- c("Driftsel" = 6, "QSTFST" = 4, "MONET" = 2)

# Create layout: 1 row for title, 2 rows x 3 columns, last spot for legend
# Layout: panel 1 = title (spans top row)
#         panels 2-6: 2=IM, 3=SS-Cline, 4=SS-Sine, 5=Hier-Grouped, 6=Hier-Swapped
layout(matrix(c(1, 1, 1,
                2, 3, 4,
                5, 6, 7), byrow = TRUE, ncol = 3),
       heights = c(0.08, 0.46, 0.46))

# Title
par(mar = c(0, 0, 1, 0))
plot.new()
text(0.5, 0.5, "Breeding tests: Number of subpopulations", cex = 3.5, font = 2)

# Plot each graph (5 graphs + 1 legend)
graph_counter <- 0
for (graph_id in 1:5) {
  graph_counter <- graph_counter + 1
  plot_data <- all_roc_data[all_roc_data$graph_id == graph_id, ]
  
  if (nrow(plot_data) == 0) {
    cat("Warning: No data for graph_id", graph_id, "\n")
    # Create empty plot to maintain layout
    par(mar = c(5, 4.5, 3, 2))
    plot.new()
    text(0.5, 0.5, paste("No data for graph", graph_id), cex = 2, col = "red")
    next
  }
  
  graph_name <- unique(plot_data$graph_name)[1]
  cat("Plotting:", graph_name, "(panel", graph_counter + 1, ")\n")
  cat("  Available pops:", unique(plot_data$pops_analyzed), "\n")
  cat("  Rows per pops: 9pops=", nrow(plot_data[plot_data$pops_analyzed == 9,]),
      "  18pops=", nrow(plot_data[plot_data$pops_analyzed == 18,]), "\n")
  
  # Set margins for this panel
  par(mar = c(5, 4.5, 3, 2))
  
  # Create empty plot with proper axes - use smaller xlim to show low FPR values
  plot.new()
  plot.window(xlim = c(1e-2, 1), ylim = c(0, 1.05), log = "x", xaxs = "i", yaxs = "i")
  
  # Add axes
  axis(1, cex.axis = axis_cex)
  axis(2, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2), 
       cex.axis = axis_cex, las = 1)
  
  # Add box around plot
  box()
  
  # Add labels and title
  title(main = graph_name, xlab = "FPR (log scale)", ylab = "TPR",
        cex.lab = label_cex, cex.main = main_cex)
  
  # Add the reference line at FPR = 0.05
  abline(v = 0.05, col = "darkgrey", lty = 2, lwd = line_width)
  
  # Now add all the lines - NO CLIPPING, use all data
  for (np in c(9, 18)) {
    tmp_np <- plot_data[plot_data$pops_analyzed == np, ]
    
    if (nrow(tmp_np) == 0) {
      cat("  Warning: No data for", np, "pops\n")
      next
    }
    
    for (method in c("Driftsel", "QSTFST", "MONET")) {
      tmp_method <- tmp_np[tmp_np$method == method, ]
      
      if (nrow(tmp_method) == 0) {
        cat("  Warning: No data for", method, "with", np, "pops\n")
        next
      }
      
      # Debug: check for degenerate data
      tpr_range <- range(tmp_method$TPR, na.rm = TRUE)
      fpr_range <- range(tmp_method$FPR_plot, na.rm = TRUE)
      
      cat("  Plotting:", method, "with", np, "pops (", nrow(tmp_method), "points)",
          "TPR range:", round(tpr_range[1], 3), "-", round(tpr_range[2], 3),
          "FPR range:", sprintf("%.2e", fpr_range[1]), "-", sprintf("%.2e", fpr_range[2]), "\n")
      
      if (nrow(tmp_method) > 0) {
        lines(tmp_method$FPR_plot, tmp_method$TPR,
              lty = linetypes[as.character(np)],
              lwd = linewidths[method],  # Use method-specific line width
              col = colours[method])
      }
    }
  }
}

# Legend in panel 7 (bottom right)
par(mar = c(2, 2, 2, 2))
plot.new()
legend("center",
       legend = c("9 pops", "18 pops", "", "Driftsel", "QSTFST", "MONET"),
       lty = c(linetypes, 0, rep(1, 3)),
       lwd = c(rep(line_width, 2), NA, rep(line_width, 3)),
       col = c(rep("black", 2), NA, colours),
       cex = 2.2, bty = "n", seg.len = 3)

dev.off()

# Save AUC results
print(all_auc_data, row.names = FALSE)
write.csv(all_auc_data, "results_December2025/AUC_PopStructure_summary_newSS.csv", row.names = FALSE)

