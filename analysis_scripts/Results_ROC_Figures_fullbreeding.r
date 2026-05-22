## ------------------ SETUP ------------------
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(pROC)
library(ggplot2)
library(dplyr)
library(stringr)

#non-interactive graphics device for headless systems
options(bitmapType = 'cairo')

df <- read.csv("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/raw_data/3methods_full_breeding_combined_november.tsv", header = TRUE, sep = "\t")

#standardize wdiff values: convert 0.0 to 0, 10.0 to 10, keep 4p6 as is
#on the new results it's in float (SS results)
df$wdiff <- as.character(df$wdiff) 
df$wdiff[df$wdiff == "0.0"] <- "0"
df$wdiff[df$wdiff == "10.0"] <- "10"

#Split my data by population structure based on model + num_pops
im18         <- df[df$model == "IM" & df$num_pops == 18, ]
ss           <- df[df$model == "SS" & df$num_pops == 18, ]
hierarchical <- df[df$model == "hierarchical" & df$num_pops == 18, ]

#Replace correlation values "cline" and "sine" with "_cline" and "_sine"
ss$correlation <- gsub("^cline$", "_cline", ss$correlation)
ss$correlation <- gsub("^sine$", "_sine", ss$correlation)

## ------------------ LABEL/STYLING HELPERS ------------------
pretty_pop <- function(x) {
  switch(x,
         "SS"            = "Stepping stones",
         "Hierarchical"  = "Hierarchical",
         "IM_18"         = "Island Model 18",
         x)
}
disp_wdiff <- function(w) {
  w <- as.character(w)
  ifelse(w == "4p6", "4.6", w)
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
# Keep the original p-threshold grid
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

## ------------------ CORE FUNCTION ------------------
process_population_structure <- function(pop_data, pop_name) {
  
  roc_data_list <- list()
  auc_results   <- list()
  
  # Enumerate scenarios (include correlation for SS/Hierarchical)
  if (pop_name %in% c("SS", "Hierarchical")) {
    scenarios <- pop_data %>%
      filter(!is.na(wdiff) & !is.na(wvar)) %>%
      mutate(correlation = na_if(str_trim(as.character(correlation)), "")) %>%
      select(wdiff, wvar, correlation) %>%
      distinct() %>% arrange(wdiff, wvar, correlation)
  } else {
    scenarios <- pop_data %>%
      filter(!is.na(wdiff) & !is.na(wvar)) %>%
      select(wdiff, wvar) %>%
      distinct() %>% arrange(wdiff, wvar)
    scenarios$correlation <- NA
  }
  
  # Build ROC per scenario
  for (i in 1:nrow(scenarios)) {
    wdiff_val <- scenarios$wdiff[i]
    wvar_val  <- scenarios$wvar[i]
    corr_val  <- scenarios$correlation[i]
    
    if (is.na(corr_val)) {
      scenario_name <- paste0("wdiff_", wdiff_val, "_wvar_", wvar_val)
      scenario_data <- pop_data[pop_data$wdiff == wdiff_val & pop_data$wvar == wvar_val, ]
    } else {
      scenario_name <- paste0("wdiff_", wdiff_val, "_wvar_", wvar_val, "_corr_", gsub("_", "", corr_val))
      scenario_data <- pop_data[pop_data$wdiff == wdiff_val &
                                  pop_data$wvar  == wvar_val &
                                  pop_data$correlation == corr_val, ]
    }
    
    cat("Processing:", scenario_name, "\n")
    if (nrow(scenario_data) == 0) next
    
    # QSTFST
    roc_qstfst <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_QSTFST", "p_value_QSTFST_Neutral", thresh)
    }))
    roc_qstfst$method <- "QSTFST"
    roc_qstfst$scenario <- scenario_name
    roc_qstfst$wdiff <- wdiff_val
    roc_qstfst$wvar <- wvar_val
    roc_qstfst$correlation <- corr_val
    roc_qstfst <- prepare_for_log_plot(roc_qstfst)
    
    # MONET
    roc_monet <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_MONET", "p_value_MONET_Neutral", thresh)
    }))
    roc_monet$method <- "MONET"
    roc_monet$scenario <- scenario_name
    roc_monet$wdiff <- wdiff_val
    roc_monet$wvar <- wvar_val
    roc_monet$correlation <- corr_val
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
    roc_driftsel$wdiff <- wdiff_val
    roc_driftsel$wvar <- wvar_val
    roc_driftsel$correlation <- corr_val
    roc_driftsel <- prepare_for_log_plot(roc_driftsel)
    
    # AUC (uses original FPR, not FPR_plot)
    auc_qstfst   <- calculate_auc(roc_qstfst$FPR,   roc_qstfst$TPR)
    auc_monet     <- calculate_auc(roc_monet$FPR,     roc_monet$TPR)
    auc_driftsel <- calculate_auc(roc_driftsel$FPR, roc_driftsel$TPR)
    
    auc_results[[scenario_name]] <- data.frame(
      population_structure = pop_name,
      scenario = scenario_name,
      wdiff = wdiff_val,
      wvar = wvar_val,
      correlation = corr_val,
      QSTFST_AUC = auc_qstfst,
      MONET_AUC   = auc_monet,
      Driftsel_AUC = auc_driftsel
    )
    
    roc_data_list[[scenario_name]] <- rbind(roc_qstfst, roc_monet, roc_driftsel)
  }
  
  all_roc_data <- do.call(rbind, roc_data_list)
  all_roc_data$wdiff <- factor(all_roc_data$wdiff, levels = c("0", "4p6", "10"))
  all_auc_data <- do.call(rbind, auc_results)
  
  ## ------------------ PLOTTING (LOG X-AXIS) ------------------
  six.panels <- pop_name %in% c("SS", "Hierarchical")
  if (six.panels) {
    png(paste0("Newlog_", pop_name, "_log.png"),
        height = 4500, width = 3000, res = 300, pointsize = 10)
    layout(matrix(c(1, 1, 2, 7, 3, 4, 5, 6), byrow = TRUE, ncol = 2),
           heights = c(1/20, 19/60, 19/60, 19/60), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  } else {
    png(paste0("Newlog_", pop_name, "_log.png"),
        height = 3000, width = 3000, res = 300, pointsize = 10)
    layout(matrix(c(1, 1, 2, 5, 3, 4), byrow = TRUE, ncol = 2),
           heights = c(1/10, 9/20, 9/20), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  }
  
  colours   <- c("Driftsel" = "#6E0D25" , "QSTFST" = "#62929E" , "MONET" =  "#F49D37")
  linetypes <- c("10" = 1, "22" = 2, "50" = 3)
  
  for (wdiff_val in levels(all_roc_data$wdiff)) {
    
    if (pop_name %in% c("SS", "Hierarchical") && wdiff_val != "0") {
      
      correlations <- unique(all_roc_data[all_roc_data$wdiff == wdiff_val &
                                            !is.na(all_roc_data$correlation), "correlation"])
      for (corr_val in correlations) {
        
        plot_data <- all_roc_data[all_roc_data$wdiff == wdiff_val &
                                    all_roc_data$correlation == corr_val &
                                    !is.na(all_roc_data$correlation), ]
        if (nrow(plot_data) == 0) next
        first <- TRUE
        
        for (wvar_val in unique(plot_data$wvar)) {
          tmp <- plot_data[plot_data$wvar == wvar_val, ]
          for (method in unique(tmp$method)) {
            tmp.method <- tmp[tmp$method == method, ]
            
            if (first) {
              par(mar = c(5.1, 4.1, 4.1, 2.1))
              plot(tmp.method$FPR_plot, tmp.method$TPR, type = "l",
                   lty = linetypes[as.character(wvar_val)],
                   lwd = line_width,
                   col = colours[method],
                   xlim = c(1e-2, 1), ylim = c(0, 1),
                   log = "x", xaxs = "i",
                   main = paste(Delta_theta, "=", disp_wdiff(wdiff_val),
                                ", correlation =", gsub("_", "", corr_val)),
                   xlab = "FPR (log scale)", ylab = "TPR",
                   cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
              abline(v = 0.05, col = "darkgrey", lty = 2, lwd = line_width)
              first <- FALSE
            } else {
              lines(tmp.method$FPR_plot, tmp.method$TPR,
                    lty = linetypes[as.character(wvar_val)],
                    lwd = line_width,
                    col = colours[method])
            }
          }
        }
      }
      
    } else {
      
      plot_data <- all_roc_data[all_roc_data$wdiff == wdiff_val, ]
      first <- TRUE
      
      for (wvar_val in unique(plot_data$wvar)) {
        tmp <- plot_data[plot_data$wvar == wvar_val, ]
        for (method in unique(tmp$method)) {
          tmp.method <- tmp[tmp$method == method, ]
          
          if (first) {
            par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(tmp.method$FPR_plot, tmp.method$TPR, type = "l",
                 lty = linetypes[as.character(wvar_val)],
                 lwd = line_width,
                 col = colours[method],
                 xlim = c(1e-2, 1), ylim = c(0, 1),
                 log = "x", xaxs = "i",
                 main = paste(Delta_theta, "=", disp_wdiff(wdiff_val)),
                 xlab = "FPR (log scale)", ylab = "TPR",
                 cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
            abline(v = 0.05, col = "darkgrey", lty = 2, lwd = line_width)
            first <- FALSE
          } else {
            lines(tmp.method$FPR_plot, tmp.method$TPR,
                  lty = linetypes[as.character(wvar_val)],
                  lwd = line_width,
                  col = colours[method])
          }
        }
      }
    }
  }
  
  # Legend
  plot(1, 1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = c(paste0(omega_chr, " = ", c(10, 22, 50)), "Driftsel", "QSTFST", "MONET"),
         ncol   = 2,
         lty    = c(linetypes, rep(0, 3)),
         lwd    = c(rep(line_width, 3), rep(NA, 3)),
         pch    = c(rep(NA, 3), rep(15, 3)),
         col    = c(rep("black", 3), colours),
         cex = 2, bty = "n")
  
  dev.off()
  
  return(list(roc_data = all_roc_data, auc_data = all_auc_data))
}

## ------------------ RUN ------------------
# Process only IM_18, SS, and Hierarchical (exclude IM_9)
results_im18 <- process_population_structure(im18, "IM_18")
results_ss   <- process_population_structure(ss,   "SS")
results_hierarchical <- process_population_structure(hierarchical, "Hierarchical")

# Combine AUCs
all_auc_results <- rbind(
  results_im18$auc_data,
  results_ss$auc_data,
  results_hierarchical$auc_data
)

print(all_auc_results, row.names = FALSE)
write.csv(all_auc_results, "AUC_summary_table.csv", row.names = FALSE)
