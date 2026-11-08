## ------------------ SETUP ------------------
setwd("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/graphing_results")

library(pROC)
library(ggplot2)
library(dplyr)
library(stringr)
options(bitmapType = "cairo")

df <- read.csv("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/raw_data/full_breeding_combined_review2026.tsv", header = TRUE, sep = "\t")
df$wdiff <- as.character(df$wdiff)

#Split by population structure based on model + num_pops
im18         <- df[df$model == "IM" & df$num_pops == 18, ]
im9          <- df[df$model == "IM" & df$num_pops == 9, ]
ss           <- df[df$model == "SS" & df$num_pops == 18, ]
hierarchical <- df[df$model == "hierarchical" & df$num_pops == 18, ]

## ------------------ LABELs ------------------
pretty_pop <- function(x) {
  switch(x,
         "SS"            = "Stepping stones",
         "Hierarchical"  = "Hierarchical",
         "IM_9"          = "Island Model 9",
         "IM_18"         = "Island Model 18",
         x)
}
disp_wdiff <- function(w) {
  w <- as.character(w)
  ifelse(w == "4p6", "4.6", w)
}

Delta_theta <- "\u0394\u03B8"         # Δθ
omega_chr  <- "ωₒ"                    # ωo using subscript o


line_width <- 3
#Text sizes 
axis_cex   <- 1.2
label_cex  <- 1.25
main_cex   <- 1.35

## ------------------ THRESHOLDS ------------------
#p-threshold grid 
p_thresholds <- sort(unique(c(10^seq(-7, -4, length.out = 1000), seq(0.0001, 1, by = 0.001))))

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

#AUC (trapezoid) — still computed and saved, though not used in the new plot
calculate_auc <- function(fpr, tpr) {
  ord <- order(fpr); fpr <- fpr[ord]; tpr <- tpr[ord]
  sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)])) / 2
}

#Prepare tables for plotting: sort by FPR (increasing)
prepare_for_linear_plot <- function(roc_df) {
  if (!nrow(roc_df)) return(roc_df)
  roc_df[order(roc_df$FPR, roc_df$threshold), ]
}

## ------------------------------------
#This version draws Threshold (y) vs FPR (x), one line per method × ωₒ
process_population_structure <- function(pop_data, pop_name) {
  
  roc_data_list <- list()
  auc_results   <- list()
  
  #(include correlation for SS/Hierarchical)
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
  
  #Build ROC-like tables per scenario (we use threshold & FPR columns)
  for (i in 1:nrow(scenarios)) {
    wdiff_val <- scenarios$wdiff[i]
    wvar_val  <- scenarios$wvar[i]
    corr_val  <- scenarios$correlation[i]
    
    if (is.na(corr_val)) {
      scenario_name <- paste0("wdiff_", wdiff_val, "_wvar_", wvar_val)
      scenario_data <- pop_data[which(pop_data$wdiff == wdiff_val & pop_data$wvar == wvar_val), ]
    } else {
      scenario_name <- paste0("wdiff_", wdiff_val, "_wvar_", wvar_val, "_corr_", gsub("_", "", corr_val))
      scenario_data <- pop_data[which(pop_data$wdiff == wdiff_val &
                                  pop_data$wvar  == wvar_val &
                                  pop_data$correlation == corr_val), ]
    }
    
    cat("Processing:", scenario_name, "\n")
    if (nrow(scenario_data) == 0) next
    
    #QSTFST
    roc_qstfst <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_QSTFST", "p_value_QSTFST_Neutral", thresh)
    }))
    roc_qstfst$method <- "QSTFST"
    roc_qstfst$scenario <- scenario_name
    roc_qstfst$wdiff <- wdiff_val
    roc_qstfst$wvar <- wvar_val
    roc_qstfst$correlation <- corr_val
    roc_qstfst <- prepare_for_linear_plot(roc_qstfst)
    
    #MONET
    roc_MONET <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_MONET", "p_value_MONET_Neutral", thresh)
    }))
    roc_MONET$method <- "MONET"
    roc_MONET$scenario <- scenario_name
    roc_MONET$wdiff <- wdiff_val
    roc_MONET$wvar <- wvar_val
    roc_MONET$correlation <- corr_val
    roc_MONET <- prepare_for_linear_plot(roc_MONET)
    
    #Driftsel (folded S thresholds; y-axis = 2 * folded S so it spans 0-1)
    roc_driftsel <- do.call(rbind, lapply(s_fold_thresholds, function(th_sfold) {
      calculate_rates_svalue_folded(
        scenario_data,
        "S_value_Driftsel", "S_value_Driftsel_Neutral",
        th_sfold
      )
    }))
    roc_driftsel$threshold <- 2 * roc_driftsel$threshold
    roc_driftsel$method <- "Driftsel"
    roc_driftsel$scenario <- scenario_name
    roc_driftsel$wdiff <- wdiff_val
    roc_driftsel$wvar <- wvar_val
    roc_driftsel$correlation <- corr_val
    roc_driftsel <- prepare_for_linear_plot(roc_driftsel)
    
    roc_htest <- do.call(rbind, lapply(s_fold_thresholds, function(th_sfold) {
      calculate_rates_svalue_folded(
        scenario_data,
        "H_value_Driftsel", "H_value_Driftsel_Neutral",
        th_sfold
      )
    }))
    roc_htest$threshold <- 2 * roc_htest$threshold
    roc_htest$method <- "H-test"
    roc_htest$scenario <- scenario_name
    roc_htest$wdiff <- wdiff_val
    roc_htest$wvar <- wvar_val
    roc_htest$correlation <- corr_val
    roc_htest <- prepare_for_linear_plot(roc_htest)
    
    #MONET w/ habitat
    roc_habitat <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "habitat_p_value", "habitat_p_value_Neutral", thresh)
    }))
    roc_habitat$method <- "MONET w/ environment"
    roc_habitat$scenario <- scenario_name
    roc_habitat$wdiff <- wdiff_val
    roc_habitat$wvar <- wvar_val
    roc_habitat$correlation <- corr_val
    roc_habitat <- prepare_for_linear_plot(roc_habitat)
    
    #AUC (still computed for your table)
    auc_qstfst   <- calculate_auc(roc_qstfst$FPR,   roc_qstfst$TPR)
    auc_MONET     <- calculate_auc(roc_MONET$FPR,     roc_MONET$TPR)
    auc_driftsel <- calculate_auc(roc_driftsel$FPR, roc_driftsel$TPR)
    auc_habitat  <- calculate_auc(roc_habitat$FPR,  roc_habitat$TPR)
    auc_htest    <- calculate_auc(roc_htest$FPR,    roc_htest$TPR)
    
    auc_results[[scenario_name]] <- data.frame(
      population_structure = pop_name,
      scenario = scenario_name,
      wdiff = wdiff_val,
      wvar = wvar_val,
      correlation = corr_val,
      QSTFST_AUC = auc_qstfst,
      MONET_AUC   = auc_MONET,
      Driftsel_AUC = auc_driftsel,
      Habitat_AUC = auc_habitat,
      Htest_AUC = auc_htest
    )
    
    roc_data_list[[scenario_name]] <- rbind(roc_qstfst, roc_MONET, roc_driftsel, roc_habitat, roc_htest)
  }
  
  all_roc_data <- do.call(rbind, roc_data_list)
  all_roc_data$wdiff <- factor(all_roc_data$wdiff, levels = c("0", "4p6", "10"))
  all_auc_data <- do.call(rbind, auc_results)
  
  ## ------------------ Threshold vs FPR ------------------------------------
  six.panels <- pop_name %in% c("SS", "Hierarchical")
  
  # Generate PNG version
  if (six.panels) {
    png(paste0("ThreshVsFPR_", pop_name, "_linear_environment.png"),
        height = 4500, width = 3000, res = 300, pointsize = 10)
    layout(matrix(c(1, 1, 2, 7, 3, 4, 5, 6), byrow = TRUE, ncol = 2),
           heights = c(1/20, 19/60, 19/60, 19/60), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  } else {
    png(paste0("ThreshVsFPR_", pop_name, "_linear_environment.png"),
        height = 3000, width = 3000, res = 300, pointsize = 10)
    layout(matrix(c(1, 1, 2, 5, 3, 4), byrow = TRUE, ncol = 2),
           heights = c(1/10, 9/20, 9/20), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  }
  
  colours   <- c("Driftsel" = "#6E0D25", "QSTFST" = "#F2C14E", "MONET" = "#7FB7E6", "MONET w/ environment" = "#053225", "H-test" = "#D1495B")
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
            d <- tmp[tmp$method == method, ]
            if (!nrow(d)) next
            
            if (first) {
              par(mar = c(5.1, 4.1, 4.1, 2.1))
              plot(d$threshold, d$FPR, type = "l",
                   lty = linetypes[as.character(wvar_val)],
                   lwd = line_width,
                   col = colours[method],
                   ylim = c(0, 1),
                   xlim = range(plot_data$threshold, na.rm = TRUE),
                   main = paste(Delta_theta, "=", disp_wdiff(wdiff_val),
                                ", correlation =", gsub("_", "", corr_val)),
                   xlab = "Threshold", ylab = "FPR",
                   cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
              abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
              first <- FALSE
            } else {
              lines(d$threshold, d$FPR,
                    lty = linetypes[as.character(wvar_val)],
                    lwd = line_width,
                    col = colours[method])
            }
          }
        }
      }
      
    } else {
      
      plot_data <- all_roc_data[all_roc_data$wdiff == wdiff_val, ]
      if (!nrow(plot_data)) next
      first <- TRUE
      
      for (wvar_val in unique(plot_data$wvar)) {
        tmp <- plot_data[plot_data$wvar == wvar_val, ]
        for (method in unique(tmp$method)) {
          d <- tmp[tmp$method == method, ]
          if (!nrow(d)) next
          
          if (first) {
            par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(d$threshold, d$FPR, type = "l",
                 lty = linetypes[as.character(wvar_val)],
                 lwd = line_width,
                 col = colours[method],
                 ylim = c(0, 1),
                 xlim = range(plot_data$threshold, na.rm = TRUE),
                 main = paste(Delta_theta, "=", disp_wdiff(wdiff_val)),
                 xlab = "Threshold", ylab = "FPR",
                 cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
            abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
            first <- FALSE
          } else {
            lines(d$threshold, d$FPR,
                  lty = linetypes[as.character(wvar_val)],
                  lwd = line_width,
                  col = colours[method])
          }
        }
      }
    }
  }
  
 
  plot(1, 1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = c("Driftsel", "QSTFST", "MONET", "MONET w/ environment", "H-test", paste0(omega_chr, " = ", c(10, 22, 50))),
         ncol   = 2,
         lty    = c(rep(0, 5), linetypes),
         lwd    = c(rep(NA, 5), rep(line_width, 3)),
         pch    = c(rep(15, 5), rep(NA, 3)),
         col    = c(colours, rep("black", 3)),
         cex = 2, bty = "n")
  
  dev.off()
  
  # Generate PDF version with same layout
  if (six.panels) {
    pdf(paste0("ThreshVsFPR_", pop_name, "_linear_environment.pdf"),
        height = 15, width = 10)
    layout(matrix(c(1, 1, 2, 7, 3, 4, 5, 6), byrow = TRUE, ncol = 2),
           heights = c(1/20, 19/60, 19/60, 19/60), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  } else {
    pdf(paste0("ThreshVsFPR_", pop_name, "_linear_environment.pdf"),
        height = 10, width = 10)
    layout(matrix(c(1, 1, 2, 5, 3, 4), byrow = TRUE, ncol = 2),
           heights = c(1/10, 9/20, 9/20), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  }
  
  colours   <- c("Driftsel" = "#6E0D25", "QSTFST" = "#F2C14E", "MONET" = "#7FB7E6", "MONET w/ environment" = "#053225", "H-test" = "#D1495B")
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
            d <- tmp[tmp$method == method, ]
            if (!nrow(d)) next
            
            if (first) {
              par(mar = c(5.1, 4.1, 4.1, 2.1))
              plot(d$threshold, d$FPR, type = "l",
                   lty = linetypes[as.character(wvar_val)],
                   lwd = line_width,
                   col = colours[method],
                   ylim = c(0, 1),
                   xlim = range(plot_data$threshold, na.rm = TRUE),
                   main = paste(Delta_theta, "=", disp_wdiff(wdiff_val),
                                ", correlation =", gsub("_", "", corr_val)),
                   xlab = "Threshold", ylab = "FPR",
                   cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
              abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
              first <- FALSE
            } else {
              lines(d$threshold, d$FPR,
                    lty = linetypes[as.character(wvar_val)],
                    lwd = line_width,
                    col = colours[method])
            }
          }
        }
      }
      
    } else {
      
      plot_data <- all_roc_data[all_roc_data$wdiff == wdiff_val, ]
      if (!nrow(plot_data)) next
      first <- TRUE
      
      for (wvar_val in unique(plot_data$wvar)) {
        tmp <- plot_data[plot_data$wvar == wvar_val, ]
        for (method in unique(tmp$method)) {
          d <- tmp[tmp$method == method, ]
          if (!nrow(d)) next
          
          if (first) {
            par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(d$threshold, d$FPR, type = "l",
                 lty = linetypes[as.character(wvar_val)],
                 lwd = line_width,
                 col = colours[method],
                 ylim = c(0, 1),
                 xlim = range(plot_data$threshold, na.rm = TRUE),
                 main = paste(Delta_theta, "=", disp_wdiff(wdiff_val)),
                 xlab = "Threshold", ylab = "FPR",
                 cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
            abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
            first <- FALSE
          } else {
            lines(d$threshold, d$FPR,
                  lty = linetypes[as.character(wvar_val)],
                  lwd = line_width,
                  col = colours[method])
          }
        }
      }
    }
  }
  
  plot(1, 1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = c("Driftsel", "QSTFST", "MONET", "MONET w/ environment", "H-test", paste0(omega_chr, " = ", c(10, 22, 50))),
         ncol   = 2,
         lty    = c(rep(0, 5), linetypes),
         lwd    = c(rep(NA, 5), rep(line_width, 3)),
         pch    = c(rep(15, 5), rep(NA, 3)),
         col    = c(colours, rep("black", 3)),
         cex = 2, bty = "n")
  
  dev.off()
  
  list(thresh_fpr_data = all_roc_data, auc_data = all_auc_data)
}

## ------------------------------------
# New function with scenario-specific SD for expected calibration
process_population_structure_with_sd <- function(pop_data, pop_name) {
  
  roc_data_list <- list()
  auc_results   <- list()
  sd_info       <- list()
  
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
  
  for (i in 1:nrow(scenarios)) {
    wdiff_val <- scenarios$wdiff[i]
    wvar_val  <- scenarios$wvar[i]
    corr_val  <- scenarios$correlation[i]
    
    if (is.na(corr_val)) {
      scenario_name <- paste0("wdiff_", wdiff_val, "_wvar_", wvar_val)
      scenario_data <- pop_data[which(pop_data$wdiff == wdiff_val & pop_data$wvar == wvar_val), ]
    } else {
      scenario_name <- paste0("wdiff_", wdiff_val, "_wvar_", wvar_val, "_corr_", gsub("_", "", corr_val))
      scenario_data <- pop_data[which(pop_data$wdiff == wdiff_val &
                                  pop_data$wvar  == wvar_val &
                                  pop_data$correlation == corr_val), ]
    }
    
    cat("Processing (SD version):", scenario_name, "\n")
    if (nrow(scenario_data) == 0) next
    
    # Calculate scenario-specific SD
    scenario_neutral_s <- scenario_data$S_value_Driftsel_Neutral
    scenario_sd <- sd(scenario_neutral_s - 0.5, na.rm = TRUE)
    cat("  Scenario SD:", scenario_sd, "\n")
    
    # Store SD info
    sd_info[[scenario_name]] <- data.frame(
      scenario = scenario_name,
      wdiff = wdiff_val,
      wvar = wvar_val,
      correlation = corr_val,
      sd = scenario_sd
    )
    
    #QSTFST
    roc_qstfst <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_QSTFST", "p_value_QSTFST_Neutral", thresh)
    }))
    roc_qstfst$method <- "QSTFST"
    roc_qstfst$scenario <- scenario_name
    roc_qstfst$wdiff <- wdiff_val
    roc_qstfst$wvar <- wvar_val
    roc_qstfst$correlation <- corr_val
    roc_qstfst <- prepare_for_linear_plot(roc_qstfst)
    
    #MONET
    roc_MONET <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "p_value_MONET", "p_value_MONET_Neutral", thresh)
    }))
    roc_MONET$method <- "MONET"
    roc_MONET$scenario <- scenario_name
    roc_MONET$wdiff <- wdiff_val
    roc_MONET$wvar <- wvar_val
    roc_MONET$correlation <- corr_val
    roc_MONET <- prepare_for_linear_plot(roc_MONET)
    
    #Driftsel (folded S thresholds using scenario-specific SD; y-axis = 2 * folded S)
    scenario_s_fold_thresholds <- S_fold_from_p(p_thresholds, s_sd = scenario_sd)
    roc_driftsel <- do.call(rbind, lapply(scenario_s_fold_thresholds, function(th_sfold) {
      calculate_rates_svalue_folded(
        scenario_data,
        "S_value_Driftsel", "S_value_Driftsel_Neutral",
        th_sfold
      )
    }))
    roc_driftsel$threshold <- 2 * roc_driftsel$threshold
    roc_driftsel$method <- "Driftsel"
    roc_driftsel$scenario <- scenario_name
    roc_driftsel$wdiff <- wdiff_val
    roc_driftsel$wvar <- wvar_val
    roc_driftsel$correlation <- corr_val
    roc_driftsel <- prepare_for_linear_plot(roc_driftsel)
    
    #MONET w/ habitat
    roc_habitat <- do.call(rbind, lapply(p_thresholds, function(thresh) {
      calculate_rates_pvalue(scenario_data, "habitat_p_value", "habitat_p_value_Neutral", thresh)
    }))
    roc_habitat$method <- "MONET w/ environment"
    roc_habitat$scenario <- scenario_name
    roc_habitat$wdiff <- wdiff_val
    roc_habitat$wvar <- wvar_val
    roc_habitat$correlation <- corr_val
    roc_habitat <- prepare_for_linear_plot(roc_habitat)
    
    #AUC
    auc_qstfst   <- calculate_auc(roc_qstfst$FPR,   roc_qstfst$TPR)
    auc_MONET     <- calculate_auc(roc_MONET$FPR,     roc_MONET$TPR)
    auc_driftsel <- calculate_auc(roc_driftsel$FPR, roc_driftsel$TPR)
    auc_habitat  <- calculate_auc(roc_habitat$FPR,  roc_habitat$TPR)
    
    auc_results[[scenario_name]] <- data.frame(
      population_structure = pop_name,
      scenario = scenario_name,
      wdiff = wdiff_val,
      wvar = wvar_val,
      correlation = corr_val,
      QSTFST_AUC = auc_qstfst,
      MONET_AUC   = auc_MONET,
      Driftsel_AUC = auc_driftsel,
      Habitat_AUC = auc_habitat
    )
    
    roc_data_list[[scenario_name]] <- rbind(roc_qstfst, roc_MONET, roc_driftsel, roc_habitat)
  }
  
  all_roc_data <- do.call(rbind, roc_data_list)
  all_roc_data$wdiff <- factor(all_roc_data$wdiff, levels = c("0", "4p6", "10"))
  all_auc_data <- do.call(rbind, auc_results)
  all_sd_data <- do.call(rbind, sd_info)
  
  ## ------------------ Threshold vs FPR with scenario-specific calibration ------------------------------------
  six.panels <- pop_name %in% c("SS", "Hierarchical")
  
  # Generate PNG version
  if (six.panels) {
    png(paste0("ThreshVsFPR_", pop_name, "_linear_environment_sd.png"),
        height = 4500, width = 3000, res = 300, pointsize = 10)
    layout(matrix(c(1, 1, 2, 7, 3, 4, 5, 6), byrow = TRUE, ncol = 2),
           heights = c(1/20, 19/60, 19/60, 19/60), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  } else {
    png(paste0("ThreshVsFPR_", pop_name, "_linear_environment_sd.png"),
        height = 3000, width = 3000, res = 300, pointsize = 10)
    layout(matrix(c(1, 1, 2, 5, 3, 4), byrow = TRUE, ncol = 2),
           heights = c(1/10, 9/20, 9/20), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  }
  
  colours   <- c("Driftsel" = "#6E0D25", "QSTFST" = "#F49D37", "MONET" = "#62929E", "MONET w/ environment" = "#053225")
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
            d <- tmp[tmp$method == method, ]
            if (!nrow(d)) next
            
            if (first) {
              par(mar = c(5.1, 4.1, 4.1, 2.1))
              plot(d$threshold, d$FPR, type = "l",
                   lty = linetypes[as.character(wvar_val)],
                   lwd = line_width,
                   col = colours[method],
                   ylim = c(0, 1),
                   xlim = range(plot_data$threshold, na.rm = TRUE),
                   main = paste(Delta_theta, "=", disp_wdiff(wdiff_val),
                                ", correlation =", gsub("_", "", corr_val)),
                   xlab = "Threshold", ylab = "FPR",
                   cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
              abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
              first <- FALSE
            } else {
              lines(d$threshold, d$FPR,
                    lty = linetypes[as.character(wvar_val)],
                    lwd = line_width,
                    col = colours[method])
            }
          }
        }
      }
      
    } else {
      
      plot_data <- all_roc_data[all_roc_data$wdiff == wdiff_val, ]
      if (!nrow(plot_data)) next
      first <- TRUE
      
      for (wvar_val in unique(plot_data$wvar)) {
        tmp <- plot_data[plot_data$wvar == wvar_val, ]
        
        for (method in unique(tmp$method)) {
          d <- tmp[tmp$method == method, ]
          if (!nrow(d)) next
          
          if (first) {
            par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(d$threshold, d$FPR, type = "l",
                 lty = linetypes[as.character(wvar_val)],
                 lwd = line_width,
                 col = colours[method],
                 ylim = c(0, 1),
                 xlim = range(plot_data$threshold, na.rm = TRUE),
                 main = paste(Delta_theta, "=", disp_wdiff(wdiff_val)),
                 xlab = "Threshold", ylab = "FPR",
                 cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
            abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
            first <- FALSE
          } else {
            lines(d$threshold, d$FPR,
                  lty = linetypes[as.character(wvar_val)],
                  lwd = line_width,
                  col = colours[method])
          }
        }
      }
    }
  }
  
  plot(1, 1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = c("Driftsel", "QSTFST", "MONET", "MONET w/ environment", paste0(omega_chr, " = ", c(10, 22, 50))),
         ncol   = 2,
         lty    = c(rep(0, 4), linetypes),
         lwd    = c(rep(NA, 4), rep(line_width, 3)),
         pch    = c(rep(15, 4), rep(NA, 3)),
         col    = c(colours, rep("black", 3)),
         cex = 2, bty = "n")
  
  dev.off()
  
  # Generate PDF version with same layout
  if (six.panels) {
    pdf(paste0("ThreshVsFPR_", pop_name, "_linear_environment_sd.pdf"),
        height = 15, width = 10)
    layout(matrix(c(1, 1, 2, 7, 3, 4, 5, 6), byrow = TRUE, ncol = 2),
           heights = c(1/20, 19/60, 19/60, 19/60), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  } else {
    pdf(paste0("ThreshVsFPR_", pop_name, "_linear_environment_sd.pdf"),
        height = 10, width = 10)
    layout(matrix(c(1, 1, 2, 5, 3, 4), byrow = TRUE, ncol = 2),
           heights = c(1/10, 9/20, 9/20), widths = c(1/2, 1/2))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(1, 1, type = "n", bty = "n", axes = FALSE)
    text(1, 1, pretty_pop(pop_name), cex = 3)
  }
  
  colours   <- c("Driftsel" = "#6E0D25", "QSTFST" = "#F49D37", "MONET" = "#62929E", "MONET w/ environment" = "#053225")
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
            d <- tmp[tmp$method == method, ]
            if (!nrow(d)) next
            
            if (first) {
              par(mar = c(5.1, 4.1, 4.1, 2.1))
              plot(d$threshold, d$FPR, type = "l",
                   lty = linetypes[as.character(wvar_val)],
                   lwd = line_width,
                   col = colours[method],
                   ylim = c(0, 1),
                   xlim = range(plot_data$threshold, na.rm = TRUE),
                   main = paste(Delta_theta, "=", disp_wdiff(wdiff_val),
                                ", correlation =", gsub("_", "", corr_val)),
                   xlab = "Threshold", ylab = "FPR",
                   cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
              abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
              first <- FALSE
            } else {
              lines(d$threshold, d$FPR,
                    lty = linetypes[as.character(wvar_val)],
                    lwd = line_width,
                    col = colours[method])
            }
          }
        }
      }
      
    } else {
      
      plot_data <- all_roc_data[all_roc_data$wdiff == wdiff_val, ]
      if (!nrow(plot_data)) next
      first <- TRUE
      
      for (wvar_val in unique(plot_data$wvar)) {
        tmp <- plot_data[plot_data$wvar == wvar_val, ]
        
        for (method in unique(tmp$method)) {
          d <- tmp[tmp$method == method, ]
          if (!nrow(d)) next
          
          if (first) {
            par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(d$threshold, d$FPR, type = "l",
                 lty = linetypes[as.character(wvar_val)],
                 lwd = line_width,
                 col = colours[method],
                 ylim = c(0, 1),
                 xlim = range(plot_data$threshold, na.rm = TRUE),
                 main = paste(Delta_theta, "=", disp_wdiff(wdiff_val)),
                 xlab = "Threshold", ylab = "FPR",
                 cex.lab = label_cex, cex.axis = axis_cex, cex.main = main_cex)
            abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
            first <- FALSE
          } else {
            lines(d$threshold, d$FPR,
                  lty = linetypes[as.character(wvar_val)],
                  lwd = line_width,
                  col = colours[method])
          }
        }
      }
    }
  }
  
  plot(1, 1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = c("Driftsel", "QSTFST", "MONET", "MONET w/ environment", paste0(omega_chr, " = ", c(10, 22, 50))),
         ncol   = 2,
         lty    = c(rep(0, 4), linetypes),
         lwd    = c(rep(NA, 4), rep(line_width, 3)),
         pch    = c(rep(15, 4), rep(NA, 3)),
         col    = c(colours, rep("black", 3)),
         cex = 2, bty = "n")
  
  dev.off()
  
  list(thresh_fpr_data = all_roc_data, auc_data = all_auc_data, sd_data = all_sd_data)
}

## ------------------ MY SANITY CHECK: realized FPR at nominal thresholds ------------------

diagnostic_fpr <- function(pop_data, pop_name) {
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
  out <- list()
  for (i in seq_len(nrow(scenarios))) {
    wd <- scenarios$wdiff[i]; wv <- scenarios$wvar[i]; cr <- scenarios$correlation[i]
    if (is.na(cr)) {
      sd_ <- pop_data[which(pop_data$wdiff == wd & pop_data$wvar == wv), ]
    } else {
      sd_ <- pop_data[which(pop_data$wdiff == wd & pop_data$wvar == wv & pop_data$correlation == cr), ]
    }
    if (!nrow(sd_)) next
    fpr_p <- function(col, th) mean(sd_[[col]] < th, na.rm = TRUE)
    s_neu_fold <- pmin(sd_$S_value_Driftsel_Neutral, 1 - sd_$S_value_Driftsel_Neutral)
    h_neu_fold <- pmin(sd_$H_value_Driftsel_Neutral, 1 - sd_$H_value_Driftsel_Neutral)
    out[[length(out) + 1]] <- data.frame(
      pop                = pop_name,
      wdiff              = wd, wvar = wv, correlation = cr,
      QSTFST_FPR_p05     = fpr_p("p_value_QSTFST_Neutral",  0.05),
      MONET_FPR_p05      = fpr_p("p_value_MONET_Neutral",    0.05),
      Habitat_FPR_p05    = fpr_p("habitat_p_value_Neutral",  0.05),
      Driftsel_FPR_sf025 = mean(s_neu_fold < 0.025, na.rm = TRUE),
      Htest_FPR_hf025    = mean(h_neu_fold < 0.025, na.rm = TRUE),
      QSTFST_FPR_p01     = fpr_p("p_value_QSTFST_Neutral",  0.01),
      MONET_FPR_p01      = fpr_p("p_value_MONET_Neutral",    0.01),
      Habitat_FPR_p01    = fpr_p("habitat_p_value_Neutral",  0.01),
      Driftsel_FPR_sf005 = mean(s_neu_fold < 0.005, na.rm = TRUE),
      Htest_FPR_hf005    = mean(h_neu_fold < 0.005, na.rm = TRUE),
      n                  = nrow(sd_)
    )
  }
  do.call(rbind, out)
}

cat("\n\n MY SANITY CHECK: realized FPR at nominal thresholds \n")
diag_all <- rbind(
  diagnostic_fpr(im18,         "IM_18"),
  diagnostic_fpr(im9,          "IM_9"),
  diagnostic_fpr(ss,           "SS"),
  diagnostic_fpr(hierarchical, "Hierarchical")
)
print(diag_all, row.names = FALSE, digits = 3)
write.csv(diag_all, "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/results_tables/FPR_at_nominal_thresholds_review2026.csv", row.names = FALSE)

## ------------------------------------------------------
# Original figures without calibration line
res_im18 <- process_population_structure(im18, "IM_18")
res_im9  <- process_population_structure(im9,  "IM_9")
res_ss   <- process_population_structure(ss,   "SS")
res_hier <- process_population_structure(hierarchical, "Hierarchical")

# SD CALIBRATION FIGURES (unsure this makes sense tbh)
##cat("\n\n GENERATING SD-SPECIFIC FIGURES \n\n")
#res_im18_sd <- process_population_structure_with_sd(im18, "IM_18")
#res_im9_sd  <- process_population_structure_with_sd(im9,  "IM_9")
#res_ss_sd   <- process_population_structure_with_sd(ss,   "SS")
#res_hier_sd <- process_population_structure_with_sd(hierarchical, "Hierarchical")


all_auc_results <- rbind(
  res_im18$auc_data,
  res_im9$auc_data,
  res_ss$auc_data,
  res_hier$auc_data
)

print(all_auc_results, row.names = FALSE)
#write.csv(all_auc_results, "AUC_summary_table.csv", row.names = FALSE)
## ------------------ COMBINED FIGURE ------------------
comb_colours   <- c("Driftsel" = "#6E0D25", "H-test" = "#D1495B", "MONET" = "#7FB7E6",
                    "MONET w/ environment" = "#053225", "QSTFST" = "#F2C14E")
comb_labels    <- c("Driftsel (S-test)", "Driftsel (H-test)", "MONET", "MONET w/ habitat", "Qst-Fst")
comb_linetypes <- c("10" = 1, "22" = 2, "50" = 3)
corr_label     <- c("_cline" = "cline", "_sine" = "parabola", "_grouped" = "grouped", "_swapped" = "swapped")

comb_heights <- c(0.25, 1, 1, 1, 1, 1)
comb_widths  <- c(0.3, 1, 1, 1)
comb_layout  <- matrix(c(1,  2,  3,  4,
                         5,  6,  7,  8,
                         9, 10, 11, 12,
                         9, 13, 14, 15,
                        16, 17, 18, 19,
                        16, 20, 21, 22), byrow = TRUE, ncol = 4)

comb_panel <- function(res, wd, corr = NULL) {
  d <- res$thresh_fpr_data
  d <- d[d$wdiff == wd, ]
  if (!is.null(corr)) d <- d[!is.na(d$correlation) & d$correlation == corr, ]
  par(mar = c(4.5, 4.5, 3, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Threshold", ylab = "FPR",
       main = if (is.null(corr)) "" else paste("Correlation =", corr_label[corr]),
       font.main = 1, cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.6)
  abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = line_width)
  for (wv in unique(d$wvar)) for (m in names(comb_colours)) {
    dm <- d[d$wvar == wv & d$method == m, ]
    lines(dm$threshold, dm$FPR, col = comb_colours[m], lty = comb_linetypes[as.character(wv)], lwd = line_width)
  }
}

comb_text <- function(lab, srt = 0, cex = 3.2, font = 2) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, lab, srt = srt, cex = cex, font = font)
}

draw_combined <- function() {
  layout(comb_layout, heights = comb_heights, widths = comb_widths)
  par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
  plot.new()
  ss_top <- 1 - sum(comb_heights[1:2]) / sum(comb_heights)
  ss_bot <- 1 - sum(comb_heights[1:4]) / sum(comb_heights)
  rect(grconvertX(0, "ndc"), grconvertY(ss_bot, "ndc"), grconvertX(1, "ndc"), grconvertY(ss_top, "ndc"),
       col = "grey96", border = NA, xpd = NA)
  for (wd in c("0", "4p6", "10")) comb_text(paste(Delta_theta, "=", disp_wdiff(wd)))

  comb_text("Island Model", srt = 90, font = 1)
  for (wd in c("0", "4p6", "10")) comb_panel(res_im18, wd)

  comb_text("Stepping stones", srt = 90, font = 1)
  comb_panel(res_ss, "0")
  comb_panel(res_ss, "4p6", "_cline")
  comb_panel(res_ss, "10",  "_cline")
  par(mar = c(4.5, 4.5, 3, 1))
  plot.new()
  legend("center",
         legend = c(comb_labels, as.expression(lapply(c(10, 22, 50), function(v) bquote("\u03C9"[o] == .(v))))),
         lty    = c(rep(0, 5), comb_linetypes),
         lwd    = c(rep(NA, 5), rep(line_width, 3)),
         pch    = c(rep(15, 5), rep(NA, 3)),
         pt.cex = 2.5, seg.len = 3,
         col    = c(comb_colours, rep("black", 3)),
         cex = 1.7, bg = "white")
  comb_panel(res_ss, "4p6", "_sine")
  comb_panel(res_ss, "10",  "_sine")

  comb_text("Hierarchical", srt = 90, font = 1)
  comb_panel(res_hier, "0")
  comb_panel(res_hier, "4p6", "_grouped")
  comb_panel(res_hier, "10",  "_grouped")
  plot.new()
  comb_panel(res_hier, "4p6", "_swapped")
  comb_panel(res_hier, "10",  "_swapped")
}

comb_w <- 15
comb_h <- comb_w * sum(comb_heights) / sum(comb_widths)
png("Figure_Results_FPR_Threshold.png", width = comb_w * 300, height = comb_h * 300, res = 300, pointsize = 12)
draw_combined()
dev.off()
cairo_pdf("Figure_Results_FPR_Threshold.pdf", width = comb_w, height = comb_h, pointsize = 12)
draw_combined()
dev.off()
