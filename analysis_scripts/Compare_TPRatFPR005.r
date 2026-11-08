# Compare the power (TPR) of the three methods (QSTFST, MONET, Driftsel)
# at an empirical FPR of 0.05, i.e. the point FPR = 0.05 of the ROC curves.
#
# For each scenario (population structure x wdiff x wvar x correlation) and
# method, the threshold is calibrated on the 500 neutral traits so that 5% of
# them are (falsely) flagged; power = proportion of selected traits flagged at
# that threshold. Statistic used (smaller = more extreme):
#   QSTFST / MONET : p-value
#   Driftsel       : folded S = min(S, 1 - S)
#
# Each replicate gives a 0/1 outcome (selected trait detected or not) per
# method; power is the mean of these. Per scenario: mean power + 95% CI per
# method, then a one-way ANOVA (TP ~ method) with Tukey HSD for the three
# pairwise comparisons (Tukey already corrects for multiple testing).

## ------------------ SETUP ------------------
df <- read.csv("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/raw_data/3methods_full_breeding_combined_november.tsv", header = TRUE, sep = "\t")

df$wdiff <- as.character(df$wdiff)
df$wdiff[df$wdiff == "0.0"]  <- "0"
df$wdiff[df$wdiff == "10.0"] <- "10"
df <- df[!is.na(df$wdiff) & !is.na(df$wvar), ]

df$pop <- NA
df$pop[df$model == "IM" & df$num_pops == 18]           <- "IM_18"
df$pop[df$model == "SS" & df$num_pops == 18]           <- "SS"
df$pop[df$model == "hierarchical" & df$num_pops == 18] <- "Hierarchical"
df <- df[!is.na(df$pop), ]

# Correlation only matters for SS / Hierarchical with wdiff > 0
df$correlation <- gsub("_", "", trimws(as.character(df$correlation)))
df$correlation[is.na(df$correlation) | df$correlation == "" |
               df$pop == "IM_18" | df$wdiff == "0"] <- "none"

target_fpr <- 0.05
methods    <- c("QSTFST", "MONET", "Driftsel")

# Statistic per method for the selected and neutral traits (smaller = more extreme)
df$sel_QSTFST   <- df$p_value_QSTFST
df$neu_QSTFST   <- df$p_value_QSTFST_Neutral
df$sel_MONET    <- df$p_value_MONET
df$neu_MONET    <- df$p_value_MONET_Neutral
df$sel_Driftsel <- pmin(df$S_value_Driftsel, 1 - df$S_value_Driftsel)
df$neu_Driftsel <- pmin(df$S_value_Driftsel_Neutral, 1 - df$S_value_Driftsel_Neutral)

# Threshold t such that mean(neutral < t) is as close as possible to, but not
# above, target_fpr (ties in the neutral values can make it slightly below)
calibrate_threshold <- function(neu, fpr = target_fpr) {
  neu <- sort(neu[!is.na(neu)])
  k   <- floor(fpr * length(neu))
  t   <- neu[k + 1]
  while (mean(neu < t) > fpr) t <- max(neu[neu < t])
  t
}

## ------------------ MEAN + 95% CI, ANOVA, TUKEY ------------------
scenarios <- unique(df[, c("pop", "wdiff", "wvar", "correlation")])
scenarios$pop   <- factor(scenarios$pop,   levels = c("IM_18", "SS", "Hierarchical"))
scenarios$wdiff <- factor(scenarios$wdiff, levels = c("0", "4p6", "10"))
scenarios <- scenarios[order(scenarios$pop, scenarios$wdiff, scenarios$wvar, scenarios$correlation), ]
scenarios$pop   <- as.character(scenarios$pop)
scenarios$wdiff <- as.character(scenarios$wdiff)

summ_list <- list(); anova_list <- list(); tukey_list <- list()

for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  d  <- df[df$pop == sc$pop & df$wdiff == sc$wdiff &
           df$wvar == sc$wvar & df$correlation == sc$correlation, ]
  lab <- data.frame(pop = sc$pop, wdiff = sc$wdiff, wvar = sc$wvar,
                    correlation = sc$correlation)

  # 0/1 detection of the selected trait at the calibrated threshold
  long <- do.call(rbind, lapply(methods, function(m) {
    t   <- calibrate_threshold(d[[paste0("neu_", m)]])
    sel <- d[[paste0("sel_", m)]]
    data.frame(method = m, threshold = t,
               realized_FPR = mean(d[[paste0("neu_", m)]] < t, na.rm = TRUE),
               TP = as.integer(sel[!is.na(sel)] < t))
  }))
  long$method <- factor(long$method, levels = methods)

  # Mean power and 95% CI (exact binomial CI, stays within [0, 1])
  summ_list[[i]] <- do.call(rbind, lapply(methods, function(m) {
    x  <- long$TP[long$method == m]
    ci <- binom.test(sum(x), length(x))$conf.int
    cbind(lab, method = m, n = length(x),
          threshold = long$threshold[long$method == m][1],
          realized_FPR = long$realized_FPR[long$method == m][1],
          power = mean(x), CI_low = ci[1], CI_high = ci[2])
  }))

  # One-way ANOVA + Tukey HSD (not testable if every replicate has the same
  # outcome for all three methods, e.g. power = 1 everywhere)
  if (length(unique(long$TP)) == 1) {
    anova_list[[i]] <- cbind(lab, F = NA, df1 = NA, df2 = NA, p_value = NA)
    tukey_list[[i]] <- cbind(lab, comparison = c("MONET-QSTFST", "Driftsel-QSTFST", "Driftsel-MONET"),
                             diff = 0, CI_low = NA, CI_high = NA, p_adj = NA)
    next
  }
  fit <- aov(TP ~ method, data = long)
  a   <- summary(fit)[[1]]
  anova_list[[i]] <- cbind(lab, F = a$`F value`[1], df1 = a$Df[1],
                           df2 = a$Df[2], p_value = a$`Pr(>F)`[1])

  tk <- TukeyHSD(fit, "method")$method
  tukey_list[[i]] <- cbind(lab, comparison = rownames(tk),
                           diff = tk[, "diff"], CI_low = tk[, "lwr"],
                           CI_high = tk[, "upr"], p_adj = tk[, "p adj"])
}

summ_tab  <- do.call(rbind, summ_list)
anova_tab <- do.call(rbind, anova_list)
tukey_tab <- do.call(rbind, tukey_list)
rownames(summ_tab) <- rownames(anova_tab) <- rownames(tukey_tab) <- NULL

cat("\n Power (TPR) at empirical FPR = 0.05 with 95% CI \n")
print(summ_tab, digits = 3, row.names = FALSE)
cat("\n One-way ANOVA: TP ~ method \n")
print(anova_tab, digits = 3, row.names = FALSE)
cat("\n Tukey HSD pairwise comparisons ====\n")
print(tukey_tab, digits = 3, row.names = FALSE)

write.csv(summ_tab,  "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results/TPRatFPR005_mean_CI.csv", row.names = FALSE)
write.csv(anova_tab, "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results/TPRatFPR005_ANOVA.csv",   row.names = FALSE)
write.csv(tukey_tab, "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results/TPRatFPR005_Tukey.csv",   row.names = FALSE)
