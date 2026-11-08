# Compare the realized FPR at the nominal 0.05 threshold of the three methods
# (QSTFST, MONET, Driftsel) when wdiff = 0, 
# i.e. the point x = 0.05 of the Threshold-vs-FPR curves drawn by FPR_by_Threshold.r.
# 
#
# A neutral trait is "rejected" (false positive) when:
#   QSTFST / MONET : p_Neutral < 0.05
#   Driftsel       : folded S = min(S, 1 - S) < 0.025  (threshold axis = 2 * folded S)
#
# Each replicate gives a 0/1 outcome per method; the FPR is the mean of these.
# Per population structure x wvar (omega_o): mean FPR + 95% CI per method,
# then a one-way ANOVA (FP ~ method) 
# with Tukey HSD for the three pairwise comparisons 

## ------------------ SETUP ------------------
df <- read.csv("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/raw_data/3methods_full_breeding_combined_november.tsv", header = TRUE, sep = "\t")

df$wdiff <- as.character(df$wdiff)
df$wdiff[df$wdiff == "0.0"] <- "0"
df <- df[!is.na(df$wdiff) & df$wdiff == "0" & !is.na(df$wvar), ]

df$pop <- NA
df$pop[df$model == "IM" & df$num_pops == 18]           <- "IM_18"
df$pop[df$model == "SS" & df$num_pops == 18]           <- "SS"
df$pop[df$model == "hierarchical" & df$num_pops == 18] <- "Hierarchical"
df <- df[!is.na(df$pop), ]

alpha   <- 0.05
methods <- c("QSTFST", "MONET", "Driftsel")

# 0/1 false-positive indicators
s_fold <- pmin(df$S_value_Driftsel_Neutral, 1 - df$S_value_Driftsel_Neutral)
df$FP_QSTFST   <- as.integer(df$p_value_QSTFST_Neutral < alpha)
df$FP_MONET    <- as.integer(df$p_value_MONET_Neutral  < alpha)
df$FP_Driftsel <- as.integer(s_fold < alpha / 2)

# Long format: one row per replicate x method
long <- do.call(rbind, lapply(methods, function(m) {
  data.frame(pop = df$pop, wvar = df$wvar, replicate = df$replicate_number,
             method = m, FP = df[[paste0("FP_", m)]])
}))
long <- long[!is.na(long$FP), ]
long$method <- factor(long$method, levels = methods)

## ------------------ MEAN + 95% CI, ANOVA, TUKEY ------------------
summ_list <- list(); anova_list <- list(); tukey_list <- list()

for (pop_name in c("IM_18", "SS", "Hierarchical")) {
 for (wv in sort(unique(long$wvar[long$pop == pop_name]))) {
  d   <- long[long$pop == pop_name & long$wvar == wv, ]
  key <- paste(pop_name, wv, sep = "_")

  # Mean FPR and 95% CI (exact binomial CI, stays within [0, 1])
  summ_list[[key]] <- do.call(rbind, lapply(methods, function(m) {
    x  <- d$FP[d$method == m]
    ci <- binom.test(sum(x), length(x))$conf.int
    data.frame(pop = pop_name, wvar = wv, method = m, n = length(x),
               FPR = mean(x), CI_low = ci[1], CI_high = ci[2])
  }))

  # One-way ANOVA + Tukey HSD
  fit <- aov(FP ~ method, data = d)
  a   <- summary(fit)[[1]]
  anova_list[[key]] <- data.frame(pop = pop_name, wvar = wv, F = a$`F value`[1],
                                       df1 = a$Df[1], df2 = a$Df[2],
                                       p_value = a$`Pr(>F)`[1])

  tk <- TukeyHSD(fit, "method")$method
  tukey_list[[key]] <- data.frame(pop = pop_name, wvar = wv, comparison = rownames(tk),
                                       diff = tk[, "diff"], CI_low = tk[, "lwr"],
                                       CI_high = tk[, "upr"], p_adj = tk[, "p adj"])
 }
}

summ_tab  <- do.call(rbind, summ_list)
anova_tab <- do.call(rbind, anova_list)
tukey_tab <- do.call(rbind, tukey_list)
rownames(summ_tab) <- rownames(anova_tab) <- rownames(tukey_tab) <- NULL

cat("\n Mean FPR at nominal 0.05 (wdiff = 0) with 95% CI \n")
print(summ_tab, digits = 3, row.names = FALSE)
cat("\n One-way ANOVA: FP ~ method \n")
print(anova_tab, digits = 3, row.names = FALSE)
cat("\n Tukey HSD pairwise comparisons \n")
print(tukey_tab, digits = 3, row.names = FALSE)

write.csv(summ_tab,  "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results/FPRat005_wdiff0_mean_CI.csv", row.names = FALSE)
write.csv(anova_tab, "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results/FPRat005_wdiff0_ANOVA.csv",   row.names = FALSE)
write.csv(tukey_tab, "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results/FPRat005_wdiff0_Tukey.csv",   row.names = FALSE)
