# Table 3: TPR at FPR 0.05 for the sampling and breeding tests (wdiff 4.6, wvar 10).
# Run from analysis_scripts/. Same TPR computation as Table 2: ROC on a p-value grid,
# TPR interpolated at FPR = 0.05.

raw <- "../raw_data"
out <- "../results"

p_thresholds <- sort(unique(c(seq(0.0001, 0.001, by = 0.00005), seq(0.001, 0.01, by = 0.0005),
                              seq(0.01, 0.1, by = 0.001), seq(0.1, 1, by = 0.005))))
s_fold_thresholds <- pmin(pmax(0.5 - qnorm(1 - p_thresholds / 2, sd = 0.10), 0), 0.5)

tpr_at_fpr <- function(sel, neu, thresholds, target = 0.05) {
  if (all(is.na(sel)) || all(is.na(neu))) return(NA)
  tpr <- sapply(thresholds, function(t) mean(sel < t, na.rm = TRUE))
  fpr <- sapply(thresholds, function(t) mean(neu < t, na.rm = TRUE))
  o <- order(fpr, tpr); fpr <- fpr[o]; tpr <- tpr[o]
  exact <- abs(fpr - target) < 1e-10
  if (any(exact)) return(mean(tpr[exact]))
  if (!any(fpr <= target) || !any(fpr >= target)) return(NA)
  f1 <- max(fpr[fpr <= target]); f2 <- min(fpr[fpr >= target])
  t1 <- tpr[fpr == f1][1];        t2 <- tpr[fpr == f2][1]
  t1 + (t2 - t1) * (target - f1) / (f2 - f1)
}

fold <- function(s) pmin(s, 1 - s)
col  <- function(d, x) if (x %in% names(d)) d[[x]] else rep(NA_real_, nrow(d))

three_methods <- function(d) c(
  QSTFST   = tpr_at_fpr(col(d, "p_value_QSTFST"), col(d, "p_value_QSTFST_Neutral"), p_thresholds),
  MONET    = tpr_at_fpr(d$p_value_MONET, d$p_value_MONET_Neutral, p_thresholds),
  Driftsel = tpr_at_fpr(fold(d$S_value_Driftsel), fold(d$S_value_Driftsel_Neutral), s_fold_thresholds)
)

read_runs <- function(f) {
  d <- read.delim(file.path(raw, f), stringsAsFactors = FALSE)
  d$wdiff <- sub("\\.0$", "", as.character(d$wdiff))
  d$wvar  <- sub("\\.0$", "", as.character(d$wvar))
  d$correlation <- sub("^_", "", d$correlation)
  d$correlation[d$correlation %in% c("", "NA")] <- NA
  d
}

pick <- function(d, model, num_pops, correlation) {
  keep <- d$model == model & d$num_pops == num_pops & d$wdiff == "4p6" & d$wvar == "10" &
    (if (is.na(correlation)) is.na(d$correlation) else d$correlation %in% correlation)
  d[which(keep), ]
}

structures <- list(
  c("island model (18)",        "IM",           NA),
  c("stepping stones cline",    "SS",           "cline"),
  c("stepping stones parabola", "SS",           "sine"),
  c("hierarchical grouped",     "hierarchical", "grouped"),
  c("hierarchical swapped",     "hierarchical", "swapped")
)

schemes <- list(
  c("Full design",       "3methods_full_breeding_combined_november.tsv"),
  c("Population number", "breeding_run1_popnumber.tsv"),
  c("F1 number",         "breeding_run2_F1number.tsv"),
  c("3x3",               "breeding_run3_3x3.tsv"),
  c("1x9",               "breeding_run4_1x9.tsv")
)

full <- read_runs(schemes[[1]][2])
rows <- list()
for (s in schemes) {
  d_scheme <- read_runs(s[2])
  for (st in structures) {
    label <- st[1]
    d <- pick(d_scheme, st[2], 18, st[3])
    if (s[1] == "Population number" && st[2] == "IM") {
      label <- "island model (9)"
      d <- pick(full, "IM", 9, NA)
    }
    rows[[length(rows) + 1]] <- data.frame(scheme = s[1], structure = label, n = nrow(d), t(three_methods(d)))
  }
}
tab <- do.call(rbind, rows)
tab$QSTFST[tab$scheme == "1x9"] <- NA

csv <- tab
csv[c("QSTFST", "MONET", "Driftsel")] <- round(csv[c("QSTFST", "MONET", "Driftsel")], 3)
write.csv(csv, file.path(out, "Table3_sampling_breeding_TPR_FPR005.csv"), row.names = FALSE)
print(csv, row.names = FALSE)

# full-design rows are rounded as in Table 2 (to 3 decimals first) so that they match it
fmt    <- function(x, full) ifelse(is.na(x), "---", sprintf("%.2f", if (full) round(x, 3) else x))
scheme <- sub("^(\\d)x(\\d)$", "$\\1 \\\\times \\2$", tab$scheme)
body <- unlist(lapply(seq_len(nrow(tab)), function(i) {
  full <- tab$scheme[i] == "Full design"
  c(if (i > 1 && tab$scheme[i] != tab$scheme[i - 1]) "\\midrule",
    sprintf("%s & %s & %s & %s & %s \\\\", scheme[i], tab$structure[i],
            fmt(tab$QSTFST[i], full), fmt(tab$MONET[i], full), fmt(tab$Driftsel[i], full)))
}))
caption <- paste(
  "True positive rate (TPR) at a false positive rate (FPR) of 0.05 for different sampling and breeding schemes and population structures ($\\Delta\\theta = 4.6$, $\\omega_o = 10$).",
  "The full design (18 subpopulations, 50 F1 individuals per subpopulation, $5 \\times 5$ crossing) is the one used in the rest of the study (same values as in Table 2); each other scheme changes one aspect of it.",
  "For the island model, the population number row is the metapopulation simulated with 9 subpopulations.",
  "Missing values indicate either that the method did not reach an FPR as low as 0.05 (Driftsel on $3 \\times 3$ for stepping stones) or that it does not accommodate the breeding scheme ($Q_{ST}$--$F_{ST}$ on $1 \\times 9$)."
)
tex <- c("\\begin{table}[htbp]", "\\centering",
         paste0("\\caption{", caption, "}"), "\\label{tab:sampling_breeding}",
         "\\begin{tabular}{llrrr}", "\\toprule",
         "Sampling or breeding scheme & Population structure & $Q_{ST}$--$F_{ST}$ & MONET: $\\log_{\\mathcal{AV}}$ & Driftsel $S$-test \\\\",
         "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(out, "Table3_sampling_breeding_TPR_FPR005.tex"))
