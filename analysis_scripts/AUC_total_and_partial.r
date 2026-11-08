##SETUP
library(pROC)

results_dir <- "/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/results"
df <- read.csv("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/raw_data/full_breeding_combined_review2026.tsv", header = TRUE, sep = "\t")
df$wdiff <- as.character(df$wdiff)
df <- df[!is.na(df$wdiff) & !is.na(df$wvar), ]

df$pop <- NA
df$pop[df$model == "IM" & df$num_pops == 18]           <- "IM_18"
df$pop[df$model == "SS" & df$num_pops == 18]           <- "SS"
df$pop[df$model == "hierarchical" & df$num_pops == 18] <- "Hierarchical"
df <- df[!is.na(df$pop), ]
df$correlation <- gsub("_", "", ifelse(is.na(df$correlation), "", df$correlation))

max_fpr <- 0.1
fold    <- function(x) pmin(x, 1 - x)
methods <- list(
  QSTFST   = list(sel = df$p_value_QSTFST,         neu = df$p_value_QSTFST_Neutral),
  MONET    = list(sel = df$p_value_MONET,          neu = df$p_value_MONET_Neutral),
  MONET_env = list(sel = df$habitat_p_value,       neu = df$habitat_p_value_Neutral),
  Driftsel = list(sel = fold(df$S_value_Driftsel), neu = fold(df$S_value_Driftsel_Neutral)),
  Htest    = list(sel = fold(df$H_value_Driftsel), neu = fold(df$H_value_Driftsel_Neutral))
)

auc_pair <- function(sel, neu) {
  sel <- sel[!is.na(sel)]; neu <- neu[!is.na(neu)]
  if (!length(sel) || !length(neu)) return(c(pAUC = NA, AUC = NA))
  r <- roc(controls = neu, cases = sel, direction = ">", quiet = TRUE)
  c(pAUC = as.numeric(auc(r, partial.auc = c(1, 1 - max_fpr), partial.auc.focus = "specificity")) / max_fpr,
    AUC  = as.numeric(auc(r)))
}

## PARTIAL AND TOTAL AUC PER SCENARIO
scenarios <- unique(df[, c("pop", "wdiff", "wvar", "correlation")])
scenarios <- scenarios[order(match(scenarios$pop, c("IM_18", "SS", "Hierarchical")),
                             match(scenarios$wdiff, c("0", "4p6", "10")),
                             scenarios$correlation, scenarios$wvar), ]

res <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  sc  <- scenarios[i, ]
  idx <- which(df$pop == sc$pop & df$wdiff == sc$wdiff & df$wvar == sc$wvar & df$correlation == sc$correlation)
  out <- sapply(methods, function(m) auc_pair(m$sel[idx], m$neu[idx]))
  cbind(sc, n = length(idx),
        setNames(as.data.frame(t(out["pAUC", ])), paste0(names(methods), "_pAUC01")),
        setNames(as.data.frame(t(out["AUC", ])),  paste0(names(methods), "_AUC")))
}))
rownames(res) <- NULL
print(res, digits = 3, row.names = FALSE)

id_cols <- c("pop", "wdiff", "wvar", "correlation", "n")
pcols   <- paste0(names(methods), "_pAUC01")
acols   <- paste0(names(methods), "_AUC")
write.csv(res[, c(id_cols, pcols)], file.path(results_dir, "partialAUC_FPR01_review2026.csv"), row.names = FALSE)
write.csv(res[, c(id_cols, acols)], file.path(results_dir, "totalAUC_FPR_review2026.csv"), row.names = FALSE)

## LATEX TABLES
pop_lab  <- c(IM_18 = "Island Model", SS = "Stepping stones", Hierarchical = "Hierarchical")
corr_lab <- c(cline = "cline", sine = "parabola", grouped = "grouped", swapped = "swapped")

write_tex <- function(cols, caption, label, file) {
  rows <- sapply(seq_len(nrow(res)), function(i) {
    v    <- round(unlist(res[i, cols]), 2)
    cell <- ifelse(is.na(v), "--", sprintf("%.2f", v))
    best <- !is.na(v) & v == max(v, na.rm = TRUE)
    cell[best] <- paste0("\\textbf{", cell[best], "}")
    first <- i == 1 || res$pop[i] != res$pop[i - 1]
    paste0(if (first && i > 1) "\\midrule\n" else "",
           if (first) pop_lab[res$pop[i]] else "", " & ",
           ifelse(res$wdiff[i] == "4p6", "4.6", res$wdiff[i]), " & ",
           ifelse(res$correlation[i] == "", "", corr_lab[res$correlation[i]]), " & ",
           res$wvar[i], " & ", paste(cell, collapse = " & "), " \\\\")
  })
  tex <- c(
    "\\begin{table}[ht]",
    "\\centering",
    paste0("\\caption{", caption,
           " ROC curves are computed from the 500 selected and 500 neutral traits of each scenario. ",
           "The highest value of each row is in bold; -- indicates that the method was not applied.}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lllrccccc}",
    "\\toprule",
    "Population structure & $\\Delta\\theta$ & Correlation & $\\omega_o$ & $Q_{ST}$-$F_{ST}$ & MONET & MONET w/ habitat & Driftsel (S-test) & Driftsel (H-test) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(tex, file.path(results_dir, file))
}

write_tex(pcols,
          paste0("Partial area under the ROC curve for FPR between 0 and ", max_fpr,
                 ", divided by ", max_fpr, " so that it equals the mean power over this FPR range ",
                 "(1 for a perfect test, ", max_fpr / 2, " for a random one)."),
          "tab:pauc_fpr01", "partialAUC_FPR01_review2026.tex")
write_tex(acols,
          "Total area under the ROC curve (FPR between 0 and 1; 1 for a perfect test, 0.5 for a random one).",
          "tab:auc_total", "totalAUC_FPR_review2026.tex")
