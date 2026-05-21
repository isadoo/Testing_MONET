# Numerical comparison of our PCA-climate scores with the ones in
# Ramirez-Valiente et al. (2023) Figure 3. RV's per-population PC scores
# are not published as a table; rv_coords.csv contains values read from
# their Figure 3 (accuracy ~ +/- 0.1 unit) and is for ranking/orientation
# checks, not for fitting.

library(vegan)  # for protest()

ours_path <- "pops_with_climate_soilpH_Im.csv"
rv_path   <- "rv_coords.csv"

ours <- read.csv(ours_path, stringsAsFactors = FALSE)
rv   <- read.csv(rv_path,   stringsAsFactors = FALSE)

# Normalize species labels to short tags ("faginea", "lusitanica")
ours$sp_short <- ifelse(grepl("faginea",    ours$Species, ignore.case = TRUE), "faginea",
                 ifelse(grepl("lusitanica", ours$Species, ignore.case = TRUE), "lusitanica",
                        NA))

for (sp in c("faginea", "lusitanica")) {
  o <- ours[ours$sp_short == sp,
            c("PopulationCode", "PC_climate1", "PC_climate2")]
  r <- rv  [rv$species == sp,
            c("pop", "pc1_rv", "pc2_rv")]
  names(r)[1] <- "PopulationCode"

  shared <- intersect(o$PopulationCode, r$PopulationCode)
  o <- o[match(shared, o$PopulationCode), ]
  r <- r[match(shared, r$PopulationCode), ]

  cat(sprintf("Q. %s : %d shared populations\n", sp, length(shared)))
  cat("Populations:", paste(shared, collapse = ", "), "\n")


  # Spearman per axis. PCA sign is arbitrary so |rho| is what matters.
  s11 <- cor(o$PC_climate1, r$pc1_rv, method = "spearman")
  s22 <- cor(o$PC_climate2, r$pc2_rv, method = "spearman")
  s12 <- cor(o$PC_climate1, r$pc2_rv, method = "spearman")
  s21 <- cor(o$PC_climate2, r$pc1_rv, method = "spearman")
  cat("\nSpearman rank correlations (sign is arbitrary in PCA):\n")
  cat(sprintf("  ours PC1 vs RV PC1 : rho = %+0.3f\n", s11))
  cat(sprintf("  ours PC2 vs RV PC2 : rho = %+0.3f\n", s22))
  cat(sprintf("  ours PC1 vs RV PC2 : rho = %+0.3f   (cross-check; should be small)\n", s12))
  cat(sprintf("  ours PC2 vs RV PC1 : rho = %+0.3f   (cross-check; should be small)\n", s21))

  # Procrustes correlation on the 2-D plane: finds the best rotation +
  # reflection + scaling aligning our (PC1, PC2) to RV's (PC1, PC2).
  X <- as.matrix(o[, c("PC_climate1", "PC_climate2")])
  Y <- as.matrix(r[, c("pc1_rv", "pc2_rv")])
  pro <- protest(X, Y, permutations = 999)
  cat(sprintf("\nProcrustes correlation (PC1+PC2 plane): r = %+0.3f, p = %0.3f\n",
              pro$t0, pro$signif))
  cat("(Procrustes handles arbitrary rotation+reflection; r close to 1 = same configuration.)\n")
}