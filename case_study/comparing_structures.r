### This script will build a figure with the five different thetaP.
#ThetaP for Q. faginea, Q. lusitanica
#ThetaP for one replicate of Island Model
#ThetaP for one replicate of Stepping stones
#ThetaP for one replicate of hierarchical structure.

library(hierfstat)
library(Matrix)
library(ggplot2)
library(corrplot)

setwd("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET") 

dat_hierarchical <- "case_study/replicates_compared_structures/hierarchical_neutral_data_g1350.dat"
dat_im <- "case_study/replicates_compared_structures/IM_neutral_data_g1000.dat"
dat_ss <- "case_study/replicates_compared_structures/SS_neutral_data_g3000.dat"

Theta_P_lusitanica <- readRDS("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study/ramirez-valiente2022_lusitanica_GH_ThetaP.rds")
Theta_P_faginea <- readRDS("/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study/ramirez-valiente2022_faginea_GH_ThetaP.rds")  

options(bitmapType = "cairo")

calculate_thetaP <- function(dat_file, n_pops = 18, n_ind_per_pop = 10, verbose = FALSE) {
  

  sim <- hierfstat::read.fstat(fname = dat_file)
  dos <- hierfstat::biall2dos(sim[, -1])  
  

  pop_id <- rep(1:n_pops, each = n_ind_per_pop)
  

  matching_matrix <- hierfstat::matching(dos)
  kinship         <- hierfstat::beta.dosage(matching_matrix, MATCHING = TRUE)
  fst             <- hierfstat::fs.dosage(matching_matrix, pop = pop_id, matching = TRUE)
  
  
  min_fst  <- min(hierfstat::mat2vec(fst$FsM))
  theta_p  <- (fst$FsM - min_fst) / (1 - min_fst)
  
  #check positive definiteness and correct if needed
  eigenvalues <- eigen(theta_p)$values
  if (any(eigenvalues < 0)) {
    if (verbose) cat("Theta.P matrix not positive definite.\n")
    if (verbose) cat("Minimum eigenvalue:", min(eigenvalues), "\n")
    theta_p     <- as.matrix(Matrix::nearPD(theta_p)$mat)
    eigenvalues <- eigen(theta_p)$values
    if (verbose) cat("Theta.P corrected. New min eigenvalue:", min(eigenvalues), "\n")
  }
  
  return(theta_p)
}

# Calculate ThetaP for the three simulated structures
cat("Calculating ThetaP for hierarchical structure...\n")
Theta_P_hierarchical <- calculate_thetaP(dat_hierarchical, n_pops = 18, n_ind_per_pop = 10, verbose = TRUE)

cat("Calculating ThetaP for Island Model...\n")
Theta_P_im <- calculate_thetaP(dat_im, n_pops = 18, n_ind_per_pop = 10, verbose = TRUE)

cat("Calculating ThetaP for Stepping Stones...\n")
Theta_P_ss <- calculate_thetaP(dat_ss, n_pops = 18, n_ind_per_pop = 10, verbose = TRUE)

# Ensure all matrices are proper matrix objects
Theta_P_hierarchical <- as.matrix(Theta_P_hierarchical)
Theta_P_im <- as.matrix(Theta_P_im)
Theta_P_ss <- as.matrix(Theta_P_ss)
Theta_P_lusitanica <- as.matrix(Theta_P_lusitanica)
Theta_P_faginea <- as.matrix(Theta_P_faginea)

# Find common color range across all five matrices
theta_min_all <- min(c(Theta_P_hierarchical, Theta_P_im, Theta_P_ss, 
                       Theta_P_lusitanica, Theta_P_faginea), na.rm = TRUE)
theta_max_all <- max(c(Theta_P_hierarchical, Theta_P_im, Theta_P_ss, 
                       Theta_P_lusitanica, Theta_P_faginea), na.rm = TRUE)

cat(sprintf("Shared ThetaP color range: [%.4f, %.4f]\n", theta_min_all, theta_max_all))

# Add small buffer to color limits
col_lim <- c(theta_min_all, theta_max_all) + c(-1e-8, 1e-8)

# Create color palette
n_cols <- 600
base_pal <- grDevices::colorRampPalette(
  c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"),
  space = "Lab"
)(n_cols)
low_emphasis_idx <- pmax(1, round((seq(0, 1, length.out = n_cols)^0.55) * (n_cols - 1)) + 1)
theta_palette <- base_pal[low_emphasis_idx]

# Create multi-panel figure
output_file <- "case_study/ThetaP_comparison_structures.pdf"
cat("Creating figure:", output_file, "\n")

grDevices::cairo_pdf(output_file, width = 18, height = 12)

# Set up layout: 
# Top row: 3 plots (hierarchical, IM, SS)
# Bottom row: 2 plots (lusitanica, faginea) with space for legend in middle
layout_matrix <- matrix(c(
  1, 1, 1, 2, 2, 2, 3, 3, 3,
  1, 1, 1, 2, 2, 2, 3, 3, 3,
  1, 1, 1, 2, 2, 2, 3, 3, 3,
  4, 4, 4, 0, 0, 0, 5, 5, 5,
  4, 4, 4, 0, 0, 0, 5, 5, 5,
  4, 4, 4, 0, 0, 0, 5, 5, 5
), nrow = 6, byrow = TRUE)

layout(layout_matrix)
par(oma = c(1, 1, 3, 1))

#Plot 1: Island Model
par(mar = c(1, 1, 3, 1))
corrplot::corrplot(
  Theta_P_im,
  method  = "color",
  is.corr = FALSE,
  col     = theta_palette,
  col.lim = col_lim,
  tl.pos  = "lt",
  tl.cex  = 1.0,
  tl.col  = "black",
  cl.pos  = "n",
  mar     = c(0, 0, 2, 0)
)
title(main = "Island Model", cex.main = 1.6)

# Plot 2: Hierarchical
par(mar = c(1, 1, 3, 1))
corrplot::corrplot(
  Theta_P_hierarchical,
  method  = "color",
  is.corr = FALSE,
  col     = theta_palette,
  col.lim = col_lim,
  tl.pos  = "lt",
  tl.cex  = 1.0,
  tl.col  = "black",
  cl.pos  = "n",  # no color legend on individual plots
  mar     = c(0, 0, 2, 0)
)
title(main = "Hierarchical Structure", cex.main = 1.6)

# Plot 3: Stepping Stones
par(mar = c(1, 1, 3, 1))
corrplot::corrplot(
  Theta_P_ss,
  method  = "color",
  is.corr = FALSE,
  col     = theta_palette,
  col.lim = col_lim,
  tl.pos  = "lt",
  tl.cex  = 1.0,
  tl.col  = "black",
  cl.pos  = "n",
  mar     = c(0, 0, 2, 0)
)
title(main = "Stepping Stones", cex.main = 1.6)

# Plot 4: Q. lusitanica
par(mar = c(1, 1, 3, 1))
corrplot::corrplot(
  Theta_P_lusitanica,
  method  = "color",
  is.corr = FALSE,
  col     = theta_palette,
  col.lim = col_lim,
  tl.pos  = "lt",
  tl.cex  = 1.0,
  tl.col  = "black",
  cl.pos  = "n",
  mar     = c(0, 0, 2, 0)
)
title(main = "Q. lusitanica", cex.main = 1.6, font.main = 4)

# Plot 5: Q. faginea with legend
par(mar = c(1, 1, 3, 3))
corrplot::corrplot(
  Theta_P_faginea,
  method  = "color",
  is.corr = FALSE,
  col     = theta_palette,
  col.lim = col_lim,
  tl.pos  = "lt",
  tl.cex  = 1.0,
  tl.col  = "black",
  cl.pos  = "r",  # legend on the right
  cl.cex  = 1.2,
  mar     = c(0, 0, 2, 0)
)
title(main = "Q. faginea", cex.main = 1.6, font.main = 4)

# Overall title
mtext(expression("Comparison of " * Theta[P] * " Matrices Across Structures"), 
      side = 3, outer = TRUE, cex = 1.8, line = 1)

dev.off()
cat("Figure saved to:", output_file, "\n")
