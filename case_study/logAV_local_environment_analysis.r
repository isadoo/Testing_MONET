# MONET analysis with local environment as fixed effects
#
# For one (species, garden) configuration, fits the monet model to each
# trait in the long-format phenotype file, including the first two
# climate principal components (PC_climate1, PC_climate2) from
# pop_id_withhabitat.csv as fixed effects.
#
# Two questions are addressed jointly:
#   1. Pattern of local vs global adaptation:
#        is the log ratio of ancestral variances significantly
#        positive (local) or negative (global)?
#   2. Environmental association:
#        are the regression coefficients (betas) of PC_climate1 and
#        PC_climate2 significantly different from zero?
#
# Configuration is taken from environment variables so the same script
# can be used by a SLURM array job:
#
#   MONET_SPECIES   one of: faginea, lusitanica
#   MONET_GARDEN    one of: GH, OC
#   MONET_BASEDIR   case_study directory (defaults to script dir)
#   MONET_TOOLSDIR  Testing_MONET/tools directory
#   MONET_OUTDIR    output directory (defaults to BASEDIR/results_logAV_env)
#
# Outputs: one .rds per trait with the monet result object, plus a
# combined per-configuration summary CSV.


suppressPackageStartupMessages({
  library(brms)
})

# ---- configuration ---------------------------------------------------
species  <- Sys.getenv("MONET_SPECIES", "faginea")
garden   <- Sys.getenv("MONET_GARDEN",  "GH")
basedir  <- Sys.getenv("MONET_BASEDIR", getwd())
toolsdir <- Sys.getenv(
  "MONET_TOOLSDIR",
  normalizePath(file.path(basedir, "..", "tools"), mustWork = FALSE)
)
outdir   <- Sys.getenv(
  "MONET_OUTDIR",
  file.path(basedir, "results_logAV_env")
)

stopifnot(species %in% c("faginea", "lusitanica"))
stopifnot(garden  %in% c("GH", "OC"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(basedir)

cat(" monet logAV + environment \n")
cat("species :", species, "\n")
cat("garden  :", garden,  "\n")
cat("basedir :", basedir, "\n")
cat("toolsdir:", toolsdir, "\n")
cat("outdir  :", outdir, "\n\n")

# ---- monet source ----------------------------------------------------
source(file.path(toolsdir, "counting_blocks_matrix.r"))
source(file.path(toolsdir, "monet.r"))

# ---- input file paths ------------------------------------------------
prefix <- sprintf("ramirez-valiente2022_%s_%s", species, garden)

# Both faginea common gardens (GH and OC) share the same wild-population
# genetic data, so ThetaP is always taken from the GH file for each species.
f_thetaP <- file.path(
  basedir,
  sprintf("ramirez-valiente2022_%s_GH_ThetaP.rds", species)
)
f_theM   <- file.path(basedir, paste0(prefix, "_TheM.rds"))
f_ped    <- file.path(basedir, paste0(prefix, "_pedigree_standardized.csv"))
f_pheno  <- file.path(basedir, paste0(prefix, "_phenotypes_standardized.csv"))
f_popenv <- file.path(basedir, "pop_id_withhabitat.csv")

missing <- !vapply(list(f_thetaP, f_theM, f_ped, f_pheno, f_popenv),
                   file.exists, logical(1))
if (any(missing)) {
  stop("Missing required input file(s):\n  ",
       paste(c(f_thetaP, f_theM, f_ped, f_pheno, f_popenv)[missing],
             collapse = "\n  "))
}

# ---- load matrices ---------------------------------------------------
Theta.P <- readRDS(f_thetaP)
The.M   <- readRDS(f_theM)
ped     <- read.csv(f_ped,   stringsAsFactors = FALSE)
pheno   <- read.csv(f_pheno, stringsAsFactors = FALSE)
popenv  <- read.csv(f_popenv, stringsAsFactors = FALSE,
                    check.names = FALSE)

stopifnot(nrow(Theta.P) == ncol(Theta.P))
stopifnot(nrow(The.M)   == ncol(The.M))
stopifnot(!is.null(rownames(Theta.P)), !is.null(rownames(The.M)))

# ---- build per-individual environment frame --------------------------
# pop_id_withhabitat.csv has one row per individual; collapse to one row
# per population for the climate PCs.
need_cols <- c("PopulationCode", "PC_climate1", "PC_climate2")
if (!all(need_cols %in% names(popenv))) {
  stop("pop_id_withhabitat.csv is missing PC_climate columns. ",
       "Run envrionmental_data.r first.")
}
pop_pc <- unique(popenv[, need_cols])
rownames(pop_pc) <- pop_pc$PopulationCode

# Restrict pedigree to individuals present in The.M, in matrix order
ped$id <- as.character(ped$id)
m_ids  <- rownames(The.M)
ped    <- ped[match(m_ids, ped$id), ]
if (any(is.na(ped$id))) {
  stop("Some The.M individuals are missing from the pedigree.")
}

# Validate populations
pop_in_matrix <- rownames(Theta.P)
pop_in_ped    <- unique(ped$sire_pop)
miss_pop      <- setdiff(pop_in_ped, pop_in_matrix)
if (length(miss_pop) > 0) {
  warning(sprintf(
    "Dropping %d population(s) from analysis (no ThetaP): %s",
    length(miss_pop), paste(miss_pop, collapse = ", ")
  ))
  ped <- ped[!ped$sire_pop %in% miss_pop, ]
  The.M <- The.M[ped$id, ped$id]
}
miss_env <- setdiff(pop_in_matrix, pop_pc$PopulationCode)
if (length(miss_env) > 0) {
  stop("Theta.P pops missing PC_climate values: ",
       paste(miss_env, collapse = ", "))
}

# Reorder pedigree so individuals appear contiguously by population, in
# the same order as Theta.P (monet expects block-structured M).
ped$pop <- factor(ped$sire_pop, levels = pop_in_matrix)
ord     <- order(ped$pop)
ped     <- ped[ord, ]
The.M   <- The.M[ped$id, ped$id]

# Per-individual PC values
ped$PC_climate1 <- pop_pc[as.character(ped$pop), "PC_climate1"]
ped$PC_climate2 <- pop_pc[as.character(ped$pop), "PC_climate2"]

# ---- iterate over traits --------------------------------------------
# Phenotype is long-format: columns include `individual`, `trait`, `trait_id`
trait_col_name <- if ("trait_id" %in% names(pheno)) "trait_id" else "trait"
if (!all(c("individual", "trait") %in% names(pheno))) {
  stop("Phenotype file must have columns 'individual' and 'trait'.")
}
pheno$individual <- as.character(pheno$individual)

trait_ids <- sort(unique(pheno[[trait_col_name]]))
cat("Found", length(trait_ids), "trait(s):\n  ",
    paste(trait_ids, collapse = ", "), "\n\n")

summary_rows <- list()

for (tid in trait_ids) {
  cat("---- trait:", tid, "----\n")
  pheno_t <- pheno[pheno[[trait_col_name]] == tid, c("individual", "trait")]
  pheno_t <- pheno_t[!is.na(pheno_t$trait), , drop = FALSE]

  # Average replicates per individual (some traits have multiple measurements)
  pheno_t <- aggregate(trait ~ individual, data = pheno_t, FUN = mean,
                       na.rm = TRUE)

  # Restrict to individuals present in The.M
  pheno_t <- pheno_t[pheno_t$individual %in% rownames(The.M), , drop = FALSE]
  if (nrow(pheno_t) < 10) {
    cat("  Skipping: only", nrow(pheno_t), "matched individuals.\n")
    next
  }

  # Build the input dataframe in the order required by The.M
  M_sub <- The.M[pheno_t$individual, pheno_t$individual]
  # Re-sort to match population block order in the subsetted set
  ped_sub <- ped[match(pheno_t$individual, ped$id), ]
  ord_sub <- order(ped_sub$pop)
  pheno_t <- pheno_t[ord_sub, ]
  ped_sub <- ped_sub[ord_sub, ]
  M_sub   <- M_sub[pheno_t$individual, pheno_t$individual]

  trait_df <- data.frame(
    id          = pheno_t$individual,
    trait       = pheno_t$trait,
    population  = as.character(ped_sub$pop),
    PC_climate1 = ped_sub$PC_climate1,
    PC_climate2 = ped_sub$PC_climate2,
    stringsAsFactors = FALSE
  )

  # Theta.P restricted to the populations actually observed in this trait
  pops_used <- intersect(rownames(Theta.P), unique(trait_df$population))
  ThetaP_sub <- Theta.P[pops_used, pops_used]
  if (length(pops_used) < 2) {
    cat("  Skipping: trait observed in <2 populations.\n")
    next
  }

  fit <- tryCatch(
    monet(
      Theta.P            = ThetaP_sub,
      M                  = M_sub,
      trait_dataframe    = trait_df,
      column_individual  = "id",
      column_trait       = "trait",
      column_population  = "population",
      formula_covariates = "PC_climate1 + PC_climate2",
      iter   = 12000,
      warmup = 4000,
      thin   = 2,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      prior = c(prior(exponential(1), class = sigma), prior(exponential(1), class = sd)),
      save_full_model = FALSE
    ),
    error = function(e) {
      cat("  monet failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) next

  # Save full result object
  out_rds <- file.path(outdir,
                       sprintf("monet_%s_%s_%s.rds", species, garden, tid))
  saveRDS(fit, out_rds)

  # Build a one-row summary
  cov_p <- fit$covariate_p_values
  pc1_p <- if (!is.null(cov_p) && "b_PC_climate1" %in% names(cov_p))
             cov_p[["b_PC_climate1"]] else NA_real_
  pc2_p <- if (!is.null(cov_p) && "b_PC_climate2" %in% names(cov_p))
             cov_p[["b_PC_climate2"]] else NA_real_

  # Posterior means/CIs of betas from the minimal samples table
  draws <- fit$sampling
  beta_summary <- function(name) {
    if (is.data.frame(draws) && name %in% names(draws)) {
      v <- draws[[name]]
      c(mean = mean(v),
        lo   = unname(quantile(v, 0.025)),
        hi   = unname(quantile(v, 0.975)))
    } else {
      c(mean = NA_real_, lo = NA_real_, hi = NA_real_)
    }
  }
  b1 <- beta_summary("b_PC_climate1")
  b2 <- beta_summary("b_PC_climate2")

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    species          = species,
    garden           = garden,
    trait            = tid,
    n_individuals    = nrow(trait_df),
    n_populations    = length(pops_used),
    log_ratio_mean   = fit$log_ratio$mean,
    log_ratio_median = fit$log_ratio$median,
    log_ratio_lo     = fit$log_ratio$ci_lower,
    log_ratio_hi     = fit$log_ratio$ci_upper,
    log_ratio_p      = fit$log_ratio$p_value,
    PC1_beta_mean    = b1["mean"],
    PC1_beta_lo      = b1["lo"],
    PC1_beta_hi      = b1["hi"],
    PC1_p            = pc1_p,
    PC2_beta_mean    = b2["mean"],
    PC2_beta_lo      = b2["lo"],
    PC2_beta_hi      = b2["hi"],
    PC2_p            = pc2_p,
    n_divergent      = fit$convergence$n_divergent,
    max_rhat         = fit$convergence$max_rhat,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  cat(sprintf("  log_ratio = %.3f  [%.3f, %.3f]  p = %.3f\n",
              fit$log_ratio$mean, fit$log_ratio$ci_lower,
              fit$log_ratio$ci_upper, fit$log_ratio$p_value))
  cat(sprintf("  PC1 beta  = %.3f  p = %.3f   PC2 beta = %.3f  p = %.3f\n\n",
              b1["mean"], pc1_p, b2["mean"], pc2_p))
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  out_csv <- file.path(outdir,
                       sprintf("summary_logAV_env_%s_%s.csv",
                               species, garden))
  write.csv(summary_df, out_csv, row.names = FALSE)
  cat("Wrote summary:", out_csv, "\n")
} else {
  cat("No traits produced a summary row.\n")
}
