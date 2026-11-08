# Rerun monet analysis for traits with poor convergence
# 
# Uses STRICT MCMC settings:
#   - warmup = 4000, iter = 12000, thin = 4  (yields ~8000 draws)
#   - adapt_delta = 0.99 (vs default 0.8)
#   - max_treedepth = 15 (vs default 10)
#
# Configuration via environment variables (same as original):
#   MONET_SPECIES, MONET_GARDEN, MONET_BASEDIR, 
#   MONET_TOOLSDIR, MONET_OUTDIR

suppressPackageStartupMessages({
  library(brms)
})

#configuration
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

cat("monet logAV rerun\n")
cat("species :", species, "\n")
cat("garden  :", garden,  "\n")
cat("basedir :", basedir, "\n")
cat("toolsdir:", toolsdir, "\n")
cat("outdir  :", outdir, "\n")
cat("\nMCMC settings: warmup=4000, iter=12000, thin=4\n")
cat("               adapt_delta=0.99, max_treedepth=15\n\n")

#load list of problematic traits
bad_traits_file <- file.path(outdir, "traits_needing_rerun.csv")
if (!file.exists(bad_traits_file)) {
  stop("Cannot find ", bad_traits_file, 
       "\nRun identify_convergence_failures.r first.")
}

bad_traits <- read.csv(bad_traits_file, stringsAsFactors = FALSE)
bad_traits <- bad_traits[bad_traits$species == species & 
                         bad_traits$garden == garden, ]

if (nrow(bad_traits) == 0) {
  cat("No problematic traits for", species, garden, "- nothing to rerun.\n")
  quit(save = "no", status = 0)
}

cat("Found", nrow(bad_traits), "trait(s) needing rerun:\n  ",
    paste(bad_traits$trait, collapse = ", "), "\n\n")

#monet source
source(file.path(toolsdir, "counting_blocks_matrix.r"))
source(file.path(toolsdir, "monet.r"))

#input file paths
prefix <- sprintf("ramirez-valiente2022_%s_%s", species, garden)

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

#load matrices
Theta.P <- readRDS(f_thetaP)
The.M   <- readRDS(f_theM)
ped     <- read.csv(f_ped,   stringsAsFactors = FALSE)
pheno   <- read.csv(f_pheno, stringsAsFactors = FALSE)
popenv  <- read.csv(f_popenv, stringsAsFactors = FALSE,
                    check.names = FALSE)

stopifnot(nrow(Theta.P) == ncol(Theta.P))
stopifnot(nrow(The.M)   == ncol(The.M))
stopifnot(!is.null(rownames(Theta.P)), !is.null(rownames(The.M)))

#build per-individual environment frame
need_cols <- c("PopulationCode", "PC_climate1", "PC_climate2")
if (!all(need_cols %in% names(popenv))) {
  stop("pop_id_withhabitat.csv is missing PC_climate columns.")
}
pop_pc <- unique(popenv[, need_cols])
rownames(pop_pc) <- pop_pc$PopulationCode

ped$id <- as.character(ped$id)
m_ids  <- rownames(The.M)
ped    <- ped[match(m_ids, ped$id), ]
if (any(is.na(ped$id))) {
  stop("Some The.M individuals are missing from the pedigree.")
}

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

ped$pop <- factor(ped$sire_pop, levels = pop_in_matrix)
ord     <- order(ped$pop)
ped     <- ped[ord, ]
The.M   <- The.M[ped$id, ped$id]

ped$PC_climate1 <- pop_pc[as.character(ped$pop), "PC_climate1"]
ped$PC_climate2 <- pop_pc[as.character(ped$pop), "PC_climate2"]

#iterate over problematic traits
trait_col_name <- if ("trait_id" %in% names(pheno)) "trait_id" else "trait"
if (!all(c("individual", "trait") %in% names(pheno))) {
  stop("Phenotype file must have columns 'individual' and 'trait'.")
}
pheno$individual <- as.character(pheno$individual)

summary_rows <- list()

for (i in seq_len(nrow(bad_traits))) {
  tid <- bad_traits$trait[i]
  cat("---- trait:", tid, "(rerun with strict settings) ----\n")
  
  pheno_t <- pheno[pheno[[trait_col_name]] == tid, c("individual", "trait")]
  pheno_t <- pheno_t[!is.na(pheno_t$trait), , drop = FALSE]

  pheno_t <- aggregate(trait ~ individual, data = pheno_t, FUN = mean,
                       na.rm = TRUE)

  pheno_t <- pheno_t[pheno_t$individual %in% rownames(The.M), , drop = FALSE]
  if (nrow(pheno_t) < 10) {
    cat("  Skipping: only", nrow(pheno_t), "matched individuals.\n")
    next
  }

  M_sub <- The.M[pheno_t$individual, pheno_t$individual]
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
      # STRICT MCMC settings
      iter   = 12000,
      warmup = 4000,
      thin   = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      save_full_model = TRUE  # Save full model to extract ESS later
    ),
    error = function(e) {
      cat("  monet failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) next

  # Save full result object (overwrites old bad one)
  out_rds <- file.path(outdir,
                       sprintf("monet_%s_%s_%s.rds", species, garden, tid))
  saveRDS(fit, out_rds)

  # Extract min ESS from full model
  min_bulk_ess <- NA_real_
  min_tail_ess <- NA_real_
  if (!is.null(fit$full_model)) {
    brms_summary <- summary(fit$full_model)
    all_ess_bulk <- c()
    all_ess_tail <- c()
    
    # Fixed effects
    if (!is.null(brms_summary$fixed)) {
      if ("Bulk_ESS" %in% colnames(brms_summary$fixed)) {
        all_ess_bulk <- c(all_ess_bulk, brms_summary$fixed[, "Bulk_ESS"])
      }
      if ("Tail_ESS" %in% colnames(brms_summary$fixed)) {
        all_ess_tail <- c(all_ess_tail, brms_summary$fixed[, "Tail_ESS"])
      }
    }
    
    # Random effects
    if (!is.null(brms_summary$random)) {
      for (re_name in names(brms_summary$random)) {
        re <- brms_summary$random[[re_name]]
        if ("Bulk_ESS" %in% colnames(re)) {
          all_ess_bulk <- c(all_ess_bulk, re[, "Bulk_ESS"])
        }
        if ("Tail_ESS" %in% colnames(re)) {
          all_ess_tail <- c(all_ess_tail, re[, "Tail_ESS"])
        }
      }
    }
    
    if (length(all_ess_bulk) > 0) min_bulk_ess <- min(all_ess_bulk, na.rm = TRUE)
    if (length(all_ess_tail) > 0) min_tail_ess <- min(all_ess_tail, na.rm = TRUE)
  }

  # Build summary row
  cov_p <- fit$covariate_p_values
  pc1_p <- if (!is.null(cov_p) && "b_PC_climate1" %in% names(cov_p))
             cov_p[["b_PC_climate1"]] else NA_real_
  pc2_p <- if (!is.null(cov_p) && "b_PC_climate2" %in% names(cov_p))
             cov_p[["b_PC_climate2"]] else NA_real_

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
    min_bulk_ess     = min_bulk_ess,
    min_tail_ess     = min_tail_ess,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  cat(sprintf("  log_ratio = %.3f  [%.3f, %.3f]  p = %.3f\n",
              fit$log_ratio$mean, fit$log_ratio$ci_lower,
              fit$log_ratio$ci_upper, fit$log_ratio$p_value))
  cat(sprintf("  PC1 beta  = %.3f  p = %.3f   PC2 beta = %.3f  p = %.3f\n",
              b1["mean"], pc1_p, b2["mean"], pc2_p))
  cat(sprintf("  Diagnostics: max_rhat=%.3f, min_bulk_ess=%.0f, min_tail_ess=%.0f\n\n",
              fit$convergence$max_rhat, min_bulk_ess, min_tail_ess))
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  out_csv <- file.path(outdir,
                       sprintf("summary_logAV_env_RERUN_%s_%s.csv",
                               species, garden))
  write.csv(summary_df, out_csv, row.names = FALSE)
  cat("Wrote rerun summary:", out_csv, "\n")
} else {
  cat("No traits produced a summary row.\n")
}
