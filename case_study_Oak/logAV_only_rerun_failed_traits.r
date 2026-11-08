# Rerun logAV_only (no covariates) monet analysis for traits with poor convergence.
#
# Identifies failing traits directly from
#   <outdir>/summary_logAV_only_<species>_<garden>.csv
# using these criteria (overridable via env vars):
#   max_rhat   > RHAT_THRESHOLD   (default 1.05)
#   n_divergent > DIV_THRESHOLD   (default 100)
#
# Uses MCMC settings:
#   warmup = 4000, iter = 12000, thin = 4   (~8000 draws)
#   adapt_delta = 0.99, max_treedepth = 15
#
# Configuration via environment variables:
#   MONET_SPECIES, MONET_GARDEN, MONET_BASEDIR,
#   MONET_TOOLSDIR, MONET_OUTDIR,
#   RHAT_THRESHOLD, DIV_THRESHOLD


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
  file.path(basedir, "results_logAV_only")
)

rhat_threshold <- as.numeric(Sys.getenv("RHAT_THRESHOLD", "1.05"))
div_threshold  <- as.numeric(Sys.getenv("DIV_THRESHOLD",  "100"))

stopifnot(species %in% c("faginea", "lusitanica"))
stopifnot(garden  %in% c("GH", "OC"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(basedir)

cat("monet logAV_only rerun\n")
cat("species :", species, "\n")
cat("garden  :", garden,  "\n")
cat("basedir :", basedir, "\n")
cat("toolsdir:", toolsdir, "\n")
cat("outdir  :", outdir, "\n")
cat(sprintf("criteria: max_rhat > %.3f OR n_divergent > %d\n",
            rhat_threshold, div_threshold))
cat("MCMC    : warmup=4000, iter=12000, thin=4\n")
cat("          adapt_delta=0.99, max_treedepth=15\n\n")

#identify failing traits from summary csv
summary_csv <- file.path(
  outdir,
  sprintf("summary_logAV_only_%s_%s.csv", species, garden)
)
if (!file.exists(summary_csv)) {
  stop("Cannot find summary file: ", summary_csv,
       "\nRun logAV_only.r first.")
}

summ <- read.csv(summary_csv, stringsAsFactors = FALSE)
needed_cols <- c("trait", "max_rhat", "n_divergent")
miss_cols <- setdiff(needed_cols, names(summ))
if (length(miss_cols) > 0) {
  stop("summary csv missing columns: ", paste(miss_cols, collapse = ", "))
}

bad <- summ[
  (!is.na(summ$max_rhat)    & summ$max_rhat    > rhat_threshold) |
  (!is.na(summ$n_divergent) & summ$n_divergent > div_threshold),
  , drop = FALSE
]

if (nrow(bad) == 0) {
  cat("All traits converged for", species, garden, "- nothing to rerun.\n")
  quit(save = "no", status = 0)
}

cat("Found", nrow(bad), "trait(s) needing rerun:\n")
for (i in seq_len(nrow(bad))) {
  cat(sprintf("  %-25s  max_rhat=%.3f  n_div=%d\n",
              bad$trait[i], bad$max_rhat[i], bad$n_divergent[i]))
}
cat("\n")

#monet source
source(file.path(toolsdir, "counting_blocks_matrix.r"))
source(file.path(toolsdir, "monet.r"))

#input file paths 

# Both faginea gardens (GH, OC) share the same ThetaP from the GH file.
f_thetaP <- file.path(
  basedir,
  sprintf("ramirez-valiente2022_%s_GH_ThetaP.rds", species)
)
f_theM   <- file.path(basedir, paste0(prefix, "_TheM.rds"))
f_ped    <- file.path(basedir, paste0(prefix, "_pedigree_standardized.csv"))
f_pheno  <- file.path(basedir, paste0(prefix, "_phenotypes_standardized.csv"))

missing <- !vapply(list(f_thetaP, f_theM, f_ped, f_pheno),
                   file.exists, logical(1))
if (any(missing)) {
  stop("Missing required input file(s):\n  ",
       paste(c(f_thetaP, f_theM, f_ped, f_pheno)[missing],
             collapse = "\n  "))
}

#load matrices
Theta.P <- readRDS(f_thetaP)
The.M   <- readRDS(f_theM)
ped     <- read.csv(f_ped,   stringsAsFactors = FALSE)
pheno   <- read.csv(f_pheno, stringsAsFactors = FALSE)

stopifnot(nrow(Theta.P) == ncol(Theta.P))
stopifnot(nrow(The.M)   == ncol(The.M))
stopifnot(!is.null(rownames(Theta.P)), !is.null(rownames(The.M)))

# ---- prepare pedigree ------------------------------------------------
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

ped$pop <- factor(ped$sire_pop, levels = pop_in_matrix)
ord     <- order(ped$pop)
ped     <- ped[ord, ]
The.M   <- The.M[ped$id, ped$id]

#iterate over problematic traits
trait_col_name <- if ("trait_id" %in% names(pheno)) "trait_id" else "trait"
if (!all(c("individual", "trait") %in% names(pheno))) {
  stop("Phenotype file must have columns 'individual' and 'trait'.")
}
pheno$individual <- as.character(pheno$individual)

summary_rows <- list()

for (i in seq_len(nrow(bad))) {
  tid <- bad$trait[i]
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

  M_sub   <- The.M[pheno_t$individual, pheno_t$individual]
  ped_sub <- ped[match(pheno_t$individual, ped$id), ]
  ord_sub <- order(ped_sub$pop)
  pheno_t <- pheno_t[ord_sub, ]
  ped_sub <- ped_sub[ord_sub, ]
  M_sub   <- M_sub[pheno_t$individual, pheno_t$individual]

  trait_df <- data.frame(
    id         = pheno_t$individual,
    trait      = pheno_t$trait,
    population = as.character(ped_sub$pop),
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
      formula_covariates = NULL,  # No covariates (logAV_only)
      # STRICT MCMC settings
      iter   = 12000,
      warmup = 4000,
      thin   = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      save_full_model = TRUE  # keep full model so ESS can be extracted
    ),
    error = function(e) {
      cat("  monet failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) next

  # Overwrite previous (bad) rds for this trait
  out_rds <- file.path(outdir,
                       sprintf("monet_%s_%s_%s.rds", species, garden, tid))
  saveRDS(fit, out_rds)

  # Extract min ESS from full model (if present)
  min_bulk_ess <- NA_real_
  min_tail_ess <- NA_real_
  if (!is.null(fit$full_model)) {
    brms_summary <- summary(fit$full_model)
    all_ess_bulk <- c()
    all_ess_tail <- c()
    if (!is.null(brms_summary$fixed)) {
      if ("Bulk_ESS" %in% colnames(brms_summary$fixed))
        all_ess_bulk <- c(all_ess_bulk, brms_summary$fixed[, "Bulk_ESS"])
      if ("Tail_ESS" %in% colnames(brms_summary$fixed))
        all_ess_tail <- c(all_ess_tail, brms_summary$fixed[, "Tail_ESS"])
    }
    if (!is.null(brms_summary$random)) {
      for (re_name in names(brms_summary$random)) {
        re <- brms_summary$random[[re_name]]
        if ("Bulk_ESS" %in% colnames(re))
          all_ess_bulk <- c(all_ess_bulk, re[, "Bulk_ESS"])
        if ("Tail_ESS" %in% colnames(re))
          all_ess_tail <- c(all_ess_tail, re[, "Tail_ESS"])
      }
    }
    if (length(all_ess_bulk) > 0) min_bulk_ess <- min(all_ess_bulk, na.rm = TRUE)
    if (length(all_ess_tail) > 0) min_tail_ess <- min(all_ess_tail, na.rm = TRUE)
  }

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
  cat(sprintf("  Diagnostics: max_rhat=%.3f, n_div=%d, min_bulk_ess=%.0f, min_tail_ess=%.0f\n\n",
              fit$convergence$max_rhat, fit$convergence$n_divergent,
              min_bulk_ess, min_tail_ess))
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  out_csv <- file.path(outdir,
                       sprintf("summary_logAV_only_RERUN_%s_%s.csv",
                               species, garden))
  write.csv(summary_df, out_csv, row.names = FALSE)
  cat("Wrote rerun summary:", out_csv, "\n")
} else {
  cat("No traits produced a summary row.\n")
}
