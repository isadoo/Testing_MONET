# MONET analysis with country (China vs USA) as a fixed effect
#
#   MONET_EXPERIMENT  one of: exp1, exp2
#   MONET_BASEDIR     case_study_Geranium data directory (defaults to script dir)
#   MONET_TOOLSDIR    Testing_MONET/tools directory
#   MONET_OUTDIR      output directory (defaults to BASEDIR/results_logAV_country)

# Outputs: - one .rds per trait with the monet result object
# - summary CSV.


suppressPackageStartupMessages({
  library(brms)
})

# parallel

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = n_cores)
cat("Using", n_cores, "cores for parallel chain sampling\n\n")
Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1",
           OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

#configuration 
experiment <- Sys.getenv("MONET_EXPERIMENT", "exp1")
trait_filter <- Sys.getenv("MONET_TRAIT", "")  # if set, run only this trait
basedir    <- Sys.getenv("MONET_BASEDIR", getwd())
toolsdir   <- Sys.getenv(
  "MONET_TOOLSDIR",
  normalizePath(file.path(basedir, "..", "..", "tools"), mustWork = FALSE)
)
outdir     <- Sys.getenv(
  "MONET_OUTDIR",
  file.path(basedir, "results_logAV_country")
)

stopifnot(experiment %in% c("exp1", "exp2"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(basedir)

cat(" monet logAV + country \n")
cat("experiment:", experiment, "\n")
cat("basedir   :", basedir, "\n")
cat("toolsdir  :", toolsdir, "\n")
cat("outdir    :", outdir, "\n\n")

#monet source 
source(file.path(toolsdir, "counting_blocks_matrix.r"))
source(file.path(toolsdir, "monet.r"))

#input file paths 
prefix <- sprintf("shirk2014_%s", experiment)

f_thetaP <- file.path(basedir, paste0(prefix, "_ThetaP.rds"))
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

#build per-individual country frame 
# Populations coded "CN-*" are Chinese; every other population is USA.
# sire_pop and dam_pop are identical for all individuals in this pedigree.
stopifnot(identical(ped$sire_pop, ped$dam_pop))

# Restrict pedigree to individuals present in The.M, in matrix order
ped$id <- as.character(ped$id)
m_ids  <- rownames(The.M)
ped    <- ped[match(m_ids, ped$id), ]
if (any(is.na(ped$id))) {
  stop("Some The.M individuals are missing from the pedigree.")
}

#check populations
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

#reorder pedigree so individuals appear in the same order as Theta.P (monet expects block-structured M).
ped$pop <- factor(ped$sire_pop, levels = pop_in_matrix)
ord     <- order(ped$pop)
ped     <- ped[ord, ]
The.M   <- The.M[ped$id, ped$id]

# Per-individual country label ("USA" is the reference level)
ped$country <- factor(
  ifelse(grepl("^CN", as.character(ped$pop)), "China", "USA"),
  levels = c("USA", "China")
)

cat("Population -> country mapping:\n")
print(table(ped$pop, ped$country))
cat("\n")

#iterate over traits
trait_col_name <- if ("trait_id" %in% names(pheno)) "trait_id" else "trait"
if (!all(c("individual", "trait") %in% names(pheno))) {
  stop("Phenotype file must have columns 'individual' and 'trait'.")
}
pheno$individual <- as.character(pheno$individual)

#derive duration traits from raw date columns 
germ_lookup <- setNames(
  pheno$trait[pheno[[trait_col_name]] == "germinated"],
  as.character(pheno$individual[pheno[[trait_col_name]] == "germinated"])
)
planted_lookup <- setNames(
  pheno$trait[pheno[[trait_col_name]] == "planted"],
  as.character(pheno$individual[pheno[[trait_col_name]] == "planted"])
)
derive_pairs <- list(
  c("dff",        "dff_days",   "germ_lookup"),
  c("juvlf_date", "juvlf_days", "planted_lookup")
)
for (pair in derive_pairs) {
  raw_name <- pair[1]; new_name <- pair[2]
  ref_lookup <- get(pair[3])
  idx <- pheno[[trait_col_name]] == raw_name
  pheno$trait[idx] <- pheno$trait[idx] -
    ref_lookup[as.character(pheno$individual[idx])]
  pheno[[trait_col_name]][idx] <- new_name
  cat("Derived", new_name, "=", raw_name, "-", pair[3], "(n =", sum(idx), ")\n")
}
pheno <- pheno[!(pheno[[trait_col_name]] %in% c("planted", "germinated")), ]
cat("\n")

#derive reproductive ratio (repwt / totalAG) 

repwt_idx      <- pheno[[trait_col_name]] == "repwt"
totalAG_lookup <- setNames(
  pheno$trait[pheno[[trait_col_name]] == "totalAG"],
  as.character(pheno$individual[pheno[[trait_col_name]] == "totalAG"])
)
ratio_vals <- pheno$trait[repwt_idx] /
  totalAG_lookup[as.character(pheno$individual[repwt_idx])]
keep <- !is.na(ratio_vals)
new_rows <- pheno[repwt_idx, ][keep, ]
new_rows$trait <- ratio_vals[keep]
new_rows[[trait_col_name]] <- "repro_ratio"
pheno <- rbind(pheno, new_rows)
cat("Derived repro_ratio = repwt / totalAG (n =", nrow(new_rows), ")\n\n")

trait_ids <- sort(unique(pheno[[trait_col_name]]))
cat("Found", length(trait_ids), "trait(s):\n  ",
    paste(trait_ids, collapse = ", "), "\n\n")

if (nzchar(trait_filter)) {
  if (!trait_filter %in% trait_ids) {
    stop("MONET_TRAIT='", trait_filter, "' not found among trait(s) in ",
         f_pheno, ". Available: ", paste(trait_ids, collapse = ", "))
  }
  trait_ids <- trait_filter
  cat("Restricting to single trait (MONET_TRAIT):", trait_filter, "\n\n")
}

#trait-specific family / standardization overrides

# totC/totN (exp2) use a Student-t response instead of gaussian: a small
# number of individuals (e.g. 343, 338, 328 in block E) have %C/%N values
# far outside the biologically plausible range, while their co-measured
# delC13/delN15 (same combustion run) look unremarkable... 
# consistent with some typoe of error feeding the %C/%N calculation
# rather than a real biological outlier. This is published third-party
# data, so rather than excluding those individuals we let a heavy-tailed
# likelihood down-weight their influence instead of removing them.
trait_family <- list(
  boltA = poisson(),
  boltB = poisson(),
  boltC = poisson(),
  totC  = student(),
  totN  = student()
)
count_traits <- c("boltA", "boltB", "boltC")  # families with no sigma param
trait_no_standardize <- c("boltA", "boltB", "boltC")

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
    id         = pheno_t$individual,
    trait      = pheno_t$trait,
    population = as.character(ped_sub$pop),
    country    = ped_sub$country,
    stringsAsFactors = FALSE
  )

  # Theta.P restricted to the populations actually observed in this trait
  pops_used <- intersect(rownames(Theta.P), unique(trait_df$population))
  ThetaP_sub <- Theta.P[pops_used, pops_used]
  if (length(pops_used) < 2) {
    cat("  Skipping: trait observed in <2 populations.\n")
    next
  }
  if (length(unique(trait_df$country)) < 2) {
    cat("  Skipping: trait observed in only one country.\n")
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
      formula_covariates = "country",
      iter   = 12000,
      warmup = 4000,
      thin   = 2,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      prior = if (tid %in% count_traits) NULL else
        c(prior(exponential(1), class = sigma), prior(exponential(1), class = sd)),
      standardize_trait = if (tid %in% trait_no_standardize) FALSE else NULL,
      family = if (!is.null(trait_family[[tid]])) trait_family[[tid]] else gaussian(),
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
                       sprintf("monet_%s_%s.rds", experiment, tid))
  saveRDS(fit, out_rds)

  # Build a one-row summary. 
  #The country dummy column name depends on brms's factor coding (e.g. "countryChina"), so find it by pattern

  cov_p <- fit$covariate_p_values
  country_col <- if (!is.null(cov_p)) grep("^country", names(cov_p), value = TRUE) else character(0)
  country_p <- if (length(country_col) == 1) cov_p[[country_col]] else NA_real_

  # Posterior means/CIs of beta from the minimal samples table
  draws <- fit$sampling
  beta_summary <- function(name) {
    if (is.data.frame(draws) && name %in% names(draws)) {
      v <- draws[[name]]
      c(mean   = mean(v),
        median = median(v),
        lo     = unname(quantile(v, 0.025)),
        hi     = unname(quantile(v, 0.975)))
    } else {
      c(mean = NA_real_, median = NA_real_, lo = NA_real_, hi = NA_real_)
    }
  }
  draw_col <- if (length(country_col) == 1) paste0("b_", country_col) else NA_character_
  bc <- if (!is.na(draw_col)) beta_summary(draw_col) else c(mean = NA_real_, median = NA_real_, lo = NA_real_, hi = NA_real_)

  row_df <- data.frame(
    experiment       = experiment,
    trait            = tid,
    n_individuals    = nrow(trait_df),
    n_populations    = length(pops_used),
    log_ratio_mean   = fit$log_ratio$mean,
    log_ratio_median = fit$log_ratio$median,
    log_ratio_lo     = fit$log_ratio$ci_lower,
    log_ratio_hi     = fit$log_ratio$ci_upper,
    log_ratio_p      = fit$log_ratio$p_value,
    country_beta_mean = bc["mean"],
    country_beta_median = bc["median"],
    country_beta_lo   = bc["lo"],
    country_beta_hi   = bc["hi"],
    country_p         = country_p,
    n_divergent      = fit$convergence$n_divergent,
    max_rhat         = fit$convergence$max_rhat,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  summary_rows[[length(summary_rows) + 1]] <- row_df

  
  trait_csv <- file.path(outdir,
                         sprintf("summary_logAV_country_%s_%s.csv", experiment, tid))
  write.csv(row_df, trait_csv, row.names = FALSE)
  cat("  Wrote:", trait_csv, "\n")

  cat(sprintf("  log_ratio = %.3f  [%.3f, %.3f]  p = %.3f\n",
              fit$log_ratio$mean, fit$log_ratio$ci_lower,
              fit$log_ratio$ci_upper, fit$log_ratio$p_value))
  cat(sprintf("  country beta = %.3f  p = %.3f\n\n", bc["mean"], country_p))
}

#(legacy)
# skipped in single-trait mode 
if (nzchar(trait_filter)) {
  cat("Single-trait mode: per-trait CSV written above; skipping combined",
      "per-experiment CSV to avoid races with other array tasks.\n")
} else if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  out_csv <- file.path(outdir,
                       sprintf("summary_logAV_country_%s.csv", experiment))
  write.csv(summary_df, out_csv, row.names = FALSE)
  cat("Wrote combined summary:", out_csv, "\n")
} else {
  cat("No traits produced a summary row.\n")
}