# MONET analysis WITHOUT environmental covariates



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
  file.path(basedir, "results_logAV_only")
)

stopifnot(species %in% c("faginea", "lusitanica"))
stopifnot(garden  %in% c("GH", "OC"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(basedir)

cat(" monet logAV only (no covariates)\n")
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

missing <- !vapply(list(f_thetaP, f_theM, f_ped, f_pheno),
                   file.exists, logical(1))
if (any(missing)) {
  stop("Missing required input file(s):\n  ",
       paste(c(f_thetaP, f_theM, f_ped, f_pheno)[missing],
             collapse = "\n  "))
}

# ---- load matrices ---------------------------------------------------
Theta.P <- readRDS(f_thetaP)
The.M   <- readRDS(f_theM)
ped     <- read.csv(f_ped,   stringsAsFactors = FALSE)
pheno   <- read.csv(f_pheno, stringsAsFactors = FALSE)

stopifnot(nrow(Theta.P) == ncol(Theta.P))
stopifnot(nrow(The.M)   == ncol(The.M))
stopifnot(!is.null(rownames(Theta.P)), !is.null(rownames(The.M)))

# ---- prepare pedigree ------------------------------------------------
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

# Reorder pedigree so individuals appear contiguously by population, in
# the same order as Theta.P (monet expects block-structured M).
ped$pop <- factor(ped$sire_pop, levels = pop_in_matrix)
ord     <- order(ped$pop)
ped     <- ped[ord, ]
The.M   <- The.M[ped$id, ped$id]

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
      formula_covariates = NULL,  # No covariates
      iter   = 5000,
      warmup = 2000,
      thin   = 2,
      save_full_model = FALSE
    ),
    error = function(e) {
      cat("  monet failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) next


  out_rds <- file.path(outdir,
                       sprintf("monet_%s_%s_%s.rds", species, garden, tid))
  saveRDS(fit, out_rds)


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
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  cat(sprintf("  log_ratio = %.3f  [%.3f, %.3f]  p = %.3f\n\n",
              fit$log_ratio$mean, fit$log_ratio$ci_lower,
              fit$log_ratio$ci_upper, fit$log_ratio$p_value))
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  out_csv <- file.path(outdir,
                       sprintf("summary_logAV_only_%s_%s.csv",
                               species, garden))
  write.csv(summary_df, out_csv, row.names = FALSE)
  cat("Wrote summary:", out_csv, "\n")
} else {
  cat("error on all traits.\n")
}
