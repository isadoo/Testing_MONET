# Marginal form (individual effects integrated out) of the logAV_country_analysis.r model:
#   Y ~ country + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = M))
#
#   MONET_EXPERIMENT  exp1 or exp2
#   MONET_TRAIT       if set, run only this trait
#   MONET_BASEDIR     data directory (defaults to working directory)
#   MONET_OUTDIR      output directory (defaults to BASEDIR/results_logAV_country_marginal)

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
outdir     <- Sys.getenv(
  "MONET_OUTDIR",
  file.path(basedir, "results_logAV_country_marginal")
)

stopifnot(experiment %in% c("exp1", "exp2"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(basedir)

cat(" monet logAV + country (marginal) \n")
cat("experiment:", experiment, "\n")
cat("basedir   :", basedir, "\n")
cat("outdir    :", outdir, "\n\n")


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
if (any(repwt_idx)) {
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
}

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

# poisson (boltA-C) and student (totC, totN) traits: use logAV_country_analysis.r
non_gaussian <- c("boltA", "boltB", "boltC", "totC", "totN")
if (trait_filter %in% non_gaussian) {
  stop("MONET_TRAIT='", trait_filter, "' is not gaussian; the marginal form needs a gaussian likelihood")
}
trait_ids <- setdiff(trait_ids, non_gaussian)

countries <- levels(ped$country)
post_summary <- function(v) {
  c(mean = mean(v), median = unname(median(v)),
    lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)),
    p = 2 * mean(sign(v) != sign(median(v))))
}

summary_rows <- list()

for (tid in trait_ids) {
  cat("---- trait:", tid, "----\n")

  pheno_t <- pheno[pheno[[trait_col_name]] == tid, c("individual", "trait")]
  pheno_t <- pheno_t[!is.na(pheno_t$trait), , drop = FALSE]
  pheno_t <- aggregate(trait ~ individual, data = pheno_t, FUN = mean, na.rm = TRUE)
  pheno_t <- pheno_t[pheno_t$individual %in% rownames(The.M), , drop = FALSE]
  if (nrow(pheno_t) < 10) {
    cat("  Skipping: only", nrow(pheno_t), "matched individuals.\n")
    next
  }

  ped_sub <- ped[match(pheno_t$individual, ped$id), ]
  ord_sub <- order(ped_sub$pop)
  pheno_t <- pheno_t[ord_sub, ]
  ped_sub <- ped_sub[ord_sub, ]
  M_sub   <- The.M[pheno_t$individual, pheno_t$individual]

  pops_used  <- intersect(rownames(Theta.P), unique(as.character(ped_sub$pop)))
  ThetaP_sub <- Theta.P[pops_used, pops_used]
  if (length(pops_used) < 2) {
    cat("  Skipping: trait observed in <2 populations.\n")
    next
  }
  if (length(unique(ped_sub$country)) < 2) {
    cat("  Skipping: trait observed in only one country.\n")
    next
  }

  Y <- pheno_t$trait
  do_standardize <- !(all(abs(Y - round(Y)) < 1e-8) && length(unique(Y)) <= 10)
  if (do_standardize) Y <- (Y - mean(Y)) / sd(Y)

  dat <- data.frame(
    Y       = Y,
    pop     = factor(as.character(ped_sub$pop), levels = pops_used),
    country = factor(ped_sub$country, levels = countries)
  )
  stopifnot(all(M_sub[outer(dat$pop, dat$pop, "!=")] == 0))

  # rotate each population block to the eigenbasis of M; w = U' 1 carries the intercept and pop effect
  rdat <- do.call(rbind, lapply(split(seq_len(nrow(dat)), dat$pop, drop = TRUE), function(i) {
    e <- eigen(M_sub[i, i], symmetric = TRUE)
    data.frame(Y = drop(crossprod(e$vectors, dat$Y[i])), w = colSums(e$vectors),
               lambda = e$values, pop = dat$pop[i[1]], country = dat$country[i[1]])
  }))
  rdat$wChina <- rdat$w * (rdat$country == "China")

  fit <- tryCatch(
    brm(
      bf(Y ~ 0 + w + wChina + (0 + w | gr(pop, cov = two.Theta.P)),
         nlf(sigma ~ sqrt(sA^2 * lambda + sE^2)),
         sA + sE ~ 1),
      data   = rdat,
      data2  = list(two.Theta.P = 2 * ThetaP_sub),
      family = brmsfamily("gaussian", link_sigma = "identity"),
      prior  = c(prior(exponential(1), class = sd),
                 prior(exponential(1), nlpar = sA, lb = 0),
                 prior(exponential(1), nlpar = sE, lb = 0)),
      iter = 12000, warmup = 4000, thin = 2,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    ),
    error = function(e) {
      cat("  brm failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) next
  print(summary(fit))

  # same column names as monet()$sampling in logAV_country_analysis.r
  d <- posterior::as_draws_df(fit)
  draws <- data.frame(
    b_Intercept       = d$b_w,
    b_countryChina    = d$b_wChina,
    sd_pop__Intercept = d$sd_pop__w,
    sd_ind__Intercept = d$b_sA_Intercept,
    sigma             = d$b_sE_Intercept
  )
  draws$V_AB      <- draws$sd_pop__Intercept^2
  draws$V_AW      <- draws$sd_ind__Intercept^2
  draws$log_ratio <- log(draws$V_AB / draws$V_AW)
  saveRDS(draws, file.path(outdir, sprintf("monet_%s_%s.rds", experiment, tid)))

  sampler <- brms::nuts_params(fit)
  lr <- post_summary(draws$log_ratio)
  cb <- post_summary(draws$b_countryChina)
  row_df <- data.frame(
    experiment        = experiment,
    trait             = tid,
    n_individuals     = nrow(dat),
    n_populations     = length(pops_used),
    log_ratio_mean    = lr[["mean"]],
    log_ratio_median  = lr[["median"]],
    log_ratio_lo      = lr[["lo"]],
    log_ratio_hi      = lr[["hi"]],
    log_ratio_p       = lr[["p"]],
    country_beta_mean = cb[["mean"]],
    country_beta_median = cb[["median"]],
    country_beta_lo   = cb[["lo"]],
    country_beta_hi   = cb[["hi"]],
    country_p         = cb[["p"]],
    n_divergent       = sum(sampler$Value[sampler$Parameter == "divergent__"]),
    max_rhat          = max(brms::rhat(fit), na.rm = TRUE),
    row.names         = NULL
  )
  summary_rows[[length(summary_rows) + 1]] <- row_df

  trait_csv <- file.path(outdir, sprintf("summary_logAV_country_%s_%s.csv", experiment, tid))
  write.csv(row_df, trait_csv, row.names = FALSE)
  cat("  Wrote:", trait_csv, "\n")
  cat(sprintf("  log_ratio = %.3f  [%.3f, %.3f]  p = %.4f\n",
              lr[["mean"]], lr[["lo"]], lr[["hi"]], lr[["p"]]))
  cat(sprintf("  country beta = %.3f  [%.3f, %.3f]  p = %.4f\n\n",
              cb[["mean"]], cb[["lo"]], cb[["hi"]], cb[["p"]]))
}

if (!nzchar(trait_filter) && length(summary_rows) > 0) {
  out_csv <- file.path(outdir, sprintf("summary_logAV_country_%s.csv", experiment))
  write.csv(do.call(rbind, summary_rows), out_csv, row.names = FALSE)
  cat("Wrote combined summary:", out_csv, "\n")
}
