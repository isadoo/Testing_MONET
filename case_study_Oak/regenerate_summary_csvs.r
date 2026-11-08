# =====================================================================
# Regenerate summary CSV files from current .rds results
# This ensures the CSV files match the current state of .rds files
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
})

results_dir <- "results_logAV_env"

# All species x garden combinations
configs <- data.frame(
  species = c("faginea", "faginea", "lusitanica"),
  garden  = c("GH", "OC", "GH"),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(configs)) {
  species <- configs$species[i]
  garden  <- configs$garden[i]
  
  cat("Processing", species, garden, "...\n")
  
  # Find all .rds files for this configuration
  pattern <- sprintf("^monet_%s_%s_.*\\.rds$", species, garden)
  rds_files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
  
  if (length(rds_files) == 0) {
    cat("  No results found, skipping.\n\n")
    next
  }
  
  summary_rows <- list()
  
  for (rds_path in rds_files) {
    # Extract trait name
    fname <- basename(rds_path)
    trait <- sub(sprintf("^monet_%s_%s_(.+)\\.rds$", species, garden), "\\1", fname)
    
    # Load result
    fit <- tryCatch(readRDS(rds_path), error = function(e) {
      cat("  Error loading", fname, ":", e$message, "\n")
      return(NULL)
    })
    if (is.null(fit)) next
    
    # Extract covariate p-values
    cov_p <- fit$covariate_p_values
    pc1_p <- if (!is.null(cov_p) && "b_PC_climate1" %in% names(cov_p))
               cov_p[["b_PC_climate1"]] else NA_real_
    pc2_p <- if (!is.null(cov_p) && "b_PC_climate2" %in% names(cov_p))
               cov_p[["b_PC_climate2"]] else NA_real_
    
    # Posterior means/CIs of betas
    beta_summary <- function(name) {
      brms_fit <- fit$sampling
      if (is.null(brms_fit) || !inherits(brms_fit, "brmsfit")) {
        return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
      }
      tryCatch({
        draws <- as.matrix(brms_fit, variable = name)
        if (ncol(draws) > 0) {
          v <- draws[, 1]
          c(mean = mean(v),
            lo   = unname(quantile(v, 0.025)),
            hi   = unname(quantile(v, 0.975)))
        } else {
          c(mean = NA_real_, lo = NA_real_, hi = NA_real_)
        }
      }, error = function(e) {
        c(mean = NA_real_, lo = NA_real_, hi = NA_real_)
      })
    }
    b1 <- beta_summary("b_PC_climate1")
    b2 <- beta_summary("b_PC_climate2")
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      species          = species,
      garden           = garden,
      trait            = trait,
      n_individuals    = fit$n_individuals,
      n_populations    = fit$n_populations,
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
  }
  
  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    # Sort by trait name
    summary_df <- summary_df[order(summary_df$trait), ]
    
    out_csv <- file.path(results_dir,
                         sprintf("summary_logAV_env_%s_%s.csv", species, garden))
    write.csv(summary_df, out_csv, row.names = FALSE)
    cat("  Wrote", nrow(summary_df), "traits to", basename(out_csv), "\n")
  } else {
    cat("  No valid results to summarize.\n")
  }
  
  cat("\n")
}

cat("Done!\n")
