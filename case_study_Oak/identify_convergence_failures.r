# =====================================================================
# Identify traits with poor convergence diagnostics
# Checks both max_rhat and min_ess from saved monet .rds files
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
})

# Configuration
results_dir <- "results_logAV_env"
rhat_threshold <- 1.05
ess_threshold <- 300

# Find all monet .rds files
rds_files <- list.files(results_dir, pattern = "^monet_.*\\.rds$", full.names = TRUE)

cat("Checking", length(rds_files), "monet results for convergence issues...\n\n")
cat("Criteria: max_rhat >", rhat_threshold, "OR min_ess <", ess_threshold, "\n\n")

bad_traits <- list()

for (rds_path in rds_files) {
  # Parse filename: monet_species_garden_trait.rds
  fname <- basename(rds_path)
  parts <- sub("^monet_(.*)_([A-Z]+)_(.+)\\.rds$", "\\1|\\2|\\3", fname)
  if (!grepl("\\|", parts)) next
  
  parts_split <- strsplit(parts, "\\|")[[1]]
  species <- parts_split[1]
  garden  <- parts_split[2]
  trait   <- parts_split[3]
  
  # Load result
  fit <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(fit)) next
  
  # Extract diagnostics
  max_rhat <- fit$convergence$max_rhat
  
  # Get min ESS from the brms samples
  # The fit$sampling contains posterior draws; check for ESS
  min_ess <- NA_real_
  if (!is.null(fit$full_model)) {
    # If full brms model was saved
    brms_summary <- summary(fit$full_model)
    if (!is.null(brms_summary$fixed)) {
      ess_vals <- c(
        if ("Bulk_ESS" %in% names(brms_summary$fixed)) brms_summary$fixed[, "Bulk_ESS"],
        if ("Tail_ESS" %in% names(brms_summary$fixed)) brms_summary$fixed[, "Tail_ESS"]
      )
      if (length(ess_vals) > 0) min_ess <- min(ess_vals, na.rm = TRUE)
    }
    # Also check random effects
    if (!is.null(brms_summary$random)) {
      for (re in brms_summary$random) {
        ess_vals <- c(
          if ("Bulk_ESS" %in% names(re)) re[, "Bulk_ESS"],
          if ("Tail_ESS" %in% names(re)) re[, "Tail_ESS"]
        )
        if (length(ess_vals) > 0) {
          re_min_ess <- min(ess_vals, na.rm = TRUE)
          if (is.na(min_ess) || re_min_ess < min_ess) min_ess <- re_min_ess
        }
      }
    }
  }
  
  # Check if problematic
  is_bad <- FALSE
  reason <- c()
  
  if (!is.na(max_rhat) && max_rhat > rhat_threshold) {
    is_bad <- TRUE
    reason <- c(reason, sprintf("Rhat=%.3f", max_rhat))
  }
  
  if (!is.na(min_ess) && min_ess < ess_threshold) {
    is_bad <- TRUE
    reason <- c(reason, sprintf("min_ESS=%.0f", min_ess))
  }
  
  if (is_bad) {
    bad_traits[[length(bad_traits) + 1]] <- data.frame(
      species = species,
      garden  = garden,
      trait   = trait,
      max_rhat = max_rhat,
      min_ess  = min_ess,
      reason   = paste(reason, collapse = ", "),
      stringsAsFactors = FALSE
    )
    cat(sprintf("✗ %s_%s_%s: %s\n", species, garden, trait, paste(reason, collapse = ", ")))
  }
}

cat("\n")

if (length(bad_traits) > 0) {
  bad_df <- do.call(rbind, bad_traits)
  
  # Save list
  out_csv <- file.path(results_dir, "traits_needing_rerun.csv")
  write.csv(bad_df, out_csv, row.names = FALSE)
  
  cat("Found", nrow(bad_df), "traits with convergence issues.\n")
  cat("Saved to:", out_csv, "\n\n")
  
  # Print summary by species/garden
  cat("Summary by configuration:\n")
  agg <- aggregate(trait ~ species + garden, data = bad_df, FUN = length)
  names(agg)[3] <- "n_bad_traits"
  print(agg, row.names = FALSE)
  
} else {
  cat("All traits converged successfully!\n")
}
