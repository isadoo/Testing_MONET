#' @title monet: Log ratios of ancestral variances
#'
#' @description Given coancestry matrices for between and within populations and a trait data frame, 
#' this function estimates the log ratio of ancestral variances using a Bayesian mixed-effects model.
#' 
#' @usage monet(Theta.P, M, trait_dataframe, column_individual = "id", column_trait = "trait", 
#'             column_population = "population", column_se = NULL, formula_covariates = NULL, 
#'             iter = 5000, warmup = 2000, thin = 2, save_full_model = FALSE,
#'             standardize_trait = NULL, ...)
#'
#' @param Theta.P A square matrix representing the coancestry matrix between populations.
#'
#' @param M A square matrix representing the kinship-based relatedness matrix for individuals.
#'
#' @param trait_dataframe A data frame containing individual IDs and trait values. 
#' The first column should be individual IDs, and the second column should be trait values.
#'
#' @param column_individual The name of the column containing individual IDs. Default is "id".
#'
#' @param column_trait The name of the column containing trait values. Default is "trait".
#'
#' @param column_population The name of the column containing population IDs. Default is "population".
#'
#' @param column_se The name of the column containing standard errors of trait values. Default is NULL (no SEs).
#' If provided, the function will incorporate measurement error in the model.
#' 
#' @param formula_covariates A character string specifying additional covariates to include in the model.
#' For example, "age + sex" would add age and sex as fixed effects. Default is NULL (no covariates).
#' 
#' @param iter Number of MCMC iterations per chain. Default is 5000.
#' 
#' @param warmup Number of warmup (burn-in) iterations per chain. Default is 2000.
#' 
#' @param thin Thinning rate for the MCMC sampler. Default is 2 (every 2nd sample is kept).
#' 
#' @param save_full_model Logical. If TRUE, saves the full brms model object in the results.
#' If FALSE (default), only minimal posterior samples (fixed effects, variance components, and log-ratio) are retained to save memory.
#'
#' @param standardize_trait Logical or NULL. Controls z-score standardization of the response trait.
#' If TRUE, always standardize (unless variance is zero/non-finite). If FALSE, does NOT standardize.
#' If NULL (default), the function auto-detects likely discrete traits (integer-like with <= 10 unique values)
#' and skips standardization for those traits.
#' 
#' @param ... Additional arguments passed to the brms function.
#'
#' @return A monet type object containing:
#' \item{sampling}{Either the full brms model object (if save_full_model = TRUE) or a data frame with minimal posterior samples including fixed effects, variance components (V_AB, V_AW), and log_ratio.}
#' \item{log_ratio}{A list containing: p_value (two-tailed test that log-ratio differs from 0), mean (mean of log-ratio posterior), median (median of log-ratio posterior), ci_lower and ci_upper (95% credible interval bounds).}
#' \item{covariate_p_values}{Named numeric vector with two-tailed posterior sign-probability p-values for fixed-effect covariates (excluding intercept). NULL when no covariates are present.}
#' \item{hypothesis}{Results from brms hypothesis test comparing population vs individual variance.}
#' \item{trait_name}{Name of the trait column analyzed.}
#' \item{formula_used}{The model formula used in the analysis.}
#' \item{convergence}{A list containing n_divergent (number of divergent transitions) and max_rhat (maximum R-hat value across parameters).}
#' 
#' @details This function standardizes trait data, constructs a Bayesian mixed-effects model 
#' using `brms`, and estimates ancestral variances. The function assumes no crosses between populations 
#' and analyzes one trait at a time.
#'
#' @author Isabela do O \email{isabela.doo@@unil.ch}
#'
#' @references 
#' 
#' - Goudet & Weir (2023)
#' - do O et al (2025)
#'
#' @export
monet <- function(Theta.P, 
                 M, 
                 trait_dataframe, 
                 column_individual = "id", 
                 column_trait = "trait", 
                 column_population = "population",
                 column_se = NULL,
                 formula_covariates = NULL,
                 iter = 5000, warmup = 2000, thin = 2,
                 save_full_model = FALSE,
                 standardize_trait = NULL,
                 ...) {
  
  #check input types and dimensions ------------------------
  if (!is.matrix(Theta.P) || !is.matrix(M)) {
    stop("Theta.P and M must be matrices.")
  }
  
  if (!is.data.frame(trait_dataframe) || ncol(trait_dataframe) < 2) {
    stop("trait_dataframe must be a data frame with at least two columns (ID and trait values).")
  }
  #------------------------------------------------------------
  
  #extract columns by name if provided ------------------------------------
  id_col <- if(is.numeric(column_individual)) column_individual else which(names(trait_dataframe) == column_individual)
  trait_col <- if(is.numeric(column_trait)) column_trait else which(names(trait_dataframe) == column_trait)
  
  #Identify populations per individual
  population_blocks_df <- counting_blocks_matrix(M) #this function counts the number of blocks of non-zero rows in a matrix
  individuals_per_population_F1 <- population_blocks_df$rows
  number_of_blocks <- length(population_blocks_df$block) 
  pop_ids <- rep(1:number_of_blocks, individuals_per_population_F1[1:number_of_blocks]) 
  number_populations <- nrow(Theta.P)
  
  if (number_of_blocks != number_populations & !(column_population %in% names(trait_dataframe))) {
    warning(paste0("Mismatch between detected groups based on M matrix (", number_of_blocks,") and populations in Theta.P dimensions (", number_populations, ")\n
                   We have have no way of knowing what is the correct number of individuals in each subpopulation."))
  }
  
  Y_raw <- trait_dataframe[, trait_col]
  
  if (is.list(Y_raw)) {
    Y_raw <- unlist(Y_raw)
  }

  # Keep original values for coercion diagnostics
  Y_raw_chr <- as.character(Y_raw)
  
  # Ensure trait is numeric for current Gaussian-style model.
  # If non-numeric values exist, fail early with a clear message.
  Y <- suppressWarnings(as.numeric(Y_raw_chr))
  coercion_failed <- !is.na(Y_raw_chr) & nzchar(trimws(Y_raw_chr)) & is.na(Y)
  if (any(coercion_failed)) {
    bad_values <- unique(Y_raw_chr[coercion_failed])
    bad_preview <- paste(utils::head(bad_values, 5), collapse = ", ")
    stop(
      paste0(
        "Trait column contains non-numeric values and cannot be modeled with the current Gaussian response. ",
        "Example offending values: ", bad_preview, "."
      )
    )
  }
  
  #remove any NAs
  valid_indices <- !is.na(Y)
  Y <- Y[valid_indices]
  
  #Filter the dataframe to match
  trait_dataframe <- trait_dataframe[valid_indices, ]
  
  if (length(Y) < 2) {
    stop("Not enough non-missing numeric trait values after filtering.")
  }

  var_Y <- var(Y)
  n_unique_y <- length(unique(Y))
  integer_like <- all(abs(Y - round(Y)) < 1e-8)
    # Validate population column against Theta.P rownames (if population column exists in dataframe)
  if (column_population %in% names(trait_dataframe)) {
    pops_in_data <- unique(na.omit(trait_dataframe[[column_population]]))
    pops_in_theta <- rownames(Theta.P)

    missing_pops <- setdiff(pops_in_data, pops_in_theta)
    if (length(missing_pops) > 0) {
      stop(
        paste0(
          "The following population IDs in column '", column_population,
          "' are not found in rownames(Theta.P): ",
          paste(missing_pops, collapse = ", ")
        )
      )
    }
  }

  # Determine whether to standardize
  if (is.null(standardize_trait)) {
    # Auto mode: skip standardization for likely discrete traits.
    do_standardize <- !(integer_like && n_unique_y <= 10)
    if (!do_standardize) {
      warning(
        "Trait appears discrete (integer-like with <= 10 unique values); skipping standardization in auto mode. ",
        "Set standardize_trait = TRUE to force z-score standardization."
      )
    }
  } else if (is.logical(standardize_trait) && length(standardize_trait) == 1 && !is.na(standardize_trait)) {
    do_standardize <- standardize_trait
  } else {
    stop("standardize_trait must be TRUE, FALSE, or NULL.")
  }

  if (do_standardize) {
    if (!is.finite(var_Y) || var_Y <= 0) {
      warning("Trait variance is non-finite or zero; skipping standardization.")
      do_standardize <- FALSE
    }
  }

  if (do_standardize) {
    Y <- Y - mean(Y)
    Y <- Y / sqrt(var_Y)
  }
  
  have_se <- !is.null(column_se) && (column_se %in% names(trait_dataframe))
  #Incorportate measurment errors (se) if provided
  if (have_se) {
    Y_se <- as.numeric(trait_dataframe[[column_se]])
    if (do_standardize) {
      Y_se <- Y_se / sqrt(var_Y)
    }
  }


  #From VB = VA*2FST
  two.Theta.P <- 2 * Theta.P

  ind_col <- trait_dataframe[,id_col]
  if (is.list(ind_col)) {
    ind_col <- unlist(ind_col)
  }
  ind_col <- as.character(ind_col)

  # Build an individual -> population lookup from M block structure.
  # If block count matches Theta.P dimension, use Theta.P rownames as labels.
  m_pop_lookup <- NULL
  if (!is.null(rownames(M))) {
    if (number_of_blocks == number_populations && !is.null(rownames(Theta.P))) {
      block_labels <- rownames(Theta.P)[seq_len(number_of_blocks)]
    } else {
      block_labels <- paste0("pop_", seq_len(number_of_blocks))
    }

    m_pop_labels <- rep(block_labels, individuals_per_population_F1[seq_len(number_of_blocks)])
    if (length(m_pop_labels) == nrow(M)) {
      names(m_pop_labels) <- rownames(M)
      m_pop_lookup <- m_pop_labels
    }
  }

  if (column_population %in% names(trait_dataframe)) {
    pop_col <- trait_dataframe[,column_population]
    if (is.list(pop_col)) {
      pop_col <- unlist(pop_col)
    }
    pop_col <- as.character(pop_col)

    overlap_with_thetap <- sum(pop_col %in% rownames(Theta.P), na.rm = TRUE)
    if (overlap_with_thetap == 0 && !is.null(m_pop_lookup)) {
      warning(
        sprintf(
          "column_population='%s' has zero overlap with rownames(Theta.P); using M-based population mapping from individual IDs instead.",
          column_population
        )
      )
      pop_col <- unname(m_pop_lookup[ind_col])
    }
  } else if (!is.null(m_pop_lookup)) {
    pop_col <- unname(m_pop_lookup[ind_col])
  } else if (length(pop_ids) == length(ind_col)) {
    pop_col <- paste0("pop_", pop_ids)
  } else {
    stop("Could not determine population labels from column_population or M block mapping.")
  }

  if (is.list(pop_col)) {
    pop_col <- unlist(pop_col)
  }
  pop_col <- as.character(pop_col)

  #Build the data frame
  dat <- data.frame(pop = pop_col, ind = ind_col, Y = Y)
  if (have_se) dat$Y_se <- Y_se

  # Align grouping levels with covariance matrix rownames.
  # brms requires the factor levels used in gr() to match the covariance matrix names exactly.
  if (is.null(rownames(M)) || is.null(rownames(Theta.P))) {
    stop("Theta.P and M must have rownames that identify grouping levels.")
  }

  keep_ind <- dat$ind %in% rownames(M)
  if (any(!keep_ind)) {
    warning(sprintf(
      "%d individuals were not represented in rownames(M), so they are not considered in the analysis.",
      sum(!keep_ind)
    ))
    dat <- dat[keep_ind, , drop = FALSE]
    trait_dataframe <- trait_dataframe[keep_ind, , drop = FALSE]
  }

  keep_pop <- dat$pop %in% rownames(Theta.P)
  if (any(!keep_pop)) {
    warning(sprintf(
      "Dropping %d rows whose pop values are not present in rownames(Theta.P).",
      sum(!keep_pop)
    ))
    dat <- dat[keep_pop, , drop = FALSE]
    trait_dataframe <- trait_dataframe[keep_pop, , drop = FALSE]
  }

  dat$ind <- factor(dat$ind, levels = rownames(M))
  dat$pop <- factor(dat$pop, levels = rownames(Theta.P))

  # Add any additional covariates that might be specified
  if (!is.null(formula_covariates)) {
    # Extract real variable names from formula terms, including random-effect grouping vars
    covariate_names <- unique(all.vars(stats::as.formula(paste("~", formula_covariates))))
    
    for (cov_name in covariate_names) {
      if (cov_name %in% names(trait_dataframe)) {
        cov_data <- trait_dataframe[, cov_name]
        if (is.list(cov_data)) cov_data <- unlist(cov_data)
        dat[[cov_name]] <- cov_data
      } else {
        stop(paste("Covariate", cov_name, "not found in trait_dataframe"))
      }
    }
  }
  
  #Build the complete formula
  base_formula <- "Y ~ 1 + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = M))"
  base_rhs <- "(1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = M))"
  rhs <- if (is.null(formula_covariates)) base_rhs else paste0(formula_covariates, " + ", base_rhs)

  if (have_se) {
    # Use measurement error on the response; still estimate residual sigma
    formula_string <- paste0("Y | se(Y_se, sigma = TRUE) ~ 1 + ", rhs)
    model_formula <- brms::bf(as.formula(formula_string))
    if (do_standardize) {
      cat("Measurement SE column: ", column_se, " (scaled to standardized Y)\n", sep = "")
    } else {
      cat("Measurement SE column: ", column_se, " (not scaled; trait not standardized)\n", sep = "")
    }
  } else {
    formula_string <- paste0("Y ~ 1 + ", rhs)
    model_formula <- as.formula(formula_string)
  }
  
  cat("Using formula:", formula_string, "\n")
  
  #Diagnostics
  n_total <- nrow(dat)
  n_complete <- sum(complete.cases(dat))
  cat("Total rows:", n_total, " ; complete rows used:", n_complete, "\n")
  
  #check for NAs
  na_counts <- colSums(is.na(dat))
  if (any(na_counts > 0)) {
    cat("Warning: NAs found in columns:\n")
    print(na_counts[na_counts > 0])
  }
  
  
  #Bayesian model - using brms package
  #Using tryCatch to handle convergence issues 
  brms_mf <- tryCatch({
    brm(
      formula = model_formula,
      data = dat,
      data2 = list(two.Theta.P = two.Theta.P, M = M), 
      iter = iter, warmup = warmup, thin = thin,
      ...
    )

  }, error = function(e) {
    cat("Error fitting model:", e$message, "\n")
    stop(e)
  })
  
  # Print summary properly
  cat("\n---- Model Summary ----\n")
  print(summary(brms_mf))
  cat("\n")
  
  # Check for convergence issues
  rhats <- brms::rhat(brms_mf)
  if (any(rhats > 1.01, na.rm = TRUE)) {
    warning("Some Rhat values > 1.01, indicating potential convergence issues")
  }

  
  
  
  # Extract sampling parameters directly from the fitted model
  iter   <- brms_mf$fit@sim$iter      # total iterations (including warmup)
  warmup <- brms_mf$fit@sim$warmup    # warmup iterations
  chains <- brms_mf$fit@sim$chains    # number of chains
  thin   <- brms_mf$fit@sim$thin      # thinning interval
  post_warmup_transitions <- (iter - warmup) * chains
  threshold_1pct <- 0.01 * post_warmup_transitions
  sampler_params <- brms::nuts_params(brms_mf)
  # get divergent transitions and check if its under 1%
  n_divergent    <- sum(subset(sampler_params, Parameter == "divergent__")$Value)
  if (n_divergent > 0) {
    warning(paste("Model had", n_divergent, "divergent transitions after warmup"))
  }

  #variance components
  var_components <- lapply(VarCorr(brms_mf, summary = FALSE), function(x) x$sd^2)
  var_df <- as.data.frame(do.call(cbind, var_components))

  hyp <- "sd_pop__Intercept^2 - sd_ind__Intercept^2 = 0"
  the_hyp <- hypothesis(brms_mf, hyp, class = NULL)

  # --- core draws we need ---
  # Get all draws once, then select columns we care about depending on save_full_model
  all_draws <- posterior::as_draws_df(brms_mf)

fe_cols <- grep("^b_", names(all_draws), value = TRUE)     # fixed effects
sd_cols <- c("sd_pop__Intercept", "sd_ind__Intercept")     # two RE SDs

have_sd <- sd_cols[sd_cols %in% names(all_draws)]

# Covariate p-values from posterior sign probability (exclude intercept)
covariate_p_values <- NULL
cov_cols <- setdiff(fe_cols, "b_Intercept")
if (length(cov_cols) > 0) {
  covariate_p_values <- sapply(cov_cols, function(col) {
    x <- all_draws[[col]]
    2 * min(mean(x >= 0, na.rm = TRUE), mean(x <= 0, na.rm = TRUE))
  })
  names(covariate_p_values) <- sub("^b_", "", cov_cols)
}

if (length(have_sd) == 2) {
  minimal_samples <- all_draws[, unique(c(fe_cols, have_sd)), drop = FALSE]
  minimal_samples$V_AB  <- minimal_samples$sd_pop__Intercept^2
  minimal_samples$V_AW  <- minimal_samples$sd_ind__Intercept^2
  minimal_samples$log_ratio <- log(minimal_samples$V_AB / minimal_samples$V_AW)
} else {
  warning("Issue with your model - missing sd_pop__Intercept and sd_ind__Intercept")
}

  # For summary stats below we still use post_samples with log_ratio
  post_samples <- minimal_samples

  # --- summaries of log-ratio ---
  if ("log_ratio" %in% names(post_samples)) {
    quant_log_med   <- stats::quantile(post_samples$log_ratio, c(0.5, 0.025, 0.975))
    mean_log_ratio  <- mean(post_samples$log_ratio)
    quant_log_ratio <- stats::quantile(post_samples$log_ratio, probs = c(0.025, 0.975))
    p_value <- 2 * mean(sign(post_samples$log_ratio) != sign(stats::median(post_samples$log_ratio)))
  } else {
    warning("No log_ratio samples found - cannot compute log-ratio summaries.")
  }

  # ------------------------------------------------------------------------------------
   # ~~~~~~~~~~~~   # preparing results object # ~~~~~~~~~~~~~~~~~~
  # ------------------------------------------------------------------------------------
  results <- list(
  sampling = if (isTRUE(save_full_model)) {
    # Return the full brms model as 'sampling'
    brms_mf
  } else {
    # Minimal sampling: fixed effects + two RE SDs + var_* + log_ratio
    minimal_samples
  },

  log_ratio = list(
    p_value = p_value,
    mean = mean_log_ratio,
    median = quant_log_med["50%"],
    ci_lower = quant_log_ratio[1],
    ci_upper = quant_log_ratio[2]
  ),

  covariate_p_values = covariate_p_values,

  hypothesis = the_hyp$hypothesis[2:5],
  trait_name = names(trait_dataframe)[trait_col],
  formula_used = formula_string,
  convergence = list(
    n_divergent = n_divergent,
    max_rhat = max(rhats, na.rm = TRUE)
  )
)
class(results) <- "monet"
return(results)
}
#' @export
print.monet <- function(x, ...) {
  cat("\n---- MONET Analysis Results ----\n\n")
  cat("Trait analyzed:", x$trait_name, "\n")
  cat("Formula used:", x$formula_used, "\n\n")
  
  cat("--- Log-Ratio of Ancestral Variances ---\n")
  cat(sprintf("  Mean: %.4f\n", x$log_ratio$mean))
  cat(sprintf("  Median: %.4f\n", x$log_ratio$median))
  cat(sprintf("  95%% Credible Interval: [%.4f, %.4f]\n", 
              x$log_ratio$ci_lower, x$log_ratio$ci_upper))
  cat(sprintf("  P-value: %.4f\n\n", x$log_ratio$p_value))
  
  cat("--- Convergence Diagnostics ---\n")
  cat(sprintf("  Max R-hat: %.4f", x$convergence$max_rhat))
  if (x$convergence$max_rhat > 1.01) {
    cat(" (WARNING: > 1.01)\n")
  }
  
  cat(sprintf("  Divergent transitions: %d", x$convergence$n_divergent))
  if (x$convergence$n_divergent > 0) {
    cat(" (WARNING: Consider increasing adapt_delta)\n")
  }
  
  cat("\nUse plot() to visualize the posterior distribution.\n\n")
  invisible(x)
}

#' @export
summary.monet <- function(object, ...) {
  cat("\n---- MONET Summary ----\n\n")
  cat(sprintf("Trait: %s\n\n", object$trait_name))
  
  cat("Log-Ratio of Ancestral Variances (log(VB/VA)):\n")
  cat(sprintf("  Mean:   %7.4f\n", object$log_ratio$mean))
  cat(sprintf("  Median: %7.4f\n", object$log_ratio$median))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
              object$log_ratio$ci_lower, object$log_ratio$ci_upper))
  cat(sprintf("  p-value: %.4f %s\n", 
              object$log_ratio$p_value,
              ifelse(object$log_ratio$p_value < 0.05, "*", "")))

  if (!is.null(object$covariate_p_values) && length(object$covariate_p_values) > 0) {
    cat("\nCovariate posterior p-values:\n")
    for (nm in names(object$covariate_p_values)) {
      pv <- object$covariate_p_values[[nm]]
      cat(sprintf("  %s: %.4g\n", nm, pv))
    }
  }

  cat("\n")
  invisible(object)
}

#' @export
plot.monet <- function(x, ...) {
  # Accept either a brmsfit (full model) or a draws data.frame in x$sampling
  if (inherits(x$sampling, "brmsfit")) {
    # compute log_ratio from the model draws
    draws <- posterior::as_draws_df(x$sampling)
    if (!all(c("sd_pop__Intercept", "sd_ind__Intercept") %in% names(draws))) {
      stop("Could not find sd_pop__Intercept/sd_ind__Intercept in model draws to compute log_ratio.")
    }
    lr <- log((draws$sd_pop__Intercept^2) / (draws$sd_ind__Intercept^2))
  } else {
    samp <- x$sampling
    if (is.null(samp) || !"log_ratio" %in% names(samp)) {
      stop("No 'log_ratio' samples found.")
    }
    lr <- samp$log_ratio
  }

  med <- stats::median(lr)
  ci  <- stats::quantile(lr, c(0.025, 0.975))

  
  d <- stats::density(lr)
  plot(d, main = "Posterior of log-ratio",
    xlab = "LogAV", ylab = "Posterior density")
    abline(v = med, lty = 2)
    abline(v = ci, lty = 3)
    legend("topright",
        legend = c(paste0("median = ", round(med, 3)),
                    paste0("95% CI [", round(ci[1], 3), ", ", round(ci[2], 3), "]")),
        lty = c(2, 3), bty = "n")

  invisible(NULL)
}