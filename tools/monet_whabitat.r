###########################################################################
#' @title monet_whabitat: Log ratios of ancestral variances with habitat information
#'
#' @description Given coancestry matrices for between and within populations, 
#' trait data frame, and habitat optima, this function estimates the log ratio 
#' of ancestral variances using a Bayesian mixed-effects model that includes 
#' environmental optima as a fixed effect.
#' 
#' @usage monet_whabitat(Theta.P, The.M, trait_dataframe, habitat_optima, 
#'                      column_individual, column_trait)
#'
#' @param Theta.P A square matrix representing the coancestry matrix between populations.
#'
#' @param The.M A square matrix representing the kinship-based relatedness matrix for individuals.
#'
#' @param trait_dataframe A data frame containing individual IDs and trait values.
#'
#' @param habitat_optima A numeric vector of environmental optima for each population.
#'
#' @param column_individual Column name or index for individual IDs (default: "id").
#'
#' @param column_trait Column name or index for trait values (default: "trait").
#'
#' @return A monet_whabitat type object containing:
#' \item{posterior_samples}{A list with posterior samples of variance components and residuals.}
#' \item{BRMS_stats}{A list with the median, lower and upper bounds of the 95% credible interval for the variance difference.}
#' \item{log_ratio}{A list with the p-value of the log ratios, the mean of the log ratio of between- and within-population variance, and confidence intervals.}
#' \item{habitat_analysis}{A list with habitat coefficient estimates, confidence intervals, and p-value.}
#' \item{BRMS_model}{A list with the fitted Bayesian model and hypothesis test results.}
#'
#' @details This function extends the standard LAVA analysis by including habitat 
#' optima as a fixed effect to test for local adaptation. The habitat coefficient 
#' tests whether populations from different environmental optima show systematic 
#' genetic differentiation in trait values.
#'
#' @author Isabela do O \email{isabela.doo@@unil.ch}
#'
#' @references 
#' - Goudet & Weir (2023)
#' - do O et al (2025)
#'
monet_whabitat <- function(Theta.P, 
                          The.M, 
                          trait_dataframe,
                          habitat_optima,
                          column_individual = "id", 
                          column_trait = "trait", ...) {
  
  # Check input types and dimensions ------------------------
  if (!is.matrix(Theta.P) || !is.matrix(The.M)) {
    stop("Theta.P and The.M must be matrices.")
  }
  
  if (!is.data.frame(trait_dataframe) || ncol(trait_dataframe) < 2) {
    stop("trait_dataframe must be a data frame with at least two columns (ID and trait values).")
  }
  
  if (!is.numeric(habitat_optima)) {
    stop("habitat_optima must be a numeric vector.")
  }
  
  # Extract columns by name if provided ------------------------------------
  id_col <- if(is.numeric(column_individual)) column_individual else which(names(trait_dataframe) == column_individual)
  trait_col <- if(is.numeric(column_trait)) column_trait else which(names(trait_dataframe) == column_trait)

  # Identify populations per individual
  source("isaChapter2/tools/counting_blocks_matrix.r")
  population_blocks_df <- counting_blocks_matrix(The.M)
  individuals_per_population_F1 <- population_blocks_df$rows
  number_of_blocks <- length(population_blocks_df$block) 
  pop_ids <- rep(1:number_of_blocks, individuals_per_population_F1[1:number_of_blocks]) 
  number_populations <- nrow(Theta.P)
  
  # Check habitat_optima length
  if (length(habitat_optima) != number_populations) {
    stop(paste0("habitat_optima length (", length(habitat_optima), 
                ") must match number of populations (", number_populations, ")."))
  }
  
  if (number_of_blocks != number_populations) {
    warning(paste0("Mismatch between detected populations based on The.M matrix (", 
                   number_of_blocks,") and Theta.P dimensions (", number_populations, ")."))
  }
  
  # Standardize trait data
  Y <- trait_dataframe[,trait_col]
  Y <- Y - mean(Y)
  var_Y <- var(Y)
  Y <- Y / sqrt(var_Y)
  
  # Labels
  row.names(Theta.P) <- colnames(Theta.P) <- paste("pop", 1:number_populations, sep = "_")
  row.names(The.M) <- colnames(The.M) <- trait_dataframe[,id_col]
  
  # From VB = VA*2FST
  two.Theta.P <- 2 * Theta.P

  pop_labels <- paste0("pop_", pop_ids)
  
  # Assign habitat optima to individuals based on their population
  individual_optima <- habitat_optima[pop_ids]
  
  # Standardize optima for better MCMC sampling
  optima_scaled <- scale(individual_optima)[,1]
  
  # Create data frame
  dat <- data.frame(
    pop = pop_labels, 
    ind = trait_dataframe[,id_col], 
    Y = Y,
    optima_scaled = optima_scaled
  )
  
  # Bayesian model with habitat as fixed effect
  cat("Fitting Bayesian model with habitat information...\n")
  brms_mf <- brm(Y ~ optima_scaled + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = The.M)), 
                 data = dat, 
                 data2 = list(two.Theta.P = two.Theta.P, The.M = The.M), 
                 family = gaussian(), 
                 chains = 8, 
                 cores = 4, 
                 iter = 3000, 
                 warmup = 1000, 
                 thin = 2, 
                 ...)
  
  # Variance components
  var_components <- lapply(VarCorr(brms_mf, summary = FALSE), function(x) x$sd^2)
  var_df <- as.data.frame(do.call(cbind, var_components))
  
  quant_med <- quantile(var_df$pop - var_df$ind, c(0.5, 0.025, 0.975))
  mean_diff <- mean(var_df$pop - var_df$ind)
  
  # Hypothesis testing for variance difference
  hyp <- "sd_pop__Intercept^2 - sd_ind__Intercept^2 = 0"
  the_hyp <- hypothesis(brms_mf, hyp, class = NULL)
  
  # Posteriors: VA,B and VA,A - estimated ancestral variances
  post_samples <- as_draws_df(brms_mf, variable = c("sd_pop__Intercept", "sd_ind__Intercept", "b_optima_scaled"))
  post_samples$log_ratio <- log(post_samples$sd_pop__Intercept^2 / post_samples$sd_ind__Intercept^2)
  
  mean_log_ratio <- mean(post_samples$log_ratio)
  quant_log_ratio <- quantile(post_samples$log_ratio, probs = c(0.025, 0.975))
  
  p_value <- 2 * mean(sign(post_samples$log_ratio) != sign(median(post_samples$log_ratio)))
  
  # Habitat-specific analysis
  optima_coefficient_mean <- mean(post_samples$b_optima_scaled)
  optima_coefficient_ci <- quantile(post_samples$b_optima_scaled, probs = c(0.025, 0.975))
  optima_p_value <- 2 * mean(sign(post_samples$b_optima_scaled) != sign(median(post_samples$b_optima_scaled)))
  
  # S3 object of class "monet_whabitat"
  results <- list(
    posteriors_samples = post_samples,
    
    # Basic statistics
    BRMS_stats = list(
      mean_diff = mean_diff,
      median_diff = quant_med[1],
      ci_lower_diff = quant_med[2],
      ci_upper_diff = quant_med[3]
    ),
    
    # Log ratio statistics
    log_ratio = list(
      mean_log_ratio = mean_log_ratio,
      log_ratio_ci_lower = quant_log_ratio[1],
      log_ratio_ci_upper = quant_log_ratio[2],
      p_value = p_value
    ),
    
    # Habitat analysis
    habitat_analysis = list(
      optima_coefficient_mean = optima_coefficient_mean,
      optima_coefficient_ci_lower = optima_coefficient_ci[1],
      optima_coefficient_ci_upper = optima_coefficient_ci[2],
      optima_p_value = optima_p_value
    ),
    
    # Model and hypothesis test results
    BRMS_model = list(
      model = brms_mf,
      hypothesis = the_hyp$hypothesis[2:5]
    ),
    
    trait_name = names(trait_dataframe)[trait_col],
    habitat_optima = habitat_optima
  )
  
  # Setting the class attribute to create an S3 object
  class(results) <- "monet_whabitat"
  
  return(results)
}

# Defining print method for the monet_whabitat S3 object
#' @export
#' @method print monet_whabitat
print.monet_whabitat <- function(x, ...) {
  # Header
  cat("\n===============================\n")
  cat("Log Ancestral Variance Analysis with Habitat (LAVA)\n")
  cat("===============================\n\n")
  
  # Key findings
  cat("Log ratio of estimated ancestral variances (between/within):\n")
  cat(sprintf("  Mean: %.4f (95%% CI: %.4f to %.4f)\n", 
              x$log_ratio$mean_log_ratio,
              x$log_ratio$log_ratio_ci_lower,
              x$log_ratio$log_ratio_ci_upper))
  cat(sprintf("  P-value: %.4f\n", x$log_ratio$p_value))
  
  cat("\nDifference in variance components (between/within):\n")
  cat(sprintf("  Mean: %.4f\n", x$BRMS_stats$mean_diff))
  cat(sprintf("  Median: %.4f (95%% CI: %.4f to %.4f)\n", 
              x$BRMS_stats$median_diff,
              x$BRMS_stats$ci_lower_diff,
              x$BRMS_stats$ci_upper_diff))
  
  cat("\nHabitat effect (test for local adaptation):\n")
  cat(sprintf("  Coefficient: %.4f (95%% CI: %.4f to %.4f)\n",
              x$habitat_analysis$optima_coefficient_mean,
              x$habitat_analysis$optima_coefficient_ci_lower,
              x$habitat_analysis$optima_coefficient_ci_upper))
  cat(sprintf("  P-value: %.4f\n", x$habitat_analysis$optima_p_value))
  
  cat("\n===============================\n")
  
  invisible(x)
}

# Plotting method for monet_whabitat S3 object
#' @export
#' @method plot monet_whabitat
plot.monet_whabitat <- function(x, which = "both", 
                               main_density = "Posterior Distribution of Log Ratio",
                               main_scatter = "Posterior Samples of Variance Components",
                               main_habitat = "Habitat Effect on Trait",
                               ...) {
  
  which <- match.arg(which, choices = c("both", "density", "scatter", "habitat", "all"))
  
  # Setting up the figure layout
  if (which == "all") {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    on.exit(par(old_par))
  } else if (which == "both") {
    old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    on.exit(par(old_par))
  }
  
  # Plot density function
  plot_density <- function() {
    density_lr <- density(x$posteriors_samples$log_ratio)
    plot(density_lr, main = main_density, 
         xlab = "Log Ratio of Ancestral Variances", ylab = "Density")
    abline(v = x$log_ratio$mean_log_ratio, col = "red", lwd = 2)
    abline(v = x$log_ratio$log_ratio_ci_lower, col = "red", lty = 2)
    abline(v = x$log_ratio$log_ratio_ci_upper, col = "red", lty = 2)
    abline(v = 0, col = "gray", lty = 2)
  }
  
  # Scatter function
  plot_scatter <- function() {
    plot(x$posteriors_samples$sd_pop__Intercept^2, 
         x$posteriors_samples$sd_ind__Intercept^2,
         xlab = "Between-population variance", 
         ylab = "Within-population variance",
         main = main_scatter,
         pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3))
    abline(a = 0, b = 1, col = "red", lwd = 2)
  }
  
  # Habitat effect plot
  plot_habitat <- function() {
    density_hab <- density(x$posteriors_samples$b_optima_scaled)
    plot(density_hab, main = main_habitat,
         xlab = "Habitat Coefficient (scaled)", ylab = "Density")
    abline(v = x$habitat_analysis$optima_coefficient_mean, col = "blue", lwd = 2)
    abline(v = x$habitat_analysis$optima_coefficient_ci_lower, col = "blue", lty = 2)
    abline(v = x$habitat_analysis$optima_coefficient_ci_upper, col = "blue", lty = 2)
    abline(v = 0, col = "gray", lty = 2)
  }
  
  # Plot based on selection
  if (which == "density") {
    plot_density()
  } else if (which == "scatter") {
    plot_scatter()
  } else if (which == "habitat") {
    plot_habitat()
  } else if (which == "both") {
    plot_density()
    plot_scatter()
  } else if (which == "all") {
    plot_density()
    plot_scatter()
    plot_habitat()
    # Fourth panel - could add population-level effects if desired
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Reserved for additional plot", cex = 1.2)
  }
  
  invisible(x)
}