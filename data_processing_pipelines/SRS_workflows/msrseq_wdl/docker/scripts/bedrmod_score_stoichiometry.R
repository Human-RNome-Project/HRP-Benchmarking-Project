# ============================================
# bedrMod Score Calculation: Stoichiometry-Based Methods
# For: GLORI, BSSeq, and similar methods
# ============================================
#
# This script calculates bedrMod confidence scores for methods that
# provide modification stoichiometry (fraction of modified reads).
#
# Statistical approach: Binomial test against background conversion rate
# H0: true modification rate = background rate
# H1: true modification rate > background rate (one-sided)
#
# Score = round(-log10(p-value)), capped at 1000
#
# Input data requirements:
#   - n_mod: number of reads supporting modification (e.g., C->T conversions)
#   - n_total: total coverage at that position
#
# ============================================

calculate_bedrmod_score_stoichiometry <- function(data, 
                                                   background_rate = NULL,
                                                   min_coverage = 10,
                                                   auto_estimate_bg = TRUE) {
  # ------------------------------------------
  # Input validation
  # ------------------------------------------
  required_cols <- c("n_mod", "n_total")
  if (!all(required_cols %in% colnames(data))) {
    stop("Data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  if (any(data$n_mod < 0, na.rm = TRUE) || any(data$n_total < 0, na.rm = TRUE)) {
    stop("n_mod and n_total cannot be negative")
  }
  
  if (any(data$n_mod > data$n_total, na.rm = TRUE)) {
    stop("n_mod cannot exceed n_total")
  }
  
  # ------------------------------------------
  # Step 1: Filter for minimum coverage
  # ------------------------------------------
  data_orig <- data
  data <- data[data$n_total >= min_coverage, ]
  
  if (nrow(data) == 0) {
    stop("No positions pass min_coverage threshold of ", min_coverage)
  }
  
  if (nrow(data) < 100) {
    warning("Less than 100 positions after filtering - p-values may be unreliable")
  }
  
  message("Retained ", nrow(data), " of ", nrow(data_orig), 
          " positions after min_coverage filter")
  
  # ------------------------------------------
  # Step 2: Calculate observed modification rate
  # ------------------------------------------
  data$mod_rate <- data$n_mod / data$n_total
  
  # ------------------------------------------
  # Step 3: Determine background rate
  # ------------------------------------------
  if (is.null(background_rate)) {
    if (auto_estimate_bg) {
      # Estimate from lower quartile (assumes most sites unmodified)
      nonzero_rates <- data$mod_rate[data$mod_rate > 0]
      
      if (length(nonzero_rates) < 50) {
        stop("Too few non-zero sites to estimate background rate")
      }
      
      bg_rate <- quantile(nonzero_rates, 0.25)
      bg_rate <- max(bg_rate, 1e-6)
      message("Auto-estimated background rate: ", round(bg_rate, 4))
    } else {
      stop("background_rate must be provided or auto_estimate_bg must be TRUE")
    }
  } else {
    bg_rate <- background_rate
    message("Using provided background rate: ", round(bg_rate, 4))
  }
  
  # Validate background rate
  if (bg_rate <= 0 || bg_rate >= 1) {
    stop("Background rate must be between 0 and 1 (exclusive), got: ", bg_rate)
  }
  
  # Warn if background rate seems unusually high
  if (bg_rate > 0.1) {
    warning("Background rate > 10% - verify this is expected for your method")
  }
  
  # ------------------------------------------
  # Step 4: Binomial test for each position
  # ------------------------------------------
  message("Calculating p-values...")
  data$pvalue <- mapply(function(k, n) {
    if (n == 0) return(1)
    binom.test(k, n, p = bg_rate, alternative = "greater")$p.value
  }, data$n_mod, data$n_total)
  
  # ------------------------------------------
  # Step 5: Calculate bedrMod score
  # ------------------------------------------
  # score = round(-log10(pvalue)), capped at 1000
  # Handle very small p-values to avoid -Inf
  data$score <- pmin(round(-log10(pmax(data$pvalue, 1e-100))), 1000)
  data$score[is.na(data$score)] <- 0  # handle any remaining NA
  
  # ------------------------------------------
  # Summary statistics
  # ------------------------------------------
  message("Score calculation complete")
  message("Score range: ", min(data$score), " to ", max(data$score))
  message("Positions with score > 20 (p < 0.01): ", sum(data$score > 20))
  message("Positions with score > 30 (p < 0.001): ", sum(data$score > 30))
  
  # Store parameters as attributes for reference
  attr(data, "background_rate") <- bg_rate
  attr(data, "min_coverage") <- min_coverage
  attr(data, "method") <- "binomial_test"
  
  return(data)
}


# ============================================
# Example usage
# ============================================

# # Generate example data
# set.seed(42)
# example_data <- data.frame(
#   position = 1:1000,
#   n_mod = rpois(1000, 5),
#   n_total = rpois(1000, 100) + 10  # ensure decent coverage
# )
# 
# # With provided background rate
# result <- calculate_bedrmod_score_stoichiometry(
#   example_data, 
#   background_rate = 0.01,
#   min_coverage = 10
# )
# 
# # With auto-estimated background rate
# result_auto <- calculate_bedrmod_score_stoichiometry(
#   example_data, 
#   auto_estimate_bg = TRUE,
#   min_coverage = 10
# )
# 
# # View results
# head(result)
# hist(result$score, breaks = 50, main = "bedrMod Score Distribution")
