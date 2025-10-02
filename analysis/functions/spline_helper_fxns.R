
#..................................................
# Script:   spline_helper_fxns.R
#
# Purpose:
#   Provide helper functions for spline-based logistic regression models
#   to compute predicted log-odds, odds ratios, and confidence intervals.
#
# Functions:
#   - compute_OR_CIs(): computes log-odds ratios, ORs, 95% CIs, and significance
#     for both unstratified and stratified reference groups
#   - fit_gam_te_model(): fits GAM model with tensor product smoothing term
#   - predict_gam(): Generate predictions from a fitted GAM model over specified 
#                     sequence of wind and distance values
#..................................................

# -----------------------------------------------------------------------------#
# FUNCTION: compute_OR_CIs
# PURPOSE:
#   Compute log-odds ratios, odds ratios (OR), 95% confidence intervals (CI),
#   and significance for predictions from a fitted logistic regression model
#   or a GAM (mgcv::gam) model. Handles both unstratified and stratified 
#   reference groups. For stratified references, the mean design vector is 
#   computed separately for each group.
#
# INPUT:
#   fitted_model   : glm or gam object; fitted model
#   unexposed_data : data frame of reference/unexposed observations
#                    - For GLM with spline bases: include manually created spline 
#                      basis columns + covariates at reference levels.
#                    - For GAM: raw observed unexposed rows + covariates at reference
#                       levels ; spline/basis handled internally by the model.
#   pred_data      : data frame of prediction points
#                    - For GLM: must include manually created spline basis + covariates
#                    - For GAM: raw data + covar at reference levels; 
#                               model handles spline internally
#                    - Must include a column named `log_odds` with model predictions
#   log_odds_ref   : numeric vector; predicted log-odds for reference/unexposed group
#   strat_var      : character or NULL; optional stratifying variable name
#
# OUTPUT:
#   pred_data: original prediction dataset augmented with:
#              - log_odds_ratio
#              - se_ratio
#              - odds_ratio
#              - or_low, or_high (95% CI)
#              - significant (logical)
# -----------------------------------------------------------------------------#
compute_OR_CIs <- function(fitted_model, 
                           unexposed_data, 
                           pred_data,
                           log_odds_ref,
                           strat_var){
  # Determine model type (gam vs. glm)
  is_gam <- inherits(fitted_model, "gam")
  
  # Design matrices
  if (is_gam) {
    # GAM: use lpmatrix to capture the evaluated smooths and parametric terms
    # at the new data points. The lpmatrix is a design matrix of the model's 
    # basis functions (including any tensor product or spline terms) evaluated 
    # at new_data
    X_pred <- predict(fitted_model, newdata = pred_data, type = "lpmatrix")
    X_ref  <- predict(fitted_model, newdata = unexposed_data, type = "lpmatrix")
  } else {
    # GLM: use model.matrix to construct the design matrix directly from the
    # predictor variables. This captures both parametric terms and any spline
    # basis terms manually created for GLM models. 
    X_pred <- model.matrix(delete.response(terms(fitted_model)), pred_data)
    X_ref  <- model.matrix(delete.response(terms(fitted_model)), unexposed_data)
  }
  
  # Coefficient covariance
  vcov_beta <- vcov(fitted_model)
  
  
  if(is.null(strat_var)){     # unstratified 
    # model matrix, reference group
    # X_ref  <- model.matrix(delete.response(terms(fitted_model)), unexposed_data) 
    # # design matrix for prediction grid
    # X_pred <- model.matrix(delete.response(terms(fitted_model)), pred_data) 

    # Average design vector, reference group
    x_ref_mean  <- colMeans(X_ref) 
    
    # coefficient covariance 
    #vcov_beta <- vcov(fitted_model) 
    
    # For each prediction distance, compute contrast = x(z) - x_ref_mean
    contrasts <- sweep(X_pred, 2, x_ref_mean, FUN = "-")  # subtract x_ref_mean from each row
    
    # Compute variance of the log-odds difference (delta)
    var_diff <- rowSums((contrasts %*% vcov_beta) * contrasts)  
    se_diff <- sqrt(var_diff)
    
    pred_data <- pred_data %>%
      mutate(
        log_odds_ratio = log_odds - log_odds_ref, 
        se_ratio       = se_diff,  
        odds_ratio     = exp(log_odds_ratio),
        or_low         = exp(log_odds_ratio - 1.96 * se_ratio),
        or_high        = exp(log_odds_ratio + 1.96 * se_ratio),
        significant    = or_low > 1 | or_high < 1
      )
  }
  if(!is.null(strat_var)){ # stratified
    # stratified case
    strat_levels <- levels(pred_data[[strat_var]])
    
    # Split reference matrix by stratifying levels
    # X_ref  <- model.matrix(delete.response(terms(fitted_model)), unexposed_data)
    # X_pred <- model.matrix(delete.response(terms(fitted_model)), pred_data)
    
    # Average design vector for reference group: Compute for each level 
    # of stratifying variable
    # (Notes: computes the mean design vector separately for each level, 
    # so each level curve is compared to the observed unexposed mean within that level)
    ref_split <- split(as.data.frame(X_ref), unexposed_data[[strat_var]])
    
    ref_means_list <- lapply(strat_levels, function(lvl) {
      if (!is.null(ref_split[[lvl]])) colMeans(as.matrix(ref_split[[lvl]])) else NA
    })
    names(ref_means_list) <- strat_levels
    ref_means_mat <- do.call(rbind, ref_means_list)
    
    # Match each prediction row to its stratifying reference mean
    ref_for_pred <- ref_means_mat[ as.character(pred_data[[strat_var]]), , drop = FALSE ]
    
    # Contrasts and variance
    #vcov_beta <- vcov(fitted_model)
    contrasts <- X_pred - ref_for_pred
    var_diff  <- rowSums((contrasts %*% vcov_beta) * contrasts)
    se_diff   <- sqrt(pmax(0, var_diff))
    
    # Compute mean log-odds for each stratified reference group
    mean_log_odds_by_group <- tapply(log_odds_ref, unexposed_data[[strat_var]], mean)
    mean_log_odds_for_pred <- mean_log_odds_by_group[ as.character(pred_data[[strat_var]]) ]
    
    # Assemble predictions
    pred_data <- pred_data %>%
      mutate(
        log_odds_ratio = log_odds - mean_log_odds_for_pred,
        se_ratio       = se_diff,
        odds_ratio     = exp(log_odds_ratio),
        or_low         = exp(log_odds_ratio - 1.96 * se_ratio),
        or_high        = exp(log_odds_ratio + 1.96 * se_ratio),
        significant    = (or_low > 1) | (or_high < 1)
      )
  }
  return(pred_data)
}


# -----------------------------------------------------------------------------#
# FUNCTION: fit_gam_te_model
# PURPOSE:
#   Fit a binomial Generalized Additive Model (GAM) with a tensor product smooth
#   term for distance and wind variables. Includes specified covariates and 
#   allows dynamic specification of knots and basis types. Computes model 
#   diagnostics including k-check and concurvity.
#
# INPUT:
#   dist_var       : character; name of the distance variable
#   wind_var       : character; name of the wind variable
#   covariates_eq  : character; formula string of covariates (e.g., "age + sex")
#   k_vec          : numeric vector of length 2; number of knots for each term in te()
#   bs_vec         : character vector of length 2; basis types for each term in te()
#   data           : data frame; dataset containing all variables
#
# OUTPUT:
#   output: list containing:
#           - model      : fitted gam object
#           - formula    : formula string used for fitting
#           - fit_info   : results of k.check diagnostics
#           - concurvity : concurvity diagnostics
# -----------------------------------------------------------------------------#
fit_gam_te_model <- function(
    dist_var,                 
    wind_var,                 
    covariates_eq,            
    k_vec = c(7, 5),          
    bs_vec = c("tp", "cr"),   
    data                      
) {
  # Dynamically construct formula
  te_term <- sprintf("te(%s, %s, k = c(%d, %d), bs = c('%s', '%s'))",
                     dist_var, wind_var,
                     k_vec[1], k_vec[2],
                     bs_vec[1], bs_vec[2])
  
  full_formula_str <- paste("case ~", te_term, "+", covariates_eq)
  full_formula     <- as.formula(full_formula_str)
  print(full_formula_str)

  # Fit model
  mod <- gam(
    formula = full_formula,
    family  = binomial,
    data    = data,
    method  = "REML"
  )
  print("Done fitting gam.")

  # Run model diagnostics
  fit_info    <- k.check(mod)
  summary_mod <- summary(mod)
  concurv     <- mgcv::concurvity(mod, full = TRUE)
  # smooth_info <- data.frame(
  #   k_prime = fit_info[,"k'"],
  #   edf = fit_info[,"edf"]
  #   k_index = fit_info[, "k-index"],
  #   k_p_value = fit_info[, "p-value"],
  #   te_p_value = summary_mod$s.table[, "p-value"],
  #   row.names = NULL
  # )
  print("Done fitting diagnostics.")

  # Save output
  output <- list(
    model = mod,
    formula = full_formula_str,
    fit_info = fit_info,
    concurvity = concurv
  )
  return(output)
}


# -----------------------------------------------------------------------------#
# FUNCTION: predict_gam
# PURPOSE:
#   Generate predictions from a fitted GAM with a tensor product smooth, 
#   compute log-odds, standard errors, odds ratios (OR), 95% confidence intervals,
#   and significance relative to an unexposed reference group. Handles setting 
#   covariates to reference levels and defining an unexposed group dynamically.
#
# INPUT:
#   mod                     : gam object; fitted GAM model
#   data                    : data frame; dataset used for GAM model fitting
#   dist_var                : character; name of the distance variable in model
#   wind_var                : character; name of the wind variable in model
#   covariates_eq           : character; formula string of covariates (e.g., "age + sex")
#   ref_covariates          : named list; reference levels for covariates
#   distance_seq            : numeric vector; distances (km) at which to predict
#   wind_seq                : numeric vector; wind values at which to predict
#   unexposed_dist_threshold: numeric; threshold distance to define unexposed group
#   unexposed_wind_value    : numeric; value of wind variable to define unexposed group
#
# OUTPUT:
#   grid_data: data frame containing prediction grid augmented with:
#              - log_odds
#              - se
#              - log_odds_ratio_group
#              - se_ratio
#              - odds_ratio_group
#              - or_low_group, or_high_group (95% CI)
#              - significant_group (logical)
# Dependencies: 
#   compute_OR_CIs() function
# -----------------------------------------------------------------------------#
predict_gam <- function(mod,
                        data,
                        dist_var,
                        wind_var,
                        covariates_eq,
                        ref_covariates,
                        distance_seq   = seq(0, 30, length.out = 100),
                        wind_seq       = seq(0, 1, length.out = 100),
                        unexposed_dist_threshold = 25, 
                        unexposed_wind_value     = 0) {

  # Create prediction grid
  grid_data    <- expand.grid(distance = distance_seq, wind = wind_seq)
  names(grid_data)[names(grid_data) == "distance"] <- dist_var
  names(grid_data)[names(grid_data) == "wind"]     <- wind_var
  
  # Add covariates at specified reference levels
  covariates <- strsplit(covariates_eq, "\\s*\\+\\s*")[[1]]
  covariates <- trimws(covariates)
  for (var in covariates) {
    grid_data[[var]] <- ref_covariates[[var]]
  }
  
  # Predict log-odds over prediction grid
  preds              <- predict(mod, newdata = grid_data, type = "link", se.fit = TRUE)
  grid_data$log_odds <- preds$fit
  
  # Filter to unexposed population (default: > 25km and never downwind)
  unexposed_data <- data %>%
    filter(!!as.name(dist_var) > unexposed_dist_threshold,
           !!as.name(wind_var) == unexposed_wind_value)
  # Add covariates at specified reference levels
  for (var in covariates) {
    unexposed_data[[var]] <- ref_covariates[[var]]
  }
  
  # Predict baseline log-odds for unexposed group
  ref_preds <- predict(mod, newdata = unexposed_data, type = "link", se.fit = TRUE)
  mean_log_odds_baseline <- mean(ref_preds$fit)

  # Compute ORs and confidence intervals
  grid_data <- compute_OR_CIs(
    fitted_model   = mod,
    unexposed_data = unexposed_data,
    pred_data      = grid_data,
    log_odds_ref   = mean_log_odds_baseline,
    strat_var      = NULL  
  )
  return(grid_data)
}


# STOP

