
#..................................................
# Script:   EpiBuffer_resultssummary_fxns.R
#
# Purpose:
#   Provide helper functions to summarize outputs from EpiBuffer models,
#   including posterior samples and summary statistics (HDI, CI, mean, median).
#
# Functions:
#   - calculate_mode(): compute the statistical mode of a vector
#   - compute_max_exp(): compute maximum possible exposure given maximum radius
#   - compute_exposure(): compute exposure for a vector of radii across 
#                         different exposure definitions (cumulative, spherical, binary)
#   - recover_true_eta(): rescale fitted eta parameters back to original scale
#   - recover_true_theta(): rescale fitted theta parameters back to original scale
#   - get_posterior_samples_radius(): apply burn-in and thinning to posterior radius samples
#   - get_posterior_samples_theta(): unscale theta, apply burn-in/thinning
#   - get_posterior_samples_othervar(): summarize posterior samples for beta, eta, gamma, rho_phi
#   - compute_posterior_summary(): summarize posterior draws with mean, median, HDI, CI
#   - get_summary_by_cluster(): compute posterior summaries at the cluster level
# Dependencies:
#   R packages: dplyr, tidyr, HDInterval
#..................................................

library(dplyr)
library(tidyr)
library(HDInterval)
# -----------------------------------------------------------------------------#
# FUNCTION: calculate_mode
# PURPOSE: Compute the mode of some data
# INPUT:
#   x: a vector
# OUTPUT: Numeric vector, mode.
# -----------------------------------------------------------------------------#
calculate_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# FUNCTION: compute_max_exp
# PURPOSE:
#   Compute the maximum exposure given the maximum radius. 
#   Used as a scaling factor in EpiBuffer models.
# INPUT: 
#   b_radius        : Integer, maximum radius
#   exposure_dists  : (n_ind_unique x m) matrix, pairwise distances
#   v               : (n_ind x n_ind_unique) matrix, location indicators
# OUTPUT:
#   Numeric scalar, maximum possible exposure
# -----------------------------------------------------------------------------#

compute_max_exp <- function(b_radius, 
                            exposure_dists, 
                            v){
  # Validate inputs
  if (!is.numeric(b_radius) || length(b_radius) != 1) {
    stop("b_radius must be a single numeric value.")
  }
  if (!is.matrix(exposure_dists) || !is.matrix(v)) {
    stop("exposure_dists and v must be matrices.")
  }
  
  # v%*%exposure_dists = (n_ind x n_ind_unique)(n_ind_unique x m) = (n_ind x m)
  v_exposure_dists <- (v%*%exposure_dists) # (n_ind x m)
  n_ind <- nrow(v_exposure_dists)
  m <- ncol(exposure_dists) # number of exposure sources
  
  # maximum possible exposure - i.e. max count of facilities in max radius
  max_radius_mat <- matrix(b_radius, nrow = n_ind, ncol = m)
  max_exp_vec <- rowSums(v_exposure_dists < max_radius_mat)
  max_exp <- max(max_exp_vec)
  return(max_exp)
}
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# FUNCTION: compute_exposure
# PURPOSE:
#   Compute exposure for each individual given a vector of radii.
# INPUT:
#   exposure_definition_indicator : Integer (0=cumulative, 1=spherical, 2=binary)
#   div_m_max                     : Logical, divide by maximum exposure if TRUE
#   m_max                         : Numeric, max possible exposure
#   radius_max                    : Numeric, max radius
#   exposure_dists                : (n_ind_unique x m) matrix, pairwise distances
#   v                             : (n_ind x n_ind_unique) matrix, location indicators
#   radius_vec                    : Numeric vector, radii per individual
# OUTPUT:
#   Numeric vector of exposures
# -----------------------------------------------------------------------------#
compute_exposure <- function(exposure_definition_indicator, 
                             div_m_max = FALSE,
                             m_max = NULL,
                             radius_max = NULL, 
                             exposure_dists, 
                             v, 
                             radius_vec
){
  # v_exposure_dists: distance between unique individuals and exposure points
  v_exposure_dists <- (v%*%exposure_dists) # (n_ind x m)
  n_ind <- nrow(v_exposure_dists) # number of individuals
  
  # compute exposure
  if(exposure_definition_indicator == 0){ # "cumulative_count"
    exposure <- rowSums(v_exposure_dists < radius_vec) 
  }
  if(exposure_definition_indicator == 1){ # "spherical"
    # radius is the distance over which correlation is nonzero
    spherical_correlation <- (1 - 
                                1.5 * (v_exposure_dists / radius_vec) + 
                                0.5 * (v_exposure_dists / radius_vec)^3)
    
    # Apply mask: Set spatial correlation to 0 where distance exceeds the radius
    mask <- v_exposure_dists < radius_vec
    spherical_correlation[!mask] <- 0 
    exposure <- rowSums(spherical_correlation) #/ max_exp
  }
  if(exposure_definition_indicator == 2){ # "binary"
    exposure <- apply(v_exposure_dists < radius_vec, 1, max)
  }
  if(div_m_max == TRUE){
    exposure <- exposure / m_max
  }
  return(exposure)
}
# -----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#
# FUNCTION: recover_true_eta
# PURPOSE:
#   Rescale fitted eta parameters back to original scale.
# INPUT:
#   etas_fitted  : matrix of fitted eta values
#   m_max        : Numeric, max exposure
#   radius_range : Numeric vector, [min, max] radii
# OUTPUT:
#   Matrix of unscaled eta parameters
# -----------------------------------------------------------------------------#
recover_true_eta <- function(etas_fitted, 
                             m_max, 
                             radius_range){
  a <- radius_range[1]
  b <- radius_range[2]
  p_d <- nrow(etas_fitted) - 1
  
  if(p_d > 2){ 
    stop("Function does not support p_d > 2. Please ensure p_d <= 2.")
  }
  etas_true <- matrix(NA, 
                      nrow = nrow(etas_fitted), # p_d
                      ncol = ncol(etas_fitted)) # mcmc_samples
  if(p_d == 0){
    
    etas_true[1,] <- (1/m_max)*etas_fitted[1,]
    
  }else if(p_d == 1){
    
    etas_true[1,] <- (1/m_max)*(etas_fitted[1,] -
                                  (a/(b - a))*etas_fitted[2,])
    
    etas_true[2,] <- (1/m_max)*((1/(b - a))*etas_fitted[2,])
    
  }else if(p_d == 2){
    
    etas_true[1,] <- (1/m_max)*(etas_fitted[1,] -
                                  (a/(b - a))*etas_fitted[2,] + 
                                  (a^2/(b - a)^2)*etas_fitted[3,]
    )
    
    etas_true[2,] <- (1/m_max)*((1/(b - a))*etas_fitted[2,] - 
                                  (2*a)/(b-a)^2*etas_fitted[3,]
    )
    
    etas_true[3,] <- (1/m_max)*(1/(b-a)^2*etas_fitted[3,])
    
  }
  return(as.matrix(etas_true))
}

# -----------------------------------------------------------------------------#
# FUNCTION: recover_true_theta
# PURPOSE:
#   Rescale fitted theta parameters back to original scale.
# INPUT:
#   thetas_fitted : matrix or vector of fitted theta values
#   m_max         : Numeric, max exposure
# OUTPUT:
#   Data frame of unscaled theta values
# -----------------------------------------------------------------------------#
recover_true_theta <- function(thetas_fitted, 
                               m_max){
  thetas_true <- thetas_fitted / m_max
  return(as.data.frame(thetas_true))
}

# -----------------------------------------------------------------------------#
# FUNCTION: get_posterior_samples_radius
# PURPOSE:
#   Apply burn-in and thinning to posterior radius samples.
# INPUT:
#   results_radius : Matrix of posterior radius samples
#   burn_in        : Numeric, iterations to discard
#   thin           : Numeric, thinning factor
# OUTPUT:
#   Data frame of thinned posterior radius samples
# -----------------------------------------------------------------------------#
get_posterior_samples_radius <- function(results_radius, 
                                         burn_in, 
                                         thin){
  # Validate inputs
  if (!is.matrix(results_radius)) {
    stop("results_radius must be a matrix with MCMC iterations as rows and individuals as columns.")
  }
  if (!is.numeric(burn_in) || burn_in < 0) {
    stop("burn_in must be a non-negative numeric value.")
  }
  if (!is.numeric(thin) || thin <= 0) {
    stop("thin must be a positive numeric value.")
  }
  
  # Reorganize radius df: Rows = MCMC iterations; columns = outcome units.
  if(ncol(results_radius) > 1){ #SpatialBuffer
    radius_df <- data.frame(t(results_radius)) 
  }else if(ncol(results_radius) == 1){ #SingleBuffer
    radius_df <- data.frame(results_radius)
    colnames(radius_df) <- c("radius")
  }
  
  radius_df$mcmc_iteration <- seq(1, nrow(radius_df))
  mcmc_samples <- max(radius_df$mcmc_iteration)
  keep_set <- seq((burn_in + 1), mcmc_samples, thin)
  radius_df <- radius_df[keep_set,] 
  
  return(radius_df)
}

# -----------------------------------------------------------------------------#
# FUNCTION: get_posterior_samples_theta
# PURPOSE:
#   Unscale theta samples and apply burn-in/thinning.
# INPUT:
#   results_theta : Matrix of posterior theta samples
#   m_max         : Numeric, max exposure
#   burn_in       : Numeric, iterations to discard
#   thin          : Numeric, thinning factor
# OUTPUT:
#   Data frame of thinned, unscaled posterior theta samples
# -----------------------------------------------------------------------------#
get_posterior_samples_theta <- function(results_theta, 
                                        m_max,
                                        burn_in, 
                                        thin){
  # Reorganize radius df: Rows = MCMC iterations; columns = outcome units.
  if(nrow(results_theta) > 1){ #SpatialBuffer
    theta_df <- data.frame(t(results_theta)) 
  }else if(nrow(results_theta) == 1){ #SingleBuffer
    theta_df <- data.frame(t(results_theta))
    colnames(theta_df) <- c("theta")
  }
  
  # unscale thetas
  theta_df <- recover_true_theta(theta_df, m_max = m_max)
  
  theta_df$mcmc_iteration <- seq(1, nrow(theta_df))
  mcmc_samples <- max(theta_df$mcmc_iteration)
  
  keep_set <- seq((burn_in + 1), mcmc_samples, thin)
  theta_df <- theta_df[keep_set, ]
  return(theta_df)
}

# -----------------------------------------------------------------------------#
# FUNCTION: get_posterior_samples_othervar ----
# PURPOSE:
#   Summarize posterior samples of beta, eta, gamma, rho_phi from EpiBuffer.
#   Apply burn-in/thinning and rescale etas/thetas.
# INPUT:
#   results           : List of posterior samples
#   model_input       : List, exposure distances and v matrix
#   exposure_indicator: Integer, exposure type
#   burn_in           : Numeric, iterations to discard
#   thin              : Numeric, thinning factor
#   radius_range      : Numeric vector, [min, max] radii
# DEPENDENCIES: 
#   calculate_mode()
# OUTPUT:
#   Dataframe of posterior samples after burnin and thinning removed
# -----------------------------------------------------------------------------#
# compute posterior summaries for beta, eta, gamma, rho_phi
get_posterior_samples_othervar <- function(results, 
                                           model_input,
                                           exposure_indicator, 
                                           burn_in, 
                                           thin, 
                                           radius_range){
  beta_df <- data.frame(t(results$beta))
  colnames(beta_df) <- paste0("beta", (seq(1:ncol(beta_df))-1))
  
  # return the etas both as they come from SpatialBuffer and unscaled
  eta_unscaled <- recover_true_eta(etas_fitted = results$eta, 
                                   m_max = results$exposure_scale, 
                                   radius_range = radius_range)
  
  eta_unscaled_df <- data.frame(t(eta_unscaled))
  colnames(eta_unscaled_df) <- paste0("eta.unscaled", (seq(1:ncol(eta_unscaled_df))-1))
  
  eta_df <- data.frame(t(results$eta))
  colnames(eta_df) <- paste0("eta.fitted", (seq(1:ncol(eta_df))-1))
  
  # exp*theta
  # Reorganize radius df: Rows = MCMC iterations; columns = outcome units.
  if(ncol(results$radius) > 1){ #SpatialBuffer
    radius_df <- data.frame(t(results$radius)) 
  }else if(ncol(results$radius) == 1){ #SingleBuffer
    radius_df <- data.frame(results$radius)
    colnames(radius_df) <- c("radius")
  }
  
  radius_df$mcmc_iteration <- seq(1, nrow(radius_df))
  mcmc_samples <- max(radius_df$mcmc_iteration)
  keep_set <- seq((burn_in + 1), mcmc_samples, thin)
  radius_df <- radius_df[keep_set,] 
  
  # Reorganize radius df: Rows = MCMC iterations; columns = outcome units.
  if(nrow(results$theta) > 1){ #SpatialBuffer
    theta_df <- data.frame(t(results$theta)) 
  }else if(nrow(results$theta) == 1){ #SingleBuffer
    theta_df <- data.frame(t(results$theta))
    colnames(theta_df) <- c("theta")
  }
  
  # unscale thetas
  theta_df <- recover_true_theta(theta_df, m_max = results$exposure_scale) %>% data.frame()
  theta_df <- theta_df[keep_set, , drop = FALSE] %>% as.matrix()
  
  # compute exposure for each individual (on each mcmc iteration)
  radius_df <- radius_df %>% dplyr::select(-mcmc_iteration)
  exposure <- apply(radius_df[1:nrow(radius_df), , drop = FALSE], 1, function(radius_row) {
    compute_exposure(
      exposure_definition_indicator = exposure_indicator, 
      div_m_max = FALSE,
      exposure_dists = model_input$exposure_dists, 
      v = model_input$v, 
      radius_vec = radius_row
    )
  }) # rows = individuals, columns = mcmc iterations
  exposure <- t(exposure) # rows = mcmc iterations, columns = individuals
  
  if(ncol(theta_df) == 1){ # 1 value per mcmc iteration
    # multpily theta vector with each row 
    exposure_times_theta <- exposure * matrix(theta_df, 
                                              nrow = nrow(exposure), 
                                              ncol = ncol(exposure), 
                                              byrow = FALSE)
  }else if(ncol(theta_df) > 1){
    exposure_times_theta <- exposure * theta_df
  }
  
  
  if(!is.null(results$gamma)){
    gamma_df <- data.frame(t(results$gamma))
    colnames(gamma_df) <- paste0("gamma", (seq(1:ncol(gamma_df))-1))
  }
  
  if(!is.null(results$rho_phi)){
    rho_phi_df <- data.frame(results$rho_phi)
    colnames(rho_phi_df) <- paste0("rhophi", (seq(1:ncol(rho_phi_df))-1))
  }
  
  if(!is.null(results$gamma)){
    allvar_df <- cbind(beta_df, eta_df, eta_unscaled_df, gamma_df, rho_phi_df)
  }else{
    allvar_df <- cbind(beta_df, eta_df, eta_unscaled_df)
  }
  
  allvar_df$mcmc_iteration <- seq(1, nrow(allvar_df))
  allvar_df <- allvar_df[keep_set, ] 
  
  return(list(allvar_thinned_samples = allvar_df, 
              exposure_thinned_samples = exposure,
              exp_theta_thinned_samples = exposure_times_theta))
}

# -----------------------------------------------------------------------------#
# FUNCTION: compute_posterior_summary ----
# PURPOSE:
#   Summarizes the posterior distribution of input data
#   - Posterior mean, HDI lower and upper 95, CI lower and upper 95
#   - If exponentiate = TRUE: take exp() of posterior values before computing 
#                             mean and intervals
# INPUT:
#   thinned_samples: Matrix, posterior samples after burn in and thinning removed 
#                    Rows = MCMC iterations. Columns = MCMC outcome units.
#   exponentiate    : Logical, exponentiate values before summary
# OUTPUT:
#   Dataframe of poster.mean, hdi.lower95, hdi.upper95, ci.lower95, ci.upper95
# -----------------------------------------------------------------------------#
compute_posterior_summary <- function(thinned_samples, 
                                      exponentiate = FALSE){
  if(exponentiate == TRUE){
    summary_df <- thinned_samples %>% 
      dplyr::select(-mcmc_iteration) %>%
      mutate(across(everything(), exp)) %>% # exponentiate before taking means and intervals
      dplyr::summarise(across(everything(), list( # apply to each column (individual)
        posterior.median = ~median(.),
        posterior.mean = ~mean(.),
        hdi.lower95 = ~hdi(., credMass = 0.95)[["lower"]],
        hdi.upper95 = ~hdi(., credMass = 0.95)[["upper"]],
        ci.lower95 = ~quantile(., 0.025),
        ci.upper95 = ~quantile(., 0.975)
      ))) %>%
      pivot_longer(cols = everything(), names_to = c("cluster", ".value"), names_sep = "_")
  }else if(exponentiate == FALSE){
    summary_df <- thinned_samples %>%
      dplyr::select(-mcmc_iteration) %>%
      dplyr::summarise(across(everything(), list( # apply to each column (individual)
        posterior.median = ~median(.),
        posterior.mean = ~mean(.),
        hdi.lower95 = ~hdi(., credMass = 0.95)[["lower"]],
        hdi.upper95 = ~hdi(., credMass = 0.95)[["upper"]],
        ci.lower95 = ~quantile(., 0.025),
        ci.upper95 = ~quantile(., 0.975)
      ))) %>%
      pivot_longer(cols = everything(), names_to = c("cluster", ".value"), names_sep = "_")
  }
  return(summary_df)
  
}

# -----------------------------------------------------------------------------#
# FUNCTION: get_summary_by_cluster
# PURPOSE:
#   Compute posterior summaries for parameters at cluster level.
#   Includes summaries for radius, theta, exp*theta, beta, gamma (if present).
# INPUT:
#   model_name         : Character, model identifier
#   all_results_summary: List of posterior summaries
#   exposure_indicator : Integer, exposure type
#   full_data          : Data frame, input data with cluster IDs
#   exponentiate_theta : Logical, exponentiate theta before summarizing
# OUTPUT:
#   List with cluster-level summaries in wide and long format
# -----------------------------------------------------------------------------#
# Reshape data to long format and save one observation per cluster
get_summary_by_cluster <- function(model_name, 
                                   all_results_summary, 
                                   exposure_indicator, 
                                   full_data, 
                                   exponentiate_theta = FALSE) {
  print(model_name)
  n_ind <- nrow(full_data)
  
  if(model_name %in% c("fixed", "single")) {
    radius <- compute_posterior_summary(all_results_summary[[model_name]]$radius_thinned_samples, 
                                        exponentiate = FALSE) %>%
      slice(rep(1, n_ind)) 
    
    theta <- compute_posterior_summary(all_results_summary[[model_name]]$theta_thinned_samples, 
                                      exponentiate = exponentiate_theta) %>%
      slice(rep(1, n_ind)) 
    
  } else {
    radius <-  compute_posterior_summary(all_results_summary[[model_name]]$radius_thinned_samples, 
                                         exponentiate = FALSE)
    theta <- compute_posterior_summary(all_results_summary[[model_name]]$theta_thinned_samples, 
                                       exponentiate = exponentiate_theta) 
  }
  exptheta_all <- data.frame(all_results_summary[[model_name]]$exp_theta_thinned)
  exptheta_all$mcmc_iteration <- rownames(exptheta_all)
  exptheta <- compute_posterior_summary(exptheta_all)
  
  exptheta_post_summary <- exptheta %>%
    mutate(cluster = str_extract(cluster, "X\\D*(\\d{3})") %>% str_remove_all("X\\D*")) %>%
    mutate(cluster = as.numeric(cluster)) %>%
    distinct(cluster, .keep_all = TRUE) %>%
    rename_with(~ paste0(.x, ".exptheta"), -cluster)
  
  radius_post_summary_cluster <- radius %>%
    mutate(cluster = as.numeric(full_data$v001)) %>%
    distinct(across(everything()), .keep_all = TRUE) %>%
    rename_with(~ paste0(.x, ".radius"), -cluster)
  
  theta_post_summary_cluster <- theta %>%
    mutate(cluster = as.numeric(full_data$v001)) %>%
    distinct(across(everything()), .keep_all = TRUE) %>%
    rename_with(~ paste0(.x, ".theta"), -cluster) 
  
  beta_post_summary <- compute_posterior_summary(
    all_results_summary[[model_name]]$othervar_thinned_samples %>% dplyr::select(starts_with("beta"), mcmc_iteration), 
    exponentiate = TRUE) %>%
    dplyr::mutate(parameter = cluster, 
                  cluster = "global", 
                  cluster_id_seq = 0, 
                  exptheta.hdi.NULL = NA)
  
  if(model_name %in% c("spatial0", "spatial1")){
    gamma_post_summary <- compute_posterior_summary(
      all_results_summary[[model_name]]$othervar_thinned_samples %>% dplyr::select(starts_with("gamma"), mcmc_iteration), 
      exponentiate = FALSE)  %>%
      dplyr::mutate(parameter = cluster, 
                    cluster = "global", 
                    cluster_id_seq = 0, 
                    exptheta.hdi.NULL = NA) 
    
  }
  
  cluster_summary <- full_join(radius_post_summary_cluster, 
                               theta_post_summary_cluster, 
                               by = "cluster") %>%
    full_join(exptheta_post_summary, by = "cluster") %>%
    mutate(cluster_id_seq = row_number())
  
  # Add indicator of exptheta posterior interval containing NULL value (0)
  cluster_summary <- cluster_summary %>%
    mutate(
      exptheta.hdi.NULL = if_else(
        hdi.lower95.exptheta <= 0 & hdi.upper95.exptheta >= 0,
        1L, 0L
      )
    )
  
  cluster_summary_long <- cluster_summary %>%
    pivot_longer(cols = -c(cluster, cluster_id_seq, exptheta.hdi.NULL), 
                 names_to = c(".value", "parameter"),
                 names_pattern = "(posterior\\.mean|posterior\\.median|hdi.lower95|hdi.upper95|ci.lower95|ci.upper95)\\.(radius|theta|exptheta)")
  
  col_order <- colnames(cluster_summary_long)
  
  beta_post_summary <- beta_post_summary %>% 
    dplyr::select(all_of(col_order)) # same col order 
  
  if(model_name %in% c("fixed", "single")){ # no gamma parameter
    cluster_summary_long <- rbind(cluster_summary_long, beta_post_summary)
  }else if(model_name %in% c("spatial0", "spatial1")){
    gamma_post_summary <- gamma_post_summary %>% 
      dplyr::select(all_of(col_order)) # same col order 
    cluster_summary_long <- rbind(cluster_summary_long, beta_post_summary, gamma_post_summary)
    
  }
  
  return(list(cluster_summary_long = cluster_summary_long, 
              cluster_summary = cluster_summary))
}

# STOP




