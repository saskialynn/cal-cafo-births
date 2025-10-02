
#..............................................................................
# Script:   EpiBuffer_mcmcdiagnostics_fxns.R
#
# Purpose:
#   Provide helper functions to calculate MCMC diagnostics from outputs of 
#   EpiBuffer models (SingleBuffer). Diagnostics include 
#   effective sample size and Geweke convergence statistics.
#
# Functions:
#   - mcmc_diagnostics_singlebuffer(): Compute diagnostics for SingleBuffer output
#
# Dependencies:
#   - R package: coda
#..............................................................................

library(coda)

# -----------------------------------------------------------------------------#
# FUNCTION: mcmc_diagnostics_singlebuffer
# PURPOSE:
#   Compute MCMC diagnostics (effective sample size, Geweke diagnostic) for the 
#   output of EpiBuffer::SingleBuffer.
#
# INPUTS:
#   results : List, output from EpiBuffer::SingleBuffer 
#   burn_in : Integer, number of initial MCMC samples to discard
#   thin    : Integer, thinning factor for MCMC samples
#
# OUTPUT:
#   List containing:
#     - n_eff             : Effective sample size for each parameter
#     - geweke            : Geweke diagnostic z-scores for each parameter
#     - n_eff_min         : Minimum effective sample size across parameters
#     - n_eff_median      : Median effective sample size across parameters
#     - geweke_range      : Range of Geweke diagnostic z-scores
#     - geweke_count_large: Count of parameters with large Geweke z-scores
# -----------------------------------------------------------------------------#
mcmc_diagnostics_singlebuffer <- function(results, 
                                          burn_in, 
                                          thin){
  # Burnin and thinning
  keep_set <- seq(from = (burn_in + 1), to = ncol(results$beta), by = thin)
  
  all <- rbind(
    results$beta[ ,keep_set], 
    results$radius[keep_set], # radius: 1 x mcmc
    results$theta[keep_set])  # theta: 1 x mcmc
  
  all <- t(all)
  
  # Compute diagnostics
  n_eff <- rep(NA, times = ncol(all))
  geweke <- rep(NA, times = ncol(all))
  
  for(j in 1:ncol(all)){
    n_eff[j] <- coda::effectiveSize(all[,j])
    geweke[j] <- coda::geweke.diag(all[,j])[[1]]
  }
  
  n_eff_min <- min(n_eff, na.rm = TRUE)
  n_eff_median <- median(n_eff, na.rm = TRUE)
  geweke_range <- range(geweke, na.rm = TRUE)
  geweke_count_large <- sum(abs(geweke) > qnorm(1.00 - (0.05/ncol(all))/2))
  
  diagnostics <- list(n_eff = n_eff, 
                      geweke = geweke, 
                      n_eff_min = n_eff_min, 
                      n_eff_median = n_eff_median,
                      geweke_range = geweke_range, 
                      geweke_count_large = geweke_count_large)
  return(diagnostics)
  
}
# -----------------------------------------------------------------------------#


# STOP
