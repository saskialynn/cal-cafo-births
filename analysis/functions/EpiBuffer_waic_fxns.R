
#..................................................
# Script:   EpiBuffer_waic_fxns.R
# Purpose:
#   Provide helper functions to calculate WAIC from outputs of EpiBuffer models.
#
# Functions:
#   - compute_waic(): calculate WAIC and effective number of parameters (pWAIC) 
#                     from posterior log-likelihood draws.
#..................................................

# -----------------------------------------------------------------------------#
# FUNCTION: compute_waic
# PURPOSE:
#   Compute WAIC and the effective number of parameters (pWAIC) based on 
#   posterior log-likelihood draws from EpiBuffer (uses formulas from Gelman (2013)).
# INPUT:
#   results : List, EpiBuffer model output containing `log_density` and `r`
#   burn_in : Numeric, number of burn-in iterations to discard
#   thin    : Numeric, thinning factor for MCMC samples
# OUTPUT:
#   List with:
#     - waic  : WAIC value
#     - pwaic : effective number of parameters
# -----------------------------------------------------------------------------#
compute_waic <- function(results, 
                         burn_in, 
                         thin){
  mcmc_samples <- length(results$r)
  
  keep_set <- seq((burn_in + 1), mcmc_samples, thin)
  
  piece <- results$log_density[,keep_set]
  
  llpd <- sum(log(rowMeans(exp(piece)))) # Gelman (2013), (eq 5)
  
  PWAIC_2 <- sum(apply(piece, 1, var))   # Gelman 2013, (eq 11)
  
  WAIC_2 <- -2*(llpd - PWAIC_2)
  return(list(waic = WAIC_2, pwaic = PWAIC_2))
}



