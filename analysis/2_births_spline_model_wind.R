#..............................................................................
# Script:   2_births_spline_model_wind.R
#
# Purpose:  Analysis of PTB as a function of distance to nearest CAFO
#           *and* proportion of gestational time spent downwind from nearest CAFO
#           Models are fit using Generalized Additive Models (GAMs) with tensor 
#           product smoothers. Predictions and marginal effects are generated 
#           across fine grids of exposure values and discrete distances.
#
# Inputs:
#   - receiver_dat_name: cleaned births dataset with distance and wind exposure
#   - receiver_filter_name: subset of births to define the analytic population
#
# Outputs:
#   - results_MODEL_*.rds: fitted GAM models by trimester and exposure variable
#   - results_PREDICTED_*.rds: grid-based predictions with confidence intervals
#   - results_MARGINAL_*.rds: marginal effect predictions at discrete distances
#
# Helper Functions:
#   - fit_gam_te_model() : fits GAM with tensor product smoothers
#   - predict_gam()      : generates predictions for full grids or marginals
#
# Dependencies:
#   - R packages: here, dplyr, splines, broom, mgcv, gratia, furrr, purrr, glue
#   - Custom scripts:
#       - code/setup_header.R
#       - code/analysis/functions/spline_helper_fxns.R
#..............................................................................
##########
# Setup
##########
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
source(file.path(here(), "cal-cafo-births/code/analysis/functions/spline_helper_fxns.R"))

library(dplyr)
library(splines)
library(mgcv)
library(furrr)
library(purrr)
library(glue)

######################
# Specify Parameters
######################
model_name               <- "windGAM"
outcome_var_raw          <- "prem3s"
dist_exposure_var_raw    <- "nearest_cafo_km" 
wind_exposure_var_suffix <- "prop_downwind"
ptb_level                <- 2
by_trimester             <- TRUE

covariates_continuous <- c("mothage")
covariates_cat        <- c("edum")
strat_var             <- c("racem", "sex")

# spline parameters
max_df       <- 10
ref_distance <- 25 # kilometers
k_vec_param  <- c(15, 15)
bs_vec_param <- c("ts", "ts")

covariates            <- c(covariates_continuous, covariates_cat)
covariates_eq         <- paste0(covariates, collapse = "+")
if(ptb_level == 1){ptb_name = "early"}
if(ptb_level == 2){ptb_name = "moderate"}

# make exposure var variables
if(by_trimester == TRUE){
  wind_exposure_var <- paste0("tri_", wind_exposure_var_suffix)
  wind_tri_var <- glue("tri{1:3}_{wind_exposure_var_suffix}")
  if(ptb_level == 1){
    wind_vars <- wind_tri_var[1:2]
  }
  if(ptb_level == 2){
    wind_vars <- wind_tri_var
  }
}
if(by_trimester == FALSE){
  if(ptb_level == 1){
    wind_exposure_var <- paste0("tri12_", wind_exposure_var_suffix)
    wind_tri_var <- glue("tri12_{wind_exposure_var_suffix}")
    wind_vars <- wind_tri_var
  }
  if(ptb_level == 2){
    wind_exposure_var <- paste0("tri123_", wind_exposure_var_suffix)
    wind_tri_var <- glue("tri123_{wind_exposure_var_suffix}")
    wind_vars <- wind_tri_var
  }
}

#################################
# Load data
#################################
receiver_dat_name    <- "births_main_analysis_df_cafo20250212_wedgeangle90_maxdist25000_NEARESTCAFO_ALLRECEIVERS_WINDEXP.csv"
receiver_filter_name <- "births_main_analysis_df.csv"
save_results_folder  <- "births_nullparous_dropduplicates_scorege50"
data_optional_name   <- NULL  # NULL #"_NC"

receiver <- fread(file.path(server_cleandata_path, receiver_dat_name))
receiver_filter <- fread(file.path(server_cleandata_path, receiver_filter_name))
receiver <- receiver %>% dplyr::filter(birth_id %in% receiver_filter$birth_id)

################
# Prepare Data
################
# drop obs that are NA for outcome of interest or any covariates
receiver <- receiver %>% filter(!is.na(!!sym(outcome_var_raw)))
receiver <- receiver %>% drop_na(all_of(covariates))
receiver <- receiver %>% drop_na(all_of(strat_var))

# Factor variables, scale continuous variables
receiver <- receiver %>%
  mutate(case = factor(case, levels = c(0, 1)),
         across(all_of(covariates_cat), as.factor), 
         across(all_of(outcome_var_raw), as.factor), 
         across(all_of(covariates_continuous), ~ as.vector(scale(.x)))) 

# Ensure reference levels for categorical covariates
receiver <- receiver %>%
  mutate(
    racem = relevel(factor(racem), ref = levels(factor(racem))[1]),
    sex   = relevel(factor(sex),   ref = levels(factor(sex))[1]),
    edum  = relevel(factor(edum),  ref = levels(factor(edum))[1])
  )

# Subset: term births and ptb category of interest
receiver_subset <- receiver %>% filter(case == 0 | prem3s == ptb_level)

# Create separate subsets per trimester
receiver_subset1 <- receiver_subset %>% 
  filter(tri1_reached == TRUE) %>%
  dplyr::select(case, 
                mothage, 
                edum, 
                all_of(wind_tri_var), 
                nearest_cafo_km)
receiver_subset2 <- receiver_subset %>% 
  filter(tri2_reached == TRUE) %>%
  dplyr::select(case, 
                mothage, 
                edum, 
                all_of(wind_tri_var), 
                nearest_cafo_km)
receiver_subset3 <- receiver_subset %>% 
  filter(tri3_reached == TRUE) %>%
  dplyr::select(case, 
                mothage, 
                edum, 
                all_of(wind_tri_var), 
                nearest_cafo_km)

# Define a vector of wind variable names
if(by_trimester == TRUE){
  if (ptb_level == 1) {
    wind_data_lookup <- list(receiver_subset1, receiver_subset2) %>%
      set_names(wind_vars)
  }
  if (ptb_level == 2) {
    wind_data_lookup <- list(receiver_subset1, receiver_subset2, receiver_subset3) %>%
      set_names(wind_vars)
  }
}
if(by_trimester == FALSE){
  if (ptb_level == 1) {
    wind_data_lookup <- list(receiver_subset2) %>%
      set_names(wind_vars)
  }
  
  if (ptb_level == 2) {
    wind_data_lookup <- list(receiver_subset3) %>%
      set_names(wind_vars)
  }
}

############################################
# Fit GAM model on each wind outcome variable
############################################
# Apply function in parallel over wind variables
plan(multisession, workers = length(wind_vars)) 
model_results_list <- future_map(
  wind_vars,
  function(wind_var_i) {
    print(glue::glue("Checking data for: {wind_var_i}"))
    print(colnames(wind_data_lookup[[wind_var_i]]))
    
    fit_gam_te_model(
      dist_var      = dist_exposure_var_raw,
      wind_var      = wind_var_i,
      covariates_eq = covariates_eq,
      k_vec         = k_vec_param,
      bs_vec        = bs_vec_param,
      data          = wind_data_lookup[[wind_var_i]]
    )
  }, 
  .options = furrr::furrr_options(seed = TRUE)
)
plan(sequential)
names(model_results_list) <- wind_vars

save_name <- paste0("results_MODEL_", model_name,
                    "_PARAM_", 
                    "k_", paste0(k_vec_param, collapse = "_"), 
                    "_bs_", paste0(bs_vec_param, collapse = "_"), 
                    "_OUTCOME_", outcome_var_raw, "_", ptb_name,  
                    "_EXP_", dist_exposure_var_raw, "_", wind_exposure_var,
                    "_COVAR_", paste0(covariates, collapse = "_"), 
                    "_REFDIST_", ref_distance, "km", data_optional_name, ".rds")

saveRDS(model_results_list, 
        file.path(transfer_modelresults_path, save_results_folder, save_name), 
        compress = TRUE)

#####################################################
# Predict using fitted GAM model (full grid predictions)
#####################################################
mothage_ref <- mean(receiver_subset$mothage, na.rm = TRUE)
edum_ref    <- factor("2", levels = levels(receiver_subset$edum))

all_pred_res <- imap(
  model_results_list,
  ~ predict_gam(
    mod                      = .x$model,
    data                     = wind_data_lookup[[.y]], # .y is the name of the list element
    dist_var                 = dist_exposure_var_raw,
    wind_var                 = .y,                     # use the name as the wind_var
    covariates_eq            = covariates_eq,
    ref_covariates           = list(mothage = mothage_ref, edum = edum_ref),
    distance_seq             = seq(0, 30, length.out = 100), # full grid
    wind_seq                 = seq(0, 1, length.out = 100),  # full grid
    unexposed_dist_threshold = 25, 
    unexposed_wind_value     = 0
  )
)

save_name <- paste0("results_PREDICTED_", model_name,
                    "_PARAM_", 
                    "k_", paste0(k_vec_param, collapse = "_"), 
                    "_bs_", paste0(bs_vec_param, collapse = "_"), 
                    "_OUTCOME_", outcome_var_raw, "_", ptb_name,  
                    "_EXP_", dist_exposure_var_raw, "_", wind_exposure_var,
                    "_COVAR_", paste0(covariates, collapse = "_"), 
                    "_REFDIST_", ref_distance, "km", data_optional_name, ".rds")

saveRDS(all_pred_res, 
        file.path(transfer_modelresults_path, save_results_folder, save_name), 
        compress = TRUE)

#####################################################
# Marginal Effects - predict effect at discrete distances
#####################################################
all_marginal_res <- imap(
  model_results_list,
  ~ predict_gam(
    mod                      = .x$model,
    data                     = wind_data_lookup[[.y]],      # .y is the name of the list element
    dist_var                 = dist_exposure_var_raw,
    wind_var                 = .y,                          # use the name as the wind_var
    covariates_eq            = covariates_eq,
    ref_covariates           = list(mothage = mothage_ref, edum = edum_ref),
    distance_seq             = c(4, 8, 12, 16, 20, 24),     # discrete distances for prediction
    wind_seq                 = seq(0, 1, length.out = 100), # full grid of wind values
    unexposed_dist_threshold = 25, 
    unexposed_wind_value     = 0
  )
)

save_name <- paste0("results_MARGINAL_", model_name,
                    "_PARAM_", 
                    "k_", paste0(k_vec_param, collapse = "_"), 
                    "_bs_", paste0(bs_vec_param, collapse = "_"), 
                    "_OUTCOME_", outcome_var_raw, "_", ptb_name,  
                    "_EXP_", dist_exposure_var_raw, "_", wind_exposure_var,
                    "_COVAR_", paste0(covariates, collapse = "_"), 
                    "_REFDIST_", ref_distance, "km", data_optional_name, ".rds")

saveRDS(all_marginal_res, 
        file.path(transfer_modelresults_path, save_results_folder, save_name), 
        compress = TRUE)


# STOP



