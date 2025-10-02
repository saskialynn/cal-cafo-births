
#...............................................................................
# Script: 3_Bayesian_SingleBR_analysis.R
#
# Purpose:
# Fit a Bayesian  model to estimate the distance range and magnitude of effects of
# ``exposure" to CAFOs and PTB; compute fit/MCMC diagnostics.
#
# Data Inputs:
# - Cleaned births dataset (case-control format) with covariates and outcomes
# - CAFO facilities data (CSV and spatial RDS formats)
#
# Outputs:
# - RDS files containing:
#     * Full MCMC results from SingleBuffer
# - WAIC values and MCMC diagnostic summaries
# - Log file capturing script output
#
# Helper Functions:
# - EpiBuffer::SingleBuffer(): fits the SingleBuffer model
# - mcmc_diagnostics_singlebuffer(): computes MCMC convergence diagnostics
# - compute_waic(): computes WAIC for model comparison
# - get_posterior_samples_*(): extracts thinned posterior samples
#
# Dependencies:
# - R packages: dplyr, sf, units, HDInterval, broom
# - Specialty packages: EpiBuffer (github.com/warrenjl/EpiBuffer )
# - Custom scripts: 
#     * code/setup_header.R
#     * EpiBuffer_waic_fxns.R
#     * EpiBuffer_mcmcdiagnostics_fxns.R
#     * EpiBuffer_resultssummary_fxns.R
#
#...............................................................................

library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
source(file.path(analysis_path, "functions/EpiBuffer_waic_fxns.R"))
source(file.path(analysis_path, "functions/EpiBuffer_mcmcdiagnostics_fxns.R"))
source(file.path(analysis_path, "functions/EpiBuffer_resultssummary_fxns.R"))

library(EpiBuffer)
library(sf)
library(units)
library(HDInterval)
library(broom)

#################################
# Specify Parameters
##################################
# General data parameters
model_name       <- "SingleBR"
outcome_var_raw  <- "prem3s"
exposure_var_raw <- "nearest_cafo_km" 
ptb_level        <- 1 # 1 = very, 2 = moderate
case_ratio       <- 1
control_ratio    <- 5
seed             <- 1234
optional_name    <- NULL

covariates_continuous <- c("mothage")
covariates_cat        <- c("edum")
strat_var             <- c("racem", "sex")
covariates    <- c(covariates_continuous, covariates_cat)
covariates_eq <- paste0(covariates, collapse = "+")

if(ptb_level == 1){ptb_name = "early"}
if(ptb_level == 2){ptb_name = "moderate"}

# SingleBuffer model parameters
likelihood_indicator          <- 0 # Binomial
exposure_definition_indicator <- 2 # presence/absence
radius_range                  <- c(0, 40) # lower and upper bound, km

# MCMC Parameters
mcmc_samples      <- 40000
burnin            <- 20000
thin              <- 2
metrop_var_radius <- 1.0

# Data file names / paths
save_results_folder <- "births_nullparous_dropduplicates_scorege50/CASECONTROL_SingleBuffer"
receiver_dat_name_main   <- "births_main_analysis_df"
receiver_dat_name_filter <- paste0(receiver_dat_name_main, 
                                   "_case", case_ratio, "control", control_ratio, 
                                   "_", outcome_var_raw, "_", ptb_name, ".csv")

###############
# Load data
##############
# Load births dataset
receiver_cc <- fread(file.path(server_cleandata_path, receiver_dat_name_filter))

# Load CAFO facilities
cafo_data_name <- "facilities20250212_cleaned.csv"
source <- fread(file.path(transfer_data_path, cafo_data_name))

# Load spatial CAFO data
source_sf <- readRDS(file.path(transfer_data_path, "facilities20250212_cleaned_sf.rds"))
source_sf <- st_transform(source_sf, crs = crs_projected_state)

#################################
# Prepare Data
#################################
# Filter CAFOs active during study period
source_select <- source %>% filter(
  construction_upper_bound < 2007 & 
    (is.na(destruction_lower_bound) | destruction_lower_bound > (2011 + 1))
)

source_select_sf <- source_sf %>% filter(
  construction_upper_bound < 2007 & 
    (is.na(destruction_lower_bound) | destruction_lower_bound > (2011 + 1))
)

####################################################
# Setup Model Inputs for SingleBuffer fxn
####################################################
# Response (y, trials)
y      <- receiver_cc$case_numeric
trials <- rep(1, times = length(y))

# Covariate matrix (X)
# convert factor variables to factors, scale continuous variables
x <- receiver_cc %>%
  dplyr::select(all_of(covariates)) %>%
  mutate(across(all_of(covariates_cat), as.factor), 
         across(all_of(covariates_continuous), ~ as.vector(scale(.x)))) 
# adds intercept column (needed for EpiBuffer models)
x_model_mat <- model.matrix( ~ ., data = x)

# Unique receiver (births) locations
n_ind            <- nrow(receiver_cc)
unique_locations <- unique(paste0(receiver_cc$lat, ",", receiver_cc$long))
n_unique         <- length(unique_locations)

# Map location to index (v vector)
location_keys  <- paste0(receiver_cc$lat, ",", receiver_cc$long)
# Map each location key to its column index
location_index <- setNames(seq_along(unique_locations), unique_locations)
v              <- location_index[location_keys]
names(v)       <- receiver_cc$birth_id

# Add numeric location ID to dataframe
receiver_cc <- receiver_cc %>%
  mutate(location = paste0(lat, ",", long),
         location_m_id = as.numeric(factor(location))) %>%
  dplyr::select(-location)

# Exposure distance matrix
receiver_sf <- st_as_sf(receiver_cc, coords = c("long", "lat"), 
                        crs = crs_unprojected_wgs84, remove = FALSE)
receiver_sf <- st_transform(receiver_sf, crs = crs_projected_state)
receiver_sf_unique <- receiver_sf %>% distinct(geometry, .keep_all = TRUE)
all(unique(location_keys) == receiver_sf_unique$loc_key) # check that ordering is same

exposure_dists <- st_distance(receiver_sf_unique, source_select_sf) / 1000
exposure_dists <- units::drop_units(exposure_dists)
filtered_exposure_dists <- exposure_dists[, !apply(exposure_dists > radius_range[2], 2, all), drop = FALSE]

# Quick nearest-CAFO check
# Get nearest source index for each receiver
nearest_idx     <- sf::st_nearest_feature(receiver_sf, source_select_sf)
nearest_sources <- source_select_sf[nearest_idx, ]
# Compute distance
distances <- sf::st_distance(receiver_sf, nearest_sources, by_element = TRUE)
receiver_sf$nearest_cafo_approx <- drop_units(distances)/1000
receiver_sf$source_id_approx <- nearest_sources$source_id
receiver_sf$nearest_cafo_diff <- receiver_sf$nearest_cafo_km - receiver_sf$nearest_cafo_approx
summary(receiver_sf$nearest_cafo_diff)
table(receiver_sf$source_id == receiver_sf$source_id_approx)
temp <- receiver_sf %>% filter(source_id != source_id_approx)

# Compile all model input 
model_input <- list(
  mcmc_samples = mcmc_samples,
  y = as.matrix(y),
  x = x_model_mat,
  v = v,
  radius_range = radius_range,
  exposure_definition_indicator = exposure_definition_indicator,
  exposure_dists = filtered_exposure_dists,
  metrop_var_radius = metrop_var_radius,
  likelihood_indicator = likelihood_indicator,
  trials = trials
)

########################################
# Setup logging print to screen results
########################################
logs_dir <- file.path(server_outdata_path, save_results_folder)
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)

# Create log file and redirect console output
log_file <- file.path(logs_dir, paste0("SINGLEBR_", 
                                       "case", case_ratio, "control", control_ratio,
                                       "_exp", exposure_definition_indicator, 
                                       "_", outcome_var_raw, "_", ptb_name, 
                                       "_maxr", radius_range[2],
                                       "_seed", seed, 
                                       optional_name,  "_log.txt"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)  # Ensure sink stops even if error occurs

print(paste0("n_ind = ", n_ind))
print(paste0("n_unique_locations = ", n_unique))
print(paste0("ncol(filtered_exposure_dists = ", ncol(filtered_exposure_dists)))

########################
# Fit SingleBuffer Model
########################
print(paste0("Fitting SingleBuffer"))
set.seed(seed)
start_time <- Sys.time()
single_results <- EpiBuffer::SingleBuffer(mcmc_samples = mcmc_samples,
                                          y = as.matrix(y),
                                          x = x_model_mat,
                                          v = v,
                                          radius_range = radius_range,
                                          exposure_definition_indicator = exposure_definition_indicator,
                                          exposure_dists = filtered_exposure_dists,
                                          metrop_var_radius = metrop_var_radius,
                                          likelihood_indicator = likelihood_indicator,
                                          trials = trials,
                                          waic_info_indicator = 0)
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("Done. Time taken:", time_taken, "\n")

####################################
# Compute WAIC and MCMC Diagnostics
####################################
print("Computing WAIC")
waic <- unlist(compute_waic(single_results, burn_in = burnin, thin = thin))
print(waic)

print("Computing MCMC Diagnostics ")
mcmc_diagnostics <- mcmc_diagnostics_singlebuffer(single_results,
                                                  burn_in = burnin,
                                                  thin = thin)
print(mcmc_diagnostics)

########################
# Save full results
########################
results_file <- file.path(server_outdata_path, save_results_folder,
                          paste0("SINGLEBR_RESULTS_", 
                                 "case", case_ratio, "control", control_ratio,
                                 "_exp", exposure_definition_indicator, 
                                 "_", outcome_var_raw, "_", ptb_name, 
                                 "_maxr", radius_range[2],
                                 "_seed", seed,
                                 optional_name,   ".rds"))

all_results <- list(single_results = single_results, 
                    time = time_taken, 
                    mcmc_diagnostics = mcmc_diagnostics, 
                    metrop_var_radius = metrop_var_radius)

saveRDS(all_results, results_file)

########################
# Thinned posterior samples
########################
# theta_keep_set <- get_posterior_samples_theta(single_results$theta,
#                                               m_max = single_results$exposure_scale, 
#                                               burn_in = burnin, 
#                                               thin = thin)
# radius_keep_set <- get_posterior_samples_radius(single_results$radius,
#                                                 burn_in = burnin, 
#                                                 thin = thin)
# other_var_keep_set <- get_posterior_samples_othervar(results = single_results,
#                                                      model_input = model_input,
#                                                      exposure_indicator = model_input$exposure_definition_indicator,
#                                                      burn_in = burnin,
#                                                      thin = thin,
#                                                      radius_range =model_input$radius_range)
# 
# results_file <- file.path(server_outdata_path, save_results_folder,
#                           paste0("SINGLEBR_THINRESULTS_", 
#                                  "case", case_ratio, "control", control_ratio,
#                                  "_exp", exposure_definition_indicator, 
#                                  "_", outcome_var_raw, "_", ptb_name, 
#                                  "_maxr", radius_range[2],
#                                  optional_name,   ".rds"))
# 
# all_results_thin <- list(theta = theta_keep_set, 
#                          radius = radius_keep_set, 
#                          other_var = other_var_keep_set)
# 
# saveRDS(all_results_thin, results_file)

# STOP


