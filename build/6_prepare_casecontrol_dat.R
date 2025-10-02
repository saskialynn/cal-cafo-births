
#...............................................................................
# Script: prepare_casecontrol_dat.R
# Purpose:      
#   Prepare case–control datasets for Bayesian sensitivity analysis
#
# Data inputs:  
#   - Births dataset with CAFO exposure and covariates
#   - Births dataset (filtered set of receivers for main analysis)
#
# Data outputs: 
#   - Case–control dataset (CSV), stratified by preterm birth category
#     Filename format: 
#     "<births_main_analysis_df>_case{ratio}control{ratio}_{outcome_var}_{ptb_name}.csv"
#
# Dependencies: 
#   - R packages: here, dplyr, data.table, tools
#   - Custom scripts: code/setup_header.R
#...............................................................................

##########################
# Setup
##########################
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
library(tools)

#################################
# Specify Parameters
##################################
outcome_var_raw <- "prem3s"
ptb_level       <- 2       # 1 = early PTB, 2 = moderate PTB
case_ratio      <- 1
control_ratio   <- 5
seed            <- 123

covariates_continuous <- c("mothage")
covariates_cat        <- c("edum")
strat_var             <- c("racem", "sex")

covariates    <- c(covariates_continuous, covariates_cat)
covariates_eq <- paste0(covariates, collapse = "+")

if(ptb_level == 1){ptb_name = "early"}
if(ptb_level == 2){ptb_name = "moderate"}

#################################
# Load data
#################################
receiver_dat_name <- "births_main_analysis_df_cafo20250212_wedgeangle90_maxdist25000_NEARESTCAFO_ALLRECEIVERS_WINDEXP.csv"
optional_name <- NULL # NULL, "_NC"
receiver_filter_name <- "births_main_analysis_df.csv"

receiver <- fread(file.path(server_cleandata_path, receiver_dat_name))
receiver_filter <- fread(file.path(server_cleandata_path, receiver_filter_name))
# Restrict to main analysis sample
receiver <- receiver %>% dplyr::filter(birth_id %in% receiver_filter$birth_id)

#################################
# Prepare Data
#################################
# Drop observations missing outcome or covariates
receiver <- receiver %>% filter(!is.na(!!sym(outcome_var_raw)))
receiver <- receiver %>% drop_na(all_of(covariates))
receiver <- receiver %>% drop_na(all_of(strat_var))

# # Create indicator variable for binary case / control status
# receiver$case <- ifelse(receiver$prem2 == 1, 1, 0)
# receiver$case <- factor(receiver$case, levels = c(0, 1))
# 
# # create numeric case variable
# receiver$case_numeric <- ifelse(receiver$case == "1", 1, 0)

# Subset to controls and PTB category of interest
receiver_subset <- receiver %>% filter(case == 0 | prem3s == ptb_level)

###########################
# Sample controls
###########################
n_cases <- sum(receiver_subset$case == 1)
n_controls <- control_ratio*n_cases

# Confirm enough controls are available (should be TRUE)
nrow(receiver_subset) - n_cases > n_controls

set.seed(seed)
receiver_controls <- receiver_subset %>% 
  filter(case == 0) %>%
  slice_sample(n = n_controls, replace = FALSE)

receiver_cases <- receiver_subset %>% filter(case == 1)

# Combine into case–control dataset
receiver_cc <- rbind(receiver_cases, receiver_controls)

# Checks
nrow(receiver_cc)
length(unique(receiver_cc$birth_id))
length(unique(receiver_cc$receiver_id))

#################################
# Save Output
#################################
save_name <- paste0(file_path_sans_ext(receiver_filter_name), 
                    "_case", case_ratio, "control", control_ratio,
                    "_", outcome_var_raw, "_", ptb_name, 
                    optional_name, ".csv")

fwrite(receiver_cc, file.path(server_cleandata_path, save_name))

# STOP




