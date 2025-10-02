
#...............................................................................
# Script: 3_prepare_births_dat.R
#
# Purpose: 
#   Clean and prepare births data for main and sensitivity analyses.
#   - Create 3-level PTB outcome variable
#   - Estimate date of conception
#   - Apply data restrictions (nulliparous / no duplicates / geocode score)
#
# Data Inputs:
#   - Births dataset (raw)
#
# Data Outputs: 
#   - "births_data.csv": Cleaned births dataset
#   - "births_main_analysis_df.csv": Main analysis dataset (covariates restricted)
#   - "births_main_analysis_nc_df.csv": Sensitivity dataset (covariates restricted)
#
# Dependencies: 
#   - R packages: here, dplyr, data.table, tidyr, stringr
#   - Custom scripts: code/setup_header.R
#...............................................................................

##########################
# Setup
##########################
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))

##########################
# Specify Data
##########################
birth_dat_name <- "births_data.csv"
cafo_dat_name  <- "facilities_2025-02-12_cleaned_sf.rds"
cafo_dat_id    <- "20250212"

# load births data
births_dat <- data.table::fread(file.path(server_rawdata_path, birth_dat_name))

##################################################
# Create 3-level PTB outcome variable (prem3s)
##################################################
# 1:	20-31 weeks (very PTB); 2:	32-36 weeks (moderate PTB); 3:	37-41 weeks (term)
table(prem3 = births_dat$prem3, prem5s = births_dat$prem5s, useNA = "always")

births_dat <- births_dat %>%
  mutate(
    prem3s = case_when(
      prem3 == 1 & !is.na(prem5s) ~ 1,
      prem3 == 2 & !is.na(prem5s) ~ 2,
      prem3 == 3 & !is.na(prem5s) ~ 3,
      TRUE ~ NA_real_
    ) %>% as.factor()
  )

table(prem3s = births_dat$prem3s, prem5s = births_dat$prem5s, useNA = "always")

###################################
# Estimate Date of Conception
###################################
# doc = bdate - (gestob x 7) + 14 
births_dat <- births_dat %>%
  mutate(
    BDATE = as.Date(BDATE, format = "%m/%d/%Y"),
    doc = BDATE - (gestob * 7) + 14
  )

###############################################
# Check geocode scores and non-unique locations
###############################################
# Count multiple "NO" matches in geocode 
births_dat %>%
  filter(`_SCORE_` >= 50) %>%
  mutate(no_count = str_count(`_NOTES_`, "NO")) %>%
  count(no_count)

births_dat %>%
  mutate(no_count = str_count(`_NOTES_`, "NO")) %>%
  count(no_count)

# Threshold for excluding poor geocodes
threshold <- births_dat %>%
  mutate(no_count = str_count(`_NOTES_`, "NO")) %>%
  filter(no_count >= 3) %>%
  summarise(max_score = max(`_SCORE_`, na.rm = TRUE)) %>%
  pull(max_score)

# Check zip code matches
zipcode_match <- births_dat %>%
  filter(`_NOTES_` == "ZC") %>%
  summarise(mean_score = mean(`_SCORE_`, na.rm = TRUE))

######################
# Rename columns
######################
births_dat <- births_dat %>%
  rename(receiver_id = `_brthidHST`, 
         birth_id = `_brthid`)

fwrite(births_dat, file.path(server_cleandata_path, birth_dat_name))

################################################################################
# Create and Save Restricted Births datasets
################################################################################

################################################
# Universal restrictions (used for all analysis)
################################################
births_restriction <- c("nulliparous", "dropduplicates", "scorege")
score_cutoff <- 50

# Filter: Nulliparous women 
if(births_restriction[1] == "nulliparous"){
  births_dat <- births_dat %>% filter(parity2gp == 1)
}

# Check for duplicate maternal id's
length(unique(births_dat$receiver_id))
# births_dat_dupes <- births_dat  %>%
#   group_by(receiver_id) %>%
#   mutate(dup_count = n()) %>%
#   filter(dup_count > 1) %>%
#   ungroup()

# Filter: Remove duplicate maternal IDs (not consistent with nulliparous)
if(births_restriction[2] == "dropduplicates"){
  births_dat <- births_dat %>% 
    group_by(receiver_id) %>%
    filter(n() == 1) %>%
    ungroup()
}

# Filter: Geocode score threshold
if(births_restriction[3] == "scorege"){
  births_dat <- births_dat %>% dplyr::filter(`_SCORE_` >= score_cutoff)
  births_restriction[3] <- paste0(births_restriction[3], score_cutoff)
}

# save_name <- paste0("births", 
#                     "_", births_restriction[1], 
#                     "_", births_restriction[2], 
#                     "_", births_restriction[3],
#                     ".csv")
# 
# fwrite(births_dat, file.path(server_cleandata_path, save_name))

##########################
# MAIN ANALYSIS DATAFRAME
##########################
covariates_continuous <- c("mothage")
covariates_cat        <- c("edum")
strat_var             <- c("racem", "sex")

covariates <- c(covariates_continuous, covariates_cat)
covariates_eq <- paste0(covariates, collapse = "+")

# Drop observations with missing outcomes or covariates
main_analysis_df <- births_dat %>% filter(!is.na(prem3s))
main_analysis_df <- main_analysis_df %>% drop_na(all_of(covariates))
main_analysis_df <- main_analysis_df %>% drop_na(all_of(strat_var))

# Binary case indicator
main_analysis_df <- main_analysis_df %>%
  mutate(
    case         = factor(ifelse(prem2 == 1, 1, 0), levels = c(0, 1)),
    case_numeric = ifelse(case == "1", 1, 0)
  )

save_name <- "births_main_analysis_df.csv"
fwrite(main_analysis_df, file.path(server_cleandata_path, save_name))

##################################################
# NO COMPLICATIONS SENSITIVITY ANALYSIS DATAFRAME
##################################################
# Quick count
temp <- births_dat %>% filter(!is.na(prem3s))
temp <- temp %>% drop_na(all_of(c("mothage", "edum", "racem", "sex")))
table(temp$prem5snc, useNA = "always")

births_restriction <- c("nulliparous", 
                        "dropduplicates", 
                        paste0("scorege", score_cutoff), 
                        "uncomplicated")

# Filter out maternal complications
if(births_restriction[4] == "uncomplicated"){
  births_dat <- births_dat %>% 
    dplyr::filter(!is.na(prem5snc))
}

# save_name <- paste0("births", 
#                     "_", births_restriction[1], 
#                     "_", births_restriction[2], 
#                     "_", births_restriction[3],
#                     "_", births_restriction[4],
#                     ".csv")
# 
# fwrite(births_dat, file.path(server_cleandata_path, save_name))

# Drop observations with missing outcomes or covariates
main_analysis_nc_df <- births_dat %>% filter(!is.na(prem3s))
main_analysis_nc_df <- main_analysis_nc_df %>% drop_na(all_of(covariates))
main_analysis_nc_df <- main_analysis_nc_df %>% drop_na(all_of(strat_var))

# Binary case indicator
main_analysis_nc_df <- main_analysis_nc_df %>%
  mutate(
    case         = factor(ifelse(prem2 == 1, 1, 0), levels = c(0, 1)),
    case_numeric = ifelse(case == "1", 1, 0)
  )

save_name <- "births_main_analysis_nc_df.csv"
fwrite(main_analysis_nc_df, file.path(server_cleandata_path, save_name))

# STOP







