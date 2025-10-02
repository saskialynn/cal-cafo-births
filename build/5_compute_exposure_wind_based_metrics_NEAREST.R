
#...............................................................................
# Script: 5_compute_exposure_wind_based_metrics_NEAREST.R
#
# Purpose: 
#   Compute wind-based exposure metrics for NEAREST CAFO to each receiver.
#
# Data Inputs:
#   - Births dataset (cleaned)
#   - CAFO dataset (cleaned)
#   - Nearest CAFO (distance & bearing) for each receiver
#   - NARR wind data (direction) for each CAFO
#
# Data Outputs:
#   - "..._NEARESTCAFO_ALLRECEIVERS_WINDEXP.csv": 
#     Births data with nearest CAFO and wind exposure metrics
#
# Helper Functions: 
#   - get_receiver_downwind_multiday_FAST()
#
# Dependencies: 
#   - R packages: here, dplyr, sf, data.table, future, furrr, tools, progressr, lubridate
#   - Custom scripts: 
#       - code/setup_header.R
#       - code/functions/get_receivers_downwind_fxn.R
#
#...............................................................................

##########################
# Setup
##########################
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
source(file.path(build_path, "functions/get_receivers_downwind_fxn.R"))

library(future)
library(furrr)
library(tools)
library(progressr)
library(lubridate)

######################
# Specify Parameters
######################
birth_dat_name <- "births_main_analysis_df.csv"
cafo_dat_name  <- "facilities20250212_cleaned_sf.rds"
dist_mat_name  <- "births_main_analysis_df_cafo20250212_NEARESTCAFO_DISTANGLE.csv"
narr_wind_name <- "narrwind_2006to2012_wedgeangle90_facilities20250212_cleaned_sf.csv"

cafo_dat_id    <- "20250212"
max_distance   <- 25000 # meters

###############
# Load Data
###############
# Receiver: Births Data
receiver <- data.table::fread(file.path(server_cleandata_path, birth_dat_name))

# Wind data for each CAFO and date
narr_wind <- fread(file.path(transfer_data_path, narr_wind_name)) 
# Extract wedge angle from file name
wedge_angle <- regmatches(narr_wind_name, 
                          regexpr("(?<=wedgeangle)\\d+", narr_wind_name, perl = TRUE))

# Source/Receiver distance & angle matrix (nearest CAFO)
source_receiver_mat <- fread(file.path(server_cleandata_path, dist_mat_name)) %>%
  as.data.table() %>% 
  filter(birth_id %in% receiver$birth_id)

################################################################################
# Define Potentially Exposed Receivers based on max_distance
################################################################################
receiver$obs_id <- paste0(receiver$receiver_id, "_", receiver$birth_id)
source_receiver_mat$obs_id <- paste0(source_receiver_mat$receiver_id, "_", source_receiver_mat$birth_id)

# Exposed population: Source/receiver pairs that are within max_distance
# (already accounting for cafo construction/destruction dates relative to birthdate)
source_receiver_mat_exposed <- source_receiver_mat %>%
  dplyr::filter(distance_m <= max_distance) %>%
  mutate(source_birth_id = paste0(source_id, "_", birth_id))
obs_id_exposed <- unique(source_receiver_mat_exposed$obs_id)

# Restrict wind data to relevant sources
narr_wind_small <- narr_wind %>% 
  dplyr::filter(source_id %in% source_receiver_mat_exposed$source_id) %>%
  dplyr::select(source_id, date, 
                wind_direction_to_wedge_calc,
                wedge_min_wedge_calc, wedge_max_wedge_calc,
                wedge_min_meteorologic, wedge_max_meteorologic) %>%
  as.data.table()

# Restrict receivers to exposed subset
receiver_exp_only <- receiver %>% 
  dplyr::filter(obs_id %in% obs_id_exposed)

###################################
# Create list of receiver tasks 
###################################
#..................ONLY RE-RUN IF NECESSARY (slow)..................#
# chunk distance / angle information into list elements by obs_id
# contains all valid (based on dates and distance cutoff) sources for each receiver

# source_receiver_mat_exposed <- source_receiver_mat_exposed %>%
#   dplyr::select(obs_id, receiver_id, birth_id, source_id,
#                 BDATE, gestob, doc,
#                 distance_m,  distance_km, bearing_source_to_receiver)
# 
# source_receiver_dist_angle_list <- split(source_receiver_mat_exposed,
#                                          by = "obs_id")
# 
# # Each task = one receiver/birth with valid CAFO(s) + daily date sequence
# obs_task_list <- purrr::map(
#   .x = unname(source_receiver_dist_angle_list),
#   .f = function(x) {
#     list(
#       obs_id = as.character(x$obs_id[1]),
#       dist_angle = as.data.table(x),
#       date_seq = seq(x$doc[1], x$BDATE[1], by = "day")
#     )
#   }
# )
# saveRDS(obs_task_list, 
#         file.path(server_cleandata_path,
#                   "main_analysis_df_obs_task_list_intermediate_for_wind_metrics.rds"))

#.....................................................................#

#..................LOAD DIRECTLY if ALREADY CREATED..................#
obs_task_list <- readRDS(file.path(server_cleandata_path,
                                   "main_analysis_df_obs_task_list_intermediate_for_wind_metrics.rds"))

#########################################
# Compute Wind Exposure for NEAREST CAFO
#########################################
plan(multisession, workers = 20)
handlers(global = TRUE)  

with_progress({
  p <- progressor(along = obs_task_list)
  
  wind_results <- future_map_dfr(
    obs_task_list,
    function(task) {
      p()  # update progress
      get_receiver_downwind_multiday_FAST(task, narr_wind_small)
    },
    .options = furrr_options(seed = TRUE)
  )
})
plan(sequential)

save_name <- paste0(tools::file_path_sans_ext(birth_dat_name),
                    "_cafo", cafo_dat_id,
                    "_wedgeangle", wedge_angle,
                    "_maxdist", max_distance, "_ALLWINDEXP.csv")

data.table::fwrite(wind_results, file.path(server_cleandata_path, save_name))

############################
# Compute Summary Statistics
###########################
downwind_df_name <- paste0(
  tools::file_path_sans_ext(birth_dat_name),
  "_cafo", cafo_dat_id,
  "_wedgeangle", wedge_angle,
  "_maxdist", max_distance,
  "_ALLWINDEXP.csv"
)
downwind_df <- data.table::fread(file.path(server_cleandata_path, downwind_df_name))

# check that downwind_df only contains valid source/receiver pairs
length(unique(paste0(downwind_df$source_id, "_", downwind_df$birth_id)))

# Compute gestational age and trimester for each date
downwind_df <- downwind_df %>%
  mutate(
    BDATE = as.Date(BDATE),
    date = as.Date(date),
    doc = as.Date(doc),

    # start week count at 1
    gestational_age_days = as.numeric(date - doc),
    gestational_age_weeks = floor(gestational_age_days / 7) + 1,

    trimester = case_when(
      gestational_age_weeks >= 1  & gestational_age_weeks <= 13 ~ "tri1",
      gestational_age_weeks >= 14 & gestational_age_weeks <= 26 ~ "tri2",
      gestational_age_weeks >= 27                               ~ "tri3",
      TRUE ~ NA_character_
    )
  )

# Trimester 3: last 30 days before delivery
downwind_df_filtered <- downwind_df %>%
  filter(
    (trimester %in% c("tri1", "tri2")) |
      (trimester == "tri3" & date >= (BDATE - 30) & date <= BDATE)
  )

# Count days of exposure by obs_id, trimester, wind exposure bin
wind_days_by_trimester <- downwind_df_filtered %>%
  count(obs_id, trimester, BINNED_wind_exposure, name = "days") %>%
  pivot_wider(
    names_from = c(trimester, BINNED_wind_exposure),
    values_from = days,
    values_fill = 0
  )

# Join exposure days back to receiver/source matrix
receiver_w_exp <- left_join(source_receiver_mat, 
                            wind_days_by_trimester, 
                            by = "obs_id")
receiver_w_exp <- receiver_w_exp %>% mutate(nearest_cafo_km = distance_m / 1000)  

# Compute downwind proportions
# NOTE: if nearest CAFO is beyond threshold, wind exposure is NA
table(wind_na = is.na(receiver_w_exp$tri1_Downwind), 
      dist_g25 = (receiver_w_exp$nearest_cafo_km > 25))

receiver_w_exp <- receiver_w_exp %>%
  mutate(
    gestation_days_total = as.numeric(as.Date(BDATE) - as.Date(doc)),
    gestation_weeks_total = floor(gestation_days_total / 7) + 1,
    
    # Logical indicators for each trimester (was there exposure opportunity?)
    tri1_reached = gestation_weeks_total >= 1,
    tri2_reached = gestation_weeks_total >= 14,
    tri3_reached = gestation_weeks_total >= 27,
    # Days in each trimester
    tri1_total_days = rowSums(across(c("tri1_Downwind", "tri1_Upwind", "tri1_AlmostDownwind",
                                       "tri1_AlmostUpwind", "tri1_Crosswind")), na.rm = TRUE),
    tri2_total_days = rowSums(across(c("tri2_Downwind", "tri2_Upwind", "tri2_AlmostDownwind",
                                       "tri2_AlmostUpwind", "tri2_Crosswind")), na.rm = TRUE),
    tri3_total_days = rowSums(across(c("tri3_Downwind", "tri3_Upwind", "tri3_AlmostDownwind",
                                       "tri3_AlmostUpwind", "tri3_Crosswind")), na.rm = TRUE),
    # Proportion of days downwind, by trimester
    tri1_prop_downwind = tri1_Downwind / tri1_total_days,
    tri2_prop_downwind = tri2_Downwind / tri2_total_days,
    tri3_prop_downwind = tri3_Downwind / tri3_total_days, 
    tri12_prop_downwind = (tri1_Downwind + tri2_Downwind) / (tri1_total_days + tri2_total_days), 
    tri123_prop_downwind = (tri1_Downwind + tri2_Downwind + tri3_Downwind) / (tri1_total_days + tri2_total_days + tri3_total_days)
  )

# Checks
table(receiver_w_exp$tri1_total_days)
table(receiver_w_exp$tri2_total_days)
table(receiver_w_exp$tri3_total_days)
summary(receiver_w_exp$tri1_prop_downwind) # NAs are people beyond threshold distance
summary(receiver_w_exp$tri2_prop_downwind)
summary(receiver_w_exp$tri3_prop_downwind) # NAs are people beyond threshold distance OR birth occurred before 3rd trimester
summary(receiver_w_exp$tri12_prop_downwind)
summary(receiver_w_exp$tri123_prop_downwind)

# Convert NA proportions to 0 (NA means they were outside distance range, not exposed)
receiver_w_exp <- receiver_w_exp %>%
  mutate(tri1_prop_downwind = ifelse(is.na(tri1_prop_downwind), 0, tri1_prop_downwind),
         tri2_prop_downwind = ifelse(is.na(tri2_prop_downwind), 0, tri2_prop_downwind), 
         tri3_prop_downwind = ifelse(is.na(tri3_prop_downwind), 0, tri3_prop_downwind),
         
         tri1_prop_upwind = tri1_Upwind/tri1_total_days, 
         tri2_prop_upwind = tri2_Upwind/tri2_total_days, 
         tri3_prop_upwind = tri3_Upwind/tri3_total_days, 
         tri1_prop_upwind = ifelse(is.na(tri1_prop_upwind), 1, tri1_prop_upwind), 
         tri2_prop_upwind = ifelse(is.na(tri2_prop_upwind), 1, tri2_prop_upwind), 
         tri3_prop_upwind = ifelse(is.na(tri3_prop_upwind), 1, tri3_prop_upwind))

# Save final dataset
save_name <- paste0(tools::file_path_sans_ext(birth_dat_name),
                    "_cafo", cafo_dat_id,
                    "_wedgeangle", wedge_angle,
                    "_maxdist", max_distance, "_NEARESTCAFO_ALLRECEIVERS_WINDEXP.csv")

data.table::fwrite(receiver_w_exp, file.path(server_cleandata_path, save_name))

# STOP

