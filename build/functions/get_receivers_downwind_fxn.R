
#...............................................................................
# Script: code/functions/get_receivers_downwind_fxn.R
#
# Purpose: 
#   Functions to compute whether a receiver (birth) is downwind of a source (CAFO) 
#   for a single day or over multiple days of gestation.
#
# Functions:
#   - get_receiver_downwind_multiday_FAST(): optimized multi-day exposure
#
# Dependencies:
#   - R packages: dplyr, data.table
#...............................................................................

# -----------------------------------------------------------------------------
# FUNCTION: get_receiver_downwind_multiday_FAST ----
# PURPOSE:
#   Compute downwind exposures over multiple days for a receiver.
#   Optimized version using data.table. 
# INPUT:
#   - task: list with obs_id, dist_angle (source-receiver info), date_seq
#   - narr_wind_small: data.table with NARR wind for relevant sources/dates
# OUTPUT:
#   - data.table with daily exposure and binned wind categories
# -----------------------------------------------------------------------------
get_receiver_downwind_multiday_FAST <- function(task, narr_wind_small) {
  dist_angle <- task$dist_angle
  date_seq <- task$date_seq
  
  # Filter narr_wind_small for relevant dates and sources
  wind_select <- narr_wind_small[
    date %in% date_seq & source_id %in% dist_angle$source_id
  ]
  
  if (nrow(wind_select) == 0) return(NULL)
  
  # Left join with data.table - Perform join with dist_angle
  setkey(dist_angle, source_id)
  joined <- merge(wind_select, dist_angle, by = "source_id", 
                  all.x = TRUE # LEFT JOIN
  )
  
  # Compute angle difference
  joined[, wind_angle_diff := ((wind_direction_to_wedge_calc - bearing_source_to_receiver + 180) %% 360) - 180]
  
  # Binary downwind flag
  joined[, BINARY_downwind_exposed := as.integer(
    bearing_source_to_receiver >= wedge_min_wedge_calc & 
      bearing_source_to_receiver <= wedge_max_wedge_calc
  )]
  
  # Categorical binned wind exposure
  joined[, BINNED_wind_exposure := 
           fifelse(wind_angle_diff > -22.5 & wind_angle_diff <= 22.5, "Downwind",
                   fifelse(wind_angle_diff > 22.5  & wind_angle_diff <= 67.5, "AlmostDownwind",
                           fifelse(wind_angle_diff > 67.5  & wind_angle_diff <= 112.5, "Crosswind",
                                   fifelse(wind_angle_diff > 112.5 & wind_angle_diff <= 157.5, "AlmostUpwind",
                                           fifelse(wind_angle_diff > 157.5 | wind_angle_diff <= -157.5, "Upwind",
                                                   fifelse(wind_angle_diff > -157.5 & wind_angle_diff <= -112.5, "AlmostUpwind",
                                                           fifelse(wind_angle_diff > -112.5 & wind_angle_diff <= -67.5, "Crosswind",
                                                                   fifelse(wind_angle_diff > -67.5 & wind_angle_diff <= -22.5, "AlmostDownwind", 
                                                                           NA_character_)
                                                           )))))))]
  
  return(joined)
}

# STOP





