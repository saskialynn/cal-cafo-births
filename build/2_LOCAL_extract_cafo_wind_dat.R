
#...............................................................................
# Script: 2_LOCAL_extract_cafo_wind_dat.R
# Purpose: Link CAFO locations to NARR grid cells and extract wind information
#          - Use NARR u and v vectors to compute wind_speed, wind_direction_to, 
#            wind_direction_from, and dispersion wedge
#
# Data inputs: 
#   - Cleaned CAFO locations (spatial data)
#   - NARR wind data (u and v components)
#
# Data outputs: 
#   - RDS and CSV of processed wind data at CAFO locations
#
# Helper functions: 
#   - extract_NARR_uv_years()
#   - uv_to_speed_direction_wedge()
#
# Dependencies:
#   - R packages: here, dplyr, sf, data.table
#   - Custom scripts: code/setup_header.R, code/functions/tidyNARRData_fxn.R
#...............................................................................

##########################
# Setup
##########################
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
source(file.path(build_path, "functions/tidyNARRData_fxn.R"))

#########################
# Load CAFO data
#########################
cafo_dat_name <- "facilities20250212_cleaned_sf" 
cafo_dat_sf <- readRDS(file.path(transfer_data_path, 
                                 paste0(cafo_dat_name, ".rds")))
source <- cafo_dat_sf %>% st_transform(crs_projected_state)

################################################################################
# Extract u/v wind vector data, compute speed, direction and wedges
################################################################################
# Get the u and v wind vectors on each day of interest at each CAFO location
start_year <- 2006
end_year   <- 2012

# Extract u component
narr_extract_u <- extract_NARR_uv_years(narr_data_path = wind_data_path, # path to wind data files
                                        source_dat_sf = source,          # source data sf object
                                        narr_variable_name = "uwnd",     
                                        narr_year_start = start_year, 
                                        narr_year_end = end_year)    
# Extract v component
narr_extract_v <- extract_NARR_uv_years(narr_data_path = wind_data_path, # path to wind data files
                                        source_dat_sf = source,          # source data sf object
                                        narr_variable_name = "vwnd",     
                                        narr_year_start = start_year, 
                                        narr_year_end = end_year)

# Compute wind speed, direction, and dispersion wedge
wedge_angle_degrees <- 90 # total wedge angle
narr_wind <- uv_to_speed_direction_wedge(
  narr_extract_u = narr_extract_u,
  narr_extract_v = narr_extract_v,
  wedge_angle   = wedge_angle_degrees / 2 # function expects half-angle
  # function expects degrees to be added to either side of the wind vector; 
  # so divide by two when inputting into function
)

# Remove duplicate columns, clean names
narr_wind <- narr_wind %>%
  dplyr::select(-ends_with(".y")) %>%
  dplyr::rename_with(~gsub(".x", "", .x), ends_with(".x")) 

##########################
# Save processed wind data
##########################
save_name <- paste0("narrwind_", 
                    start_year, "to", end_year, "_", 
                    "wedgeangle", wedge_angle_degrees, "_", 
                    cafo_dat_name)

saveRDS(narr_wind, file.path(transfer_data_path, paste0(save_name, ".rds")))
fwrite(narr_wind, file.path(transfer_data_path, paste0(save_name, ".csv")), row.names = FALSE)

# STOP

