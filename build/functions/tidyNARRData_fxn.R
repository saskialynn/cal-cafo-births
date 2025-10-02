
#..............................................................................
# Script:  tidyNARRData_fxn.R
#
# Purpose: Helper functions to extract, process, and compute wind data 
#          from NARR u/v vectors for CAFO locations.
#
# Dependencies: raster, purrr, sf, dplyr, tidyr
#..............................................................................

# Load Package Dependencies
library(raster)
library(purrr)
library(sf)
library(dplyr)

#-------------------------------------------------------------------------------
# Function: extract_NARR_uv()
# Purpose: Extract a single year of NARR u or v wind data at CAFO locations
# Inputs:
#   narr_data_path: path to wind data file (.nc)
#   source_dat_sf: sf object of source locations (must have source_id column)
#   narr_variable_name: "uwnd" or "vwnd"
#   narr_year: year to extract (depends on which data you have obtained)
# Output:
#   narr_data_long: tidied long-format dataframe
# Dependencies: raster
#-------------------------------------------------------------------------------
extract_NARR_uv <- function(narr_data_path,
                            source_dat_sf,
                            narr_variable_name, 
                            narr_year) { 
  
  # Import NARR data as raster stack
  narr_stack <- raster::stack(narr_data_path)
  
  # Extract narr_data values for each raster at each source site and tidy
  narr_extract <- raster::extract(narr_stack, 
                                  st_transform(source_dat_sf, projection(narr_stack)), 
                                  cellnumbers = TRUE) %>%
    as_tibble() %>% 
    dplyr::rename(narr_cells = cells) %>%
    dplyr::rename_with(~gsub("X", "", .x), starts_with("X")) # remove leading X from date
  
  # Extract NARR grid cooordinates (centroid coordinates for each cell)
  cell_coords <- raster::xyFromCell(narr_stack, narr_extract$narr_cells) %>%
    as_tibble() %>%
    rename(narr_lon = x, narr_lat = y)
  
  # Extract source coordinates from geometry column
  source_coords <- st_coordinates(source_dat_sf) %>%
    as_tibble() %>%
    rename(source_lon = X, source_lat = Y) # X = Longitude, Y = Latitude
  
  # Project NARR coordinates to source CRS
  narr_sf_proj <- st_as_sf(cell_coords, 
                           coords = c("narr_lon", "narr_lat"), 
                           crs = st_crs(narr_stack)) %>% 
    st_transform(crs(source_dat_sf)) 
  
  # Compute distance from each source to NARR cell
  cell_to_source_dist <- st_distance(source_dat_sf, narr_sf_proj, by_element = TRUE) %>% data.frame()
  
  cell_coords_proj <- narr_sf_proj %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(narr_lon = X, narr_lat = Y)
  
  # Combine all data into single data frame
  narr_data <- cbind(source_id = source_dat_sf$source_id, 
                     source_coords, 
                     cell_coords_proj, 
                     cell_to_source_dist_m = cell_to_source_dist[,1], 
                     narr_extract) 
  
  
  
  # Tidy output: convert data from wide to long format
  narr_data_long <- narr_data %>%
    tidyr::pivot_longer(-c(source_id, narr_cells, 
                           narr_lon, narr_lat, source_lon, source_lat, cell_to_source_dist_m), 
                        names_to  = "full_date", 
                        values_to = as.character(narr_variable_name)) %>%
    dplyr::mutate(date = as.Date(full_date, "%Y.%m.%d")) %>% # rename the columns to Y.m.d format
    dplyr::select(-full_date)
  
  return(narr_data_long)
}

#-------------------------------------------------------------------------------
# Function: extract_NARR_uv_years()
# Purpose: Extract multiple years of NARR u and v wind data
# Inputs: 
#   narr_data_path:     path to wind data directory
#   source_dat_sf:      source location sf object; must have source_id column
#   narr_variable_name: variable to extract, "uwnd" or "vwnd"
#   narr_year_start:    first year to extract
#   narr_year_end:      last year to extract
# Output:
#   narr_wind_temp: combined tidied dataframe for all years
# Dependencies:
#   extract_NARR_uv()
#-------------------------------------------------------------------------------
extract_NARR_uv_years <- function(narr_data_path, 
                                  source_dat_sf, 
                                  narr_variable_name, 
                                  narr_year_start, 
                                  narr_year_end){
  narr_wind_temp <- tibble()
  
  # Loop through all years and tidy the NARR data
  for (year in c(narr_year_start:narr_year_end)) {
    
    if(narr_variable_name == "uwnd"){
      narr_data_path_year <- paste0(narr_data_path, 
                                    "/wind_u_1999_2012/uwnd.10m.",
                                    as.character(year),
                                    ".nc",
                                    sep = "")
    }else if(narr_variable_name == "vwnd"){
      narr_data_path_year <- paste0(narr_data_path, 
                                    "/wind_v_1999_2012/vwnd.10m.",
                                    as.character(year),
                                    ".nc",
                                    sep = "")
    }
    
    narr_wind_temp <- bind_rows(narr_wind_temp,
                                extract_NARR_uv(narr_data_path_year,
                                                source_dat_sf,
                                                narr_variable_name,
                                                year)
    )
  }
  return(narr_wind_temp)
}

#-------------------------------------------------------------------------------
# Function: fix_angle_range()
# Purpose: Constrain angle ranges to [-180, +180] relative to North
#   (i.e. using positive y-axis as reference (0 degrees; North))
# Input: 
#   degree: angle degree in range [0, 360]
# Output:
#   degree: angle constrained to [+180, -180]
#-------------------------------------------------------------------------------

fix_angle_range <- function(degree){ # zero degrees is positive y-axis (North)
  if(degree < -180){ # if rotated > 180 degrees in counterclockwise direction, convert to clockwise rotation
    # Quadrants 1 and 2
    while(degree < -180){
      degree <- degree + 360
    }
  } else if(degree > 180){ # if rotated > 180 degrees in clockwise direction, convert to counterclockwise rotation
    # Quadrants 3 and 4
    while(degree > 180){
      degree <- degree - 360
    }
  }
  return(degree)
}

#-------------------------------------------------------------------------------
# Function: uv_to_speed_direction()
# Purpose: Compute wind speed and direction from u and v vector data
# Inputs:
#   narr_extract_u: dataframe with columns narr_cells, source_id, uwnd, date  
#   narr_extract_v: dataframe with columns narr_cells, source_id, vwnd, date
#   wedge_angle: half-angle (degrees) for wedge calculation
# Output:
#   narr_wind: merged dataframe with uwnd and vwnd vectors at each source_id and date, 
#              with computed wind direction and speed metrics added: 
#               narr_wind_speed
#               narr_wind_direction_from: meteorologic convention, direction wind is blowing from
#               narr_wind_direction_t: meterologic convention, direction wind is blowing to
# Dependencies: dplyr, purrr
#-------------------------------------------------------------------------------
uv_to_speed_direction_wedge <- function(narr_extract_u, 
                                        narr_extract_v, 
                                        wedge_angle){
  narr_wind <- left_join(narr_extract_u, narr_extract_v, by = c("source_id", "date")) %>%
    dplyr::mutate(
      narr_wind_speed = sqrt(narr_extract_u$uwnd^2 + narr_extract_v$vwnd^2),
      
      # STANDARD METEROLOGIC CONVENTION 
      # atan2(v, u) = arctangent of v / u: computes angle in radians between (u,v) vector and positive x-axis
      # 180/pi converts radians to degrees
      # (270 - ...) aligns with meteorological convention for direction from
      # modulo 360 to ensure values are always in 0 to 360 degree range
      narr_wind_direction_from = (
        270 - ((180 / pi) * (atan2(narr_extract_v$vwnd, narr_extract_u$uwnd))))%%360, 
      narr_wind_direction_to = (narr_wind_direction_from + 180)%%360,
      # compute wedge min / max using meteorologic convention
      wedge_max_meteorologic = (narr_wind_direction_to + wedge_angle)%%360,
      wedge_min_meteorologic = (narr_wind_direction_to - wedge_angle)%%360, 
      
      # MATHEMATICAL CONVENTION 
      narr_wind_direction_to_math = 270 - narr_wind_direction_from, 
      wedge_max_math = (narr_wind_direction_to_math + wedge_angle)%%360, 
      wedge_min_math = (narr_wind_direction_to_math - wedge_angle)%%360,
      
      # ADJUSTED METEROLOGIC CONVENTION
      wind_direction_to_wedge_calc = atan2(narr_extract_u$uwnd,
                                           narr_extract_v$vwnd)*(180/pi),
      
      wind_direction_from_wedge_calc = (wind_direction_to_wedge_calc + 180) %>% 
        purrr::map(fix_angle_range) %>% unlist(), 
      # compute wedge min / max allowing counterclockwise rotation from 0 degrees
      wedge_max_wedge_calc = (wind_direction_to_wedge_calc + wedge_angle) %>% 
        purrr::map(fix_angle_range) %>% unlist(),
      wedge_min_wedge_calc = (wind_direction_to_wedge_calc - wedge_angle) %>% 
        purrr::map(fix_angle_range) %>% unlist()
    )
  return(narr_wind)
}

#..............................................................................
# NOTES on WIND DIRECTION CONVENTIONS
# Standard meteorologic convention:
#   angles are relative to y-axis (North = zero degrees) and increase clockwise (always positive)
#     North = 0 degrees, East = 90 degrees, South = 180 degrees, West = 270 degrees

# Mathematical convention: 
#   angles are relative to the x-axis (angle = 0) (so that tan(angle) = v/u)
#   and increase (+ direction) counterclockwise moving from east (x-axis, 0) to north (y-axis, 90)
#   Use: math convention is needed for sin/cos calculations 
#         and is mostly useful for visualization/plotting on axes

# Adjusted meteorologic reference: 
#   Meteorologic convention adjusted to allow for counter-clockwise rotation (negative angles)
#   this is useful for later wedge calculations and avoids issues in 4th quadrant
#   for 1st and 2nd quadrants, equals meteorologic convention (positive angles, clockwise rotation from +y axis)
#   for 3rd and 4th quadrants, angles are negative and rotate counterclockwise from +y axis
#   clockwise direction is positive; counterclockwise is negative
# Constrained to between +180 and -180: 
#     North (+y axis) =   0 degrees
#     East (+x axis)  = +90 degrees
#     South(-y axis)  = +180 degrees
#     West (-x axis)  = -90 degrees

# atan2(y, x): returns the angle (in radians) between the positive x-axis and 
#              the vector from the origin to (x, y). East-Counterclockwise convention

# REFERENCES for wind speed + direction calculation: 
# Meteorologic convention direction calculation: 
#     www.eol.ucar.edu/content/wind-direction-quick-reference
#     confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
# Converting between Math and Meteorologic convention: 
# www.e-education.psu.edu/meteo300/node/719#:~:text=The%20math%20angle%20is%20measured,degrees%20minus%20the%20math%20angle.
# Calculator for quick checks: www.cactus2000.de/uk/unit/masswin.shtml
# varying conventions and applications of atan2: https://en.wikipedia.org/wiki/Atan2
#..............................................................................







