
#..............................................................................
# Script:  get_pairwise_distance_angles_fxn.R
#
# Purpose: Functions to compute geographic relationships between source (CAFO) 
#          and receiver (maternal residence) including distances and angles.
#
# Dependencies: sf, dplyr, nngeo, furrr, geosphere
#..............................................................................

library(sf)
library(dplyr)
library(nngeo)
library(furrr)
library(geosphere)

# -----------------------------------------------------------------------------#
# FUNCTION: adjust_azimuth_to_atan ----
# PURPOSE:
#   Convert azimuth angle [0, 360] to atan equivalent [-180, 180]
#   Azimuth angle is a positive clockwise angle referenced from the positive 
#     Y axis (geometry) or the North meridian (geography): North = 0
#     Azimuth has range [0, 360]
#   Function adjusts azimuth angle to (-180 to 180 degrees) range
#   where 3rd and 4th quadrants are negative rotation
# 
# INPUT:
#   azimuth_deg: numeric; value of azimuth angle
# OUTPUT:
#   atan_deg: numeric; value of atan angle
#   returns equivalent of: atan2((birth_lat-cafo_lat),(birth_lon-cafo_lon))*180/pi
# -----------------------------------------------------------------------------#

adjust_azimuth_to_atan <- function(azimuth_deg) {
  if(azimuth_deg > 180) { # 3rd or 4th quadrant: convert clockwise (+) rotation to (-) CCW rotation
    atan_deg <- azimuth_deg - 360
  }else if(azimuth_deg < -180) { # azimuth angle will actually never be negative, so this condition isn't needed
    atan_deg <- azimuth_deg + 360
  }else{
    atan_deg <- azimuth_deg
  }
  return(atan_deg)
}


# -----------------------------------------------------------------------------#
# FUNCTION: get_pairwise_distances_angles 
# PURPOSE:
#   For each source/receiver pair, calculate:
#     - Distance from source to receiver
#     - Azimuth angle (North = +Y = 0; [0, 360])
#     - Atan angle (North = +Y = 0; [-180, 180])
#     - Bearing angle using geosphere package
# INPUT:
#   source_dat_sf: sf object with source_id and geometry (point locations)
#   receiver_dat_sf: sf object with receiver_id and geometry (point locations)
#   source_dat_buffer_sf: sf object source buffers for pre-filtering
# OUTPUT:
#   source_receiver_mat: Dataframe with distances and angles for each source-receiver pair
#         azimuth angle (North = +y = 0; [0, 360] range)
#         atan angle (North = +y = 0; [-180, 180] range)
#         bearing angle: same as atan, computed with built in R function
# -----------------------------------------------------------------------------#
get_pairwise_distances_angles <- function(
    source_dat_sf_proj,
    receiver_dat_sf_proj,  
    source_dat_buffer_sf
){ 
  
  if(sf::st_crs(receiver_dat_sf_proj)$epsg != sf::st_crs(source_dat_sf_proj)$epsg){
    print("Source and Receiver SF dat must have matching CRS.")
    return(NULL)
  }
  
  # Filter to receivers within buffer distance of source
  points_within_buffer <- st_within(receiver_dat_sf_proj, source_dat_buffer_sf, 
                                    sparse = TRUE) %>% data.frame()
  receiver_select_proj <- receiver_dat_sf_proj[points_within_buffer$row.id, ]
  source_select_proj <- source_dat_sf_proj[points_within_buffer$col.id,]
  
  # Calculate pairwise distances
  print("calculating source/receiver distances")
  source_receiver_mat <- sf::st_distance(source_select_proj, 
                                         receiver_select_proj, 
                                         by_element = TRUE) %>% 
    as.data.frame()  %>%
    mutate(source_id = source_select_proj$source_id,
           receiver_id = receiver_select_proj$receiver_id, 
           birth_id = receiver_select_proj$birth_id) %>%
    rename(distance = '.')
  
  # Calculate angles (where source is reference - Angle FROM SOURCE to RECEIVER)
  
  # st_azimuth returns the planar azimuth in degrees of the target point from the origin point
  # Azimuth angle is a positive clockwise angle referenced from the 
  # positive Y axis (geometry) or the North meridian (geography):
  # North (+y) = 0; ranging between 0 and 360 clockwise from north.
  print("calculating angles")
  
  # st_azimuth does not check crs, but expects geographic coordinates
  # values returned by st_azimuth differ slightly (usually 1-5 degrees)
  # depending on whether you provide geographic or projected coordinates
  receiver_dat_sf_lonlat <- st_transform(receiver_dat_sf_proj, crs = crs_unprojected_wgs84)
  source_dat_sf_lonlat <- st_transform(source_dat_sf_proj, crs = crs_unprojected_wgs84)
  
  receiver_select_lonlat <- receiver_dat_sf_lonlat[points_within_buffer$row.id, ]
  source_select_lonlat <- source_dat_sf_lonlat[points_within_buffer$col.id,]
  
  angle_azimuth_matrix <- nngeo::st_azimuth(source_select_lonlat, receiver_select_lonlat)
  angle_atan_matrix <- sapply(angle_azimuth_matrix, adjust_azimuth_to_atan)
  
  source_receiver_mat$angle_azimuth <- angle_azimuth_matrix
  source_receiver_mat$angle_atan <- angle_atan_matrix
  
  rm(angle_azimuth_matrix)
  rm(angle_atan_matrix)
  
  # Calculate bearing (from point 1 (source) to point 2 (receiver))
  # input points should be in longitude / latitude
  # Returns angles where: North = +y = 0 degrees; [-180, 180] range
  source_coords <- st_coordinates(source_select_lonlat)
  receiver_coords <- st_coordinates(receiver_select_lonlat)
  
  bearing_matrix <- geosphere::bearing(source_coords, receiver_coords)
  source_receiver_mat$bearing_source_to_receiver <- bearing_matrix
  
  rm(bearing_matrix)
  
  return(source_receiver_mat)
} 


# -----------------------------------------------------------------------------#
# FUNCTION: get_nearest_source_distances_angles 
# PURPOSE:
#   Compute distance and angles for each receiver to the nearest active source
# INPUT:
#   source_dat_sf: sf object with source_id and geometry (point locations)
#   receiver_dat_sf: sf object with receiver_id and geometry (point locations)
#   crs_unprojected_wgs84: CRS for longitude/latitude projection
# OUTPUT:
#   source_receiver_mat: tibble with nearest source, distance from source to receiver,
#     angle from source to receiver - 
#         azimuth angle (North = +y = 0; [0, 360] range)
#         atan angle (North = +y = 0; [-180, 180] range)
#         bearing angle: same as atan, computed with built in R function
# -----------------------------------------------------------------------------#
get_nearest_active_source_distance_bearing <- function(
    source_dat_sf_proj,
    receiver_dat_sf_proj,
    crs_unprojected_wgs84
) {
  # Check CRS
  if (sf::st_crs(receiver_dat_sf_proj)$epsg != sf::st_crs(source_dat_sf_proj)$epsg) {
    stop("Source and Receiver SF data must have matching CRS.")
  }
  
  # Nearest source index for each receiver
  nearest_idx <- sf::st_nearest_feature(receiver_dat_sf_proj, source_dat_sf_proj)
  nearest_sources <- source_dat_sf_proj[nearest_idx, ]
  
  # Compute distance
  distances <- sf::st_distance(receiver_dat_sf_proj, nearest_sources, by_element = TRUE)
  
  # Transform to lon/lat
  receiver_lonlat <- sf::st_transform(receiver_dat_sf_proj, crs = crs_unprojected_wgs84)
  source_lonlat <- sf::st_transform(nearest_sources, crs = crs_unprojected_wgs84)
  
  # Compute azimuth + atan angle
  azimuths <- nngeo::st_azimuth(source_lonlat, receiver_lonlat)
  atan_angles <- sapply(azimuths, adjust_azimuth_to_atan)
  
  # Compute bearing
  receiver_coords <- sf::st_coordinates(receiver_lonlat)
  source_coords <- sf::st_coordinates(source_lonlat)
  bearings <- geosphere::bearing(source_coords, receiver_coords)
  
  # Return result
  tibble::tibble(
    receiver_id                = receiver_dat_sf_proj$receiver_id,
    birth_id                   = receiver_dat_sf_proj$birth_id,
    birthyear                  = receiver_dat_sf_proj$birthyear,
    source_id                  = nearest_sources$source_id,
    distance_m                 = as.numeric(distances),
    distance_km                = distance_m/1000,
    angle_azimuth              = azimuths,
    angle_atan                 = atan_angles,
    bearing_source_to_receiver = bearings
  )
}

# STOP

