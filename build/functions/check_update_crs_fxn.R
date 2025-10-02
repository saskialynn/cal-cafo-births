# function to check the crs of an sf object and convert it to the desired crs
check_update_crs <- function(source, # sf object
                              receiver, # sf object
                              projected_crs = crs_projected, # project projected crs
                              geographic_crs = crs_nad83, # project geographic crs
                              output_crs) { # desired crs for both sf objects
  
  # Check if both source and receiver have CRS information
  if (is.null(st_crs(source)) || is.null(st_crs(receiver))) {
    stop("Both source and receiver objects must have defined CRS.")
  }
  
  # Extract EPSG codes
  epsg_source <- sf::st_crs(source)$epsg
  epsg_receiver <- sf::st_crs(receiver)$epsg
  
  # Check if EPSG codes match
  if (epsg_source != epsg_receiver) {
    stop("Source and receiver objects must have the same CRS.")
  }
  
  # Check if EPSG code is projected or geographic
  is_geographic <- st_is_longlat(source) # TRUE: geographic crs ; FALSE: projected crs
  want_geographic <- st_is_longlat(output_crs) # TRUE: geographic crs ; FALSE: projected crs
  
  # Convert CRS if needed
  if (is_geographic && !want_geographic) { # is geographic, want projected
    source <- st_transform(source, output_crs)
    receiver <- st_transform(receiver, output_crs)
    print("converted geographic to projected coordinates")
  } else if (!is_geographic && want_geographic) { # is projected, want geographic
    source <- st_transform(source, output_crs)
    receiver <- st_transform(receiver, output_crs)
    print("converted projected to geographic coordinates")
  }else{ # is geographic, want geographic OR is projected, want projected
    print("no change in crs needed")
  }
  
  # Return the updated source and receiver
  return(list(source = source, receiver = receiver))
}

# sample function use
# temp <- check_update_crs(source = source, 
#                          receiver = receiver, 
#                          output_crs = crs_nad83)
# source <- temp$source
# receiver <- temp$receiver