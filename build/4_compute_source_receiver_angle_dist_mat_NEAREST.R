
#...............................................................................
# Script: 4_compute_source_receiver_angle_dist_mat_NEAREST.R
#
# Purpose: 
#   Compute source/receiver angles and distances.
#   - Identify nearest active CAFO (during birth year) for each birth
#   - Compute distance and bearing between source and receiver
#
# Data Inputs: 
#   - Births dataset (cleaned)
#   - CAFO dataset (cleaned)
#
# Data Outputs: 
#   - "..._NEARESTCAFO_DISTANGLE.csv": Births data with nearest CAFO info added
#
# Helper Functions: 
#   - check_update_crs()
#   - get_nearest_active_source_distance_bearing()
#
# Dependencies: 
#   - R packages: here, dplyr, sf, data.table, future, furrr, tools
#   - Custom scripts: 
#       - code/setup_header.R
#       - code/functions/get_pairwise_distance_angles_fxn.R
#       - code/functions/check_update_crs_fxn.R
#...............................................................................


##########################
# Setup
##########################
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
source(file.path(build_path, "functions/check_update_crs_fxn.R"))
source(file.path(build_path, "functions/get_pairwise_distance_angles_fxn.R"))

# parallel processing 
library(future)
library(furrr)
library(tools)

#######################################
# Specify Parameters & Input Data
#######################################
birth_dat_name <- "births_main_analysis_df.csv" 
cafo_dat_name  <- "facilities20250212_cleaned_sf.rds"
cafo_dat_id    <- "20250212"

# Maximum buffer size (meters)
max_buffer <- 250000 

######################################
# Link Source and Receiver Locations 
######################################
# Receiver: Maternal residential addresses from CA births data
receiver <- data.table::fread(file.path(server_cleandata_path,  birth_dat_name))
length(unique(paste0(receiver$receiver_id, receiver$birth_id))) 

receiver_sf <- st_as_sf(receiver, 
                        coords = c("long", "lat"), 
                        crs = crs_unprojected_wgs84, 
                        remove = FALSE) 
receiver_sf_lonlat <- receiver_sf

# Source: CAFO locations
source_sf <- readRDS(file.path(transfer_data_path, cafo_dat_name))
source_sf_lonlat <- source_sf

# Ensure matching (unprojected) CRS
stopifnot(st_crs(source_sf) == st_crs(receiver_sf))
stopifnot(st_crs(source_sf)$epsg == st_crs(receiver_sf)$epsg)

# Convert to projected CRS for buffering operations
temp <- check_update_crs(source         = source_sf,
                         receiver       = receiver_sf,
                         projected_crs  = crs_projected_state,
                         geographic_crs = crs_unprojected_wgs84,
                         output_crs     = crs_projected_state)

source_sf <- temp$source
receiver_sf <- temp$receiver
st_is_longlat(source_sf)
st_is_longlat(receiver_sf)

#########################################
# NEAREST SOURCE TO RECEIVER ANALYSIS
#########################################
# Split receivers by birth year
receiver_year_chunks <- receiver_sf %>%
  group_split(birthyear, .keep = TRUE)
names(receiver_year_chunks) <- sapply(receiver_year_chunks, 
                                      function(chunk) unique(chunk$birthyear))

# Further split into smaller chunks (while keeping birth years together) for parallel processing
num_cores <- 20
receiver_chunks <- list()

for(y in 1:length(receiver_year_chunks)){
  receiver_chunks[[y]] <- receiver_year_chunks[[y]] %>%
    dplyr::select(receiver_id, birth_id, birthyear, geometry) %>%  
    mutate(chunk_id = cut(
      row_number(),
      breaks = num_cores/length(receiver_year_chunks),
      labels = FALSE,
      include.lowest = TRUE
    )) %>%
    group_split(chunk_id, .keep = TRUE) %>%
    setNames(paste0("year_", unique(receiver_year_chunks[[y]]$birthyear), 
                    "_chunk_", seq_len(length(.))))
}
names(receiver_chunks) <- names(receiver_year_chunks)
all_chunks <- flatten(receiver_chunks)

# Prepare list of active CAFOs for each birth year
source_by_birthyear <- list()
for (year in min(receiver$birthyear):max(receiver$birthyear)) {
  active_sources <- source_sf %>%
    # constructed before birth year, destroyed after birth year = present
    filter(year > construction_upper_bound, year < destruction_lower_bound_temp) %>% 
    dplyr:: select(source_id, construction_upper_bound, destruction_lower_bound_temp, geometry)
  
  source_by_birthyear[[as.character(year)]] <- active_sources
}

# Pair receiver chunks with valid source data
paired_chunks <- list()
for (chunk in all_chunks) {
  year <- unique(chunk$birthyear)
  source_chunk <- source_by_birthyear[[as.character(year)]]
  paired_chunks[[length(paired_chunks) + 1]] <- list(
    receiver = chunk,
    source = source_chunk
  )
}

# Get nearest active source (CAFO) for each receiver (birth)
plan(multisession, workers = num_cores) 

system.time(
  nearest_results <- furrr::future_map_dfr(
    paired_chunks,
    function(chunk_pair) {
      get_nearest_active_source_distance_bearing(
        source_dat_sf_proj = chunk_pair$source,
        receiver_dat_sf_proj = chunk_pair$receiver,
        crs_unprojected_wgs84 = crs_unprojected_wgs84
      )
    },
    .options = furrr::furrr_options(seed = 123)
  )
)
plan(sequential) 

# Add source and receiver information to results
nearest_results_save <- nearest_results %>%
  left_join(source_sf_lonlat, by = "source_id") %>% 
  dplyr::select(-geometry) %>%
  left_join(receiver_sf_lonlat, by = c("receiver_id", "birth_id", "birthyear")) %>% 
  dplyr::select(-geometry)

save_name <- paste0(tools::file_path_sans_ext(birth_dat_name),  
                    "_cafo", cafo_dat_id, 
                    "_NEARESTCAFO_DISTANGLE.csv")

data.table::fwrite(nearest_results_save, file.path(server_cleandata_path, save_name))

# STOP




