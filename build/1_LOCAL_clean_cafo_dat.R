
#..............................................................................
# Script: 1_LOCAL_clean_cafo_dat.R
# Purpose: Clean and tidy CAFO facility data 
#          - Create animal type indicators
#          - Clean and simplify date fields
#          - Add county information
#
# Data inputs: 
#   - raw CAFO facility data, from CalFF
#
# Data outputs: 
#   - [datainput]_cleaned.csv    (cleaned tabular data)
#   - [datainput]_cleaned_sf.rds (cleaned spatial data)
#
# Dependencies:
#   - R packages: here, dplyr, lubridate, sf
#   - Custom: code/setup_header.R, code/build/county_lookup.R
#..............................................................................

##########################
# Setup
##########################
library(here)
source(file.path(here(),"cal-cafo-births/code/setup_header.R"))
source(file.path(here(), "cal-cafo-births/code/build/county_lookup.R")) 
cafo_dat_name <- "facilities_2025-02-12.csv"

##########################
# Load CAFO data
##########################
cafo_dat <- read.csv(file.path(cafo_data_path, cafo_dat_name)) %>%
  rename(source_id = facility_id)

####################################
# Tidy data
####################################
# Create indicator columns of animal_type
cafo_dat <- cafo_dat %>%
  mutate(
    swine   = ifelse(animal_type == "swine", 1, 0),
    cattle  = ifelse(animal_type == "cattle", 1, 0),
    dairy   = ifelse(animal_type == "dairy", 1, 0),
    poultry = ifelse(animal_type == "poultry", 1, 0),
    other   = ifelse(animal_type %in% c("goats", "sheep", "two or more", "unknown"), 1, 0)
  )


# Clean date formatting
# Format dates appropriately: Convert columns to POSIXct and extract only the year
# (year is the only meaningful part of the date info)
cafo_dat <- cafo_dat %>%
  mutate(
    across(
      c(construction_lower_bound, construction_upper_bound, 
        destruction_lower_bound, destruction_upper_bound),
      ~ year(as.POSIXct(.x, format = "%Y-%m-%d %H:%M:%S"))
    )
  )

# Handle source construction/destruction dates
# Note: all facilities have a construction_upper_bound (never NA)
# if a facility was never destroyed, then destruction_lower_bound is NA
# if destruction_lower_bound is NA, set an arbitrary high value
cafo_dat <- cafo_dat %>%
  mutate(destruction_lower_bound_temp = ifelse(is.na(destruction_lower_bound), 
                                               2050, destruction_lower_bound))

# Add county codes 
cafo_dat$county_name <- cafo_dat$county
cafo_dat <- merge(cafo_dat, county_lookup, 
                  by.x = "county", 
                  by.y = "county_name",
                  all.x = TRUE)

# Convert to spatial object
cafo_dat_sf <- st_as_sf(cafo_dat, 
                        coords = c("longitude", "latitude"), 
                        crs = 4326)

# Save to transfer path for use on server
write.csv(cafo_dat, 
          file.path(transfer_data_path, 
                    paste0(file_path_sans_ext(cafo_dat_name), "_cleaned.csv")), 
          row.names = FALSE)
saveRDS(cafo_dat_sf, 
        file.path(transfer_data_path, 
                  paste0(file_path_sans_ext(cafo_dat_name), "_cleaned_sf.rds")))

# STOP
