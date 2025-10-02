##################
# Packages
##################
library(tidyverse) # ggplot, dplyr, tidyr, purrr, readr, tibble, stringr, forcats, lubridate
library(data.table)
library(sf)
library(here)
##################
# File Paths
##################
# LOCAL COMPUTER
root <- here() #

# data folders
transfer_data_path         <- file.path(root, "data")
transfer_output_path       <- file.path(root, "output")
transfer_modelresults_path <- file.path(transfer_output_path, "model_results")
server_cleandata_path      <- file.path(server_root, "data/clean")
server_rawdata_path        <- file.path(server_root, "data/raw")
server_outdata_path        <- file.path(server_root, "data/output")
wind_data_path             <- file.path(root, "raw_cafo_wind_data/WindData")
cafo_data_path             <- file.path(root, "raw_cafo_wind_data/CA_CAFO_Data")

# code paths
build_path    <- file.path(root, "build")
analysis_path <- file.path("analysis")

# output paths
plots_dir  <- file.path(root, "output/plots")
maps_dir   <- file.path(root, "output/maps")
tables_dir <- file.path(root, "output/tables")

######################################################
# Project coordinate reference system (CRS) 
######################################################
# UNPROJECTED CRS: for geographic data (latitude / longitude)
crs_unprojected_nad83 <- st_crs(4269) # units = degrees 
crs_unprojected_wgs84 <- st_crs(4326) # units = degrees 

# PROJECTED CRS: for creating buffers (planar operations)
crs_projected_state <- st_crs(3310) # units = meters


