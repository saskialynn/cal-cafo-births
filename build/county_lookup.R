#..............................................................................
# Script: county_lookup.R
# Purpose: Provides a lookup table linking California county names to numeric codes
#
# Data inputs: None (hardcoded)
# Data outputs: `county_lookup` dataframe
#
# Dependencies: None
#..............................................................................

county_lookup <- data.frame(
  county_code = 1:58,
  county_name = c("Alameda", "Alpine", "Amador", 
                  "Butte", 
                  "Calaveras", "Colusa", 
                  "Contra Costa", 
                  "Del Norte", 
                  "El Dorado", 
                  "Fresno", 
                  "Glenn", 
                  "Humboldt", 
                  "Imperial", "Inyo", 
                  "Kern", "Kings", 
                  "Lake", "Lassen", "Los Angeles", 
                  "Madera", "Marin", "Mariposa", "Mendocino", "Merced", "Modoc", "Mono", "Monterey", 
                  "Napa", "Nevada", 
                  "Orange", 
                  "Placer", "Plumas", 
                  "Riverside", 
                  "Sacramento", "San Benito", "San Bernardino", "San Diego", "San Francisco", 
                  "San Joaquin", "San Luis Obispo", "San Mateo", "Santa Barbara", "Santa Clara", 
                  "Santa Cruz", "Shasta", "Sierra", "Siskiyou", "Solano", "Sonoma", "Stanislaus", "Sutter", 
                  "Tehama", "Trinity", "Tulare", "Tuolumne", 
                  "Ventura", 
                  "Yolo", 
                  "Yuba")
)