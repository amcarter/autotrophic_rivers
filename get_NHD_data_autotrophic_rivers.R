# Script to extract data from the NHD 
# variables of interest: slope, upstream distance, accumulated upslope area
#                       nearest upstream lake or reservoir
# A Carter
# 9/23/2021

# Setup ####
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(nhdplusTools)
library(sf)
# source('C:/Users/Alice Carter/git/oxygen_proc_model/src/batch_summary_nhd_streamcat.R')
# dataframe of daily data from all autotrophic rivers in streampulse database
dat <- read_csv('data/autotrophic_siteyears_daily_filled_data.csv')
sites <- select(dat, -year, -Date, -trophic, -DOY, -GPP_C_filled, -ER_C_filled) %>%
  mutate(site_no = substr(sitecode, 6, nchar(sitecode)))
sites <- sites[!duplicated(sites),]  


# get slopes (don't need to run this again, they are now saved in the dataframe)
# nhd <- as_tibble(get_nhdplus(comid = sites$COMID))
# dat <- nhd %>%
#   select(COMID = comid, lengthkm, slope) %>%
#   right_join(dat, by = 'COMID') %>%
#   mutate(slope_src = 'NHDV2')
# write_csv(dat, 'data/autotrophic_siteyears_daily_filled_data.csv')

# trace upstream flowlines ####
nldi_nwis <- list(featureSource = "comid", 
                  featureID = sites$COMID[1])

flowline <- navigate_nldi(nldi_feature = nldi_nwis, mode = "UT")

subset_file <- tempfile(fileext = ".gpkg")
subset <- subset_nhdplus(comids = flowline$UT$nhdplus_comid,
                         output_file = subset_file,
                         nhdplus_data = "download", 
                         flowline_only = FALSE,
                         return_data = TRUE, overwrite = TRUE)
flowline <- sf::read_sf(subset_file, "NHDFlowline_Network")
catchment <- sf::read_sf(subset_file, "CatchmentSP")
waterbody <- sf::read_sf(subset_file, "NHDWaterbody")
 # I'm pausing here because I haven't been able to get this to work and I am 
# going to change direction for a while. 
# When I come back:
# I'm looking for a way to find 
# 1. total upstream flowpath length
# 2. Distance upstream to waterbodies, are there any tribs in the way, etc.

