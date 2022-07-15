# Download and format snake river data from nwis
# A Carter
# 11/2021

# Setup ####
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(dataRetrieval)
library(streamMetabolizer)
library(lubridate)

# There are four sites along the snake river that have oxygen data 
# The first site, at Adrian, OR, is the one that I've been working with so far
site_nos <- c('13173600', '13213100', '13093383', '13013650')
pcodes <- select(parameterCdFile, parm_cd = parameter_cd, 
                 parameter_nm, parameter_units) %>%
  as_tibble()
# get info for Snake River at Adrian, OR:
readNWISsite(site_nos[1])
av <- whatNWISdata(siteNumber = site_nos[1], service = 'uv') %>%
  select(parm_cd) %>%
  left_join(pcodes) %>%
  mutate(variable = c('temp_water_C', 'discharge_m3s', 'spc_uscm', 'DO_mgl',
                      'pH', 'chl_ugl', 'turb_FNU'))

# download raw daily data
dat <- readNWISuv(siteNumbers = site_nos[1],
                  parameterCd = av$parm_cd)

dat <- dat %>% 
  as.tibble() %>%
  rename(datetime_UTC = dateTime) %>%
  select(-ends_with('cd')) %>%
  pivot_longer(starts_with('X')) %>%
  mutate(parm_cd = substr(name, 3, 7)) %>%
  left_join(av) %>%
  select(-name, -parm_cd, -parameter_nm, -parameter_units) %>%
  pivot_wider(names_from = 'variable', values_from = 'value')

plot(dat$datetime_UTC, dat$DO_mgl)
par(new = T)
plot(dat$datetime_UTC, dat$discharge_m3s, log = 'y', type = 'l', col = 4)
par(new = T)
plot(dat$datetime_UTC, dat$chl_ugl, type = 'l', col = 3)
# it looks like only one year of this data was included in the powell center
# because discharge is missing from much of 2009. I will remodel it for the 
# entire period of available data, but I need the depth ests from the pw data

# Pair with powell center data to get depths ####
snake <- read_csv('../loticlentic_synthesis/data/powell_data_import/model_inputs/all_powell_data09.csv')
colnames(snake) <- c('regionID', 'sitecode', 'datetime_UTC', 'variable', 'value')
snake <- snake %>% 
  filter(sitecode == paste0('nwis-', site_nos[1])) %>%
  pivot_wider(values_from = value, names_from = variable)
plot(snake$datetime_UTC, snake$DO_mgL)
plot(snake$datetime_UTC, snake$Discharge_m3s)
plot(snake$Depth_m, snake$Discharge_m3s, log = 'xy')

# Prep snake river data for running stream metabolizer ####

#look at required inputs for a bayesian model in stream Metabolizer
metab_inputs(type="bayes", input="data")
mdat <- readNWISsite(site_nos[1])
lat <- mdat$dec_lat_va
lon <- mdat$dec_long_va

# remove leading and ending NAs
w <- which(!is.na(snake$DO_mgL))
snake <- snake[min(w):max(w), ]
w <- which(!is.na(snake$Discharge_m3s))
snake <- snake[min(w):max(w), ]

#Convert datetime to solar time 
snake$DateTime_MST <- with_tz(snake$datetime_UTC, tz="MST") # convert to PST timezone
snake$solar.time <- calc_solar_time(snake$DateTime_MST, longitude=lon)
snake$light <- calc_light(snake$solar.time, latitude=lat,longitude=lon)
  
snake <- snake %>%
    select(-DateTime_MST, -Light_PAR) %>%
    rename(DO.obs = DO_mgL, 
           DO.sat = satDO_mgL,
           temp.water = WaterTemp_C,
           depth = Depth_m,
           discharge = Discharge_m3s)

  

write_csv(snake, 'data/raw_snake.csv')
