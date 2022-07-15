# Download and format data from nwis
# A Carter
# 11/2021

# Setup ####
library(tidyverse)
library(dataRetrieval)
library(lubridate)

# NWIS site nos to download
site_nos <- c('04067500', '14206950')
pcodes <- select(parameterCdFile, parm_cd = parameter_cd, 
                 parameter_nm, parameter_units) %>%
  as_tibble()

# look at what info is available and select parameter codes:
readNWISsite(site_nos)
av <- whatNWISdata(siteNumber = site_nos, service = 'uv') %>%
  select(parm_cd) %>%
  left_join(pcodes) %>%
  distinct()

#select desired parameters
vars <- av[c(1,2,5),] %>%
  mutate(variable = c('temp_C', 'discharge_cfs', 'DO_mgl'))

# download raw daily data
dat <- readNWISuv(siteNumbers = site_nos,
                  parameterCd = vars$parm_cd, startDate = '2019-01-01', endDate = '2019-01-04')

dat <- dat %>% 
  as_tibble() %>%
  rename(datetime_UTC = dateTime) %>%
  select(-ends_with('cd')) %>%
  pivot_longer(starts_with('X')) %>%
  mutate(parm_cd = substr(name, 3, 7)) %>%
  left_join(vars) %>%
  select(-name, -parm_cd, -parameter_nm, -parameter_units) %>%
  pivot_wider(names_from = 'variable', values_from = 'value')

write_csv(dat, 'raw_nwis_data.csv')
