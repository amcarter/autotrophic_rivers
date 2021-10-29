# Pair autotrophic river dataset with Chla an other NWIS available data

# A Carter
# 9/23/2021

# Setup ####
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(dataRetrieval)
library(beepr)
# load in a list of the autotrophic siteyears from streampulse database
dat <- read_csv('data/autotrophic_siteyears_daily_filled_data.csv') %>%
  select(sitecode, year)
dat <- dat[!duplicated(dat),]

# load instantaneous site data:
inst_files <- dir('data/inst', full.names = TRUE)
inst <- read_csv(inst_files[1])
has_vars <- unique(inst$variable)

# See what data is available for these sites ####
pcodes <- select(parameterCdFile, parameter_cd, 
                 parameter_nm, parameter_units) %>%
  as_tibble()
available_data <- data.frame()

for(s in unique(dat$sitecode)){
  site <- substr(s, 6, nchar(s))
  tmp <- whatNWISdata(siteNumber = site, service = c('dv', 'uv')) %>% 
    select(site_no, parameter_cd = parm_cd, stat_cd, begin_date, end_date, 
           count_nu, data_type_cd) %>%
    mutate(site_no = as.numeric(site_no))%>%
    left_join(pcodes, by = 'parameter_cd') %>%
    as_tibble()
  available_data <- bind_rows(available_data, tmp)
  print(s)
}

av_parms <- unique(available_data[available_data$data_type_cd == 'uv',]$parameter_nm)
has_vars
av_parms

codes <- available_data %>% 
  filter(parameter_nm %in% av_parms[c(5, 8, 10, 11, 13, 16:18, 21, 22, 24)], 
         data_type_cd == 'uv') %>%
  select(parameter_cd, parameter_units, parameter_nm)
codes <- codes[!duplicated(codes),]
codes$variable <- c('spc_uScm', 'turbidity_FNU', 'chl_grab_ugl', 
                    'nitrate_nitrite_mgNl', 'fchl_ugl', 'tss_mgl', 
                    'vel_sensor_fs', 'vel_comp_fs', 'nitrate_mgNl',
                    'chl_ugl', 'fdom_ugl')
get_parms_uv <- unique(codes$parameter_cd)
# read in desired datasets: ####
# DO, Temp, Q, spc, chl, turb, 
for(i in 2:length(inst_files)){
  inst <- read_csv(inst_files[i]) %>%
    select(sitecode = siteID, dateTimeUTC, variable, value) 
  nwis_data <- data.frame()
  for(s in unique(inst$sitecode)){
    print(s)
    site <- substr(s, 6, nchar(s))
    error_encountered <- FALSE
    dates <- filter(inst, sitecode == s) %>%
      summarize(start = as.Date(min(dateTimeUTC))-1,
                end = as.Date(max(dateTimeUTC))+1)
    tmp <- tryCatch(readNWISuv(siteNumbers = site, parameterCd = get_parms_uv,
                               startDate = dates$start[1], 
                               endDate = dates$end[1]), 
                  error=function(e) error_encountered <- TRUE)
  if(error_encountered | nrow(tmp) == 0) next
  
  tmp <- tmp %>%
      mutate(sitecode = paste('nwis', site_no, sep = '_'))%>%
      select(-site_no, -ends_with('cd')) %>%
      pivot_longer(starts_with('X'), names_to = 'var', values_to = 'value') %>%
      mutate(parameter_cd = str_match(var, '([0-9]+)_.....$')[, 2]) %>%
      left_join(codes, by = 'parameter_cd') %>%
      select(sitecode, dateTimeUTC = dateTime, variable, value) %>%
      filter(!is.na(value))%>%
      group_by(sitecode, dateTimeUTC, variable) %>%
      summarize(value = mean(value, na.rm = T)) %>%
      ungroup() 
  nwis_data <- bind_rows(nwis_data, tmp)
  }
  inst <- bind_rows(inst, nwis_data)
  write_csv(inst, paste0('data/inst/auto_siteyears_inst_nwis_',i,'.csv'))
}

# Combine instantaneous files into a single, daily file with summarized data
# Note: This needs to be updated so that it summarizes based on local date, not UTC
daily <- data.frame()
for(i in 1:3){
  inst <- read_csv(paste0('data/inst/auto_siteyears_inst_nwis_', i, '.csv'))
  tmp <- inst %>%
    mutate(date = as.Date(dateTimeUTC)) %>%
    group_by(sitecode, date, variable) %>%
    summarize(mean = mean(value, na.rm = T),
              median = median(value, na.rm = T),
              sd = sd(value, na.rm = T)) %>%
    ungroup()
  daily <- bind_rows(daily, tmp)
}

dat <- read_csv('data/autotrophic_siteyears_daily_filled_2.csv')
daily %>%
  ggplot(aes(date, mean, col = sitecode)) +
  geom_point() +
  facet_wrap(.~variable, scales = 'free_y')

daily <- daily %>%
  mutate(mean = case_when(mean < -1000 ~ NA_real_,
                           TRUE ~ mean))
# I am going to group values of nitrate + nitrite with nitrate and 
  # values of chl with fchl for now. Revisit this if that seems problematic
daily <- daily %>%
  mutate(variable = case_when(variable == 'nitrate_nitrite_mgNl' ~ 'nitrate_mgNl',
                              variable == 'fchl_ugl' ~ 'chl_ugl',
                              variable == 'vel_comp_fs' ~ 'velocity_fs',
                              TRUE ~ variable))
write_csv(daily, 'daily_sensor_and_nwis_data_auto.csv') 
