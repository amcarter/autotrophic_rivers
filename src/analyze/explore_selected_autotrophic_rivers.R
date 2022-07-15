# Examine site specific trends in autotrophy
# A Carter
# 10/14/2021

# setup ####
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(xts)
library(dygraphs)

# load data from three selected autotrophic rivers 
#   (generated in ar1_model.rmd)
dat <- read_csv('data/selected_autotrophic_rivers_daily.csv')
# sites <- unique(dat$sitecode)
# dat_iv <- data.frame()
# for(i in 1:3){
#   inst <- read_csv(paste0('data/inst/auto_siteyears_inst_nwis_', i, '.csv')) %>%
#     filter(sitecode %in% sites)
#   dat_iv <- bind_rows(dat_iv, inst)
# }
# write_csv(dat_iv, 'data/selected_autotrophic_rivers_inst.csv')
dat_iv <- read_csv('data/selected_autotrophic_rivers_inst.csv')

# ggplot(dat_iv, aes(dateTimeUTC, value, color = sitecode)) +
#   geom_point() +
#   facet_wrap(~variable)

dat %>% select(sitecode, Name) %>%
  distinct()

# East Canyon Creek: ####
# this creek flows from east canyon reservoir near Park City, UT. 
# sattelite imagery shows little riparian vegetation, the usgs site is 
# ~ 20 km downstream of the reservoir. It is monitored due to concerns about DO.
ec_code <- c('nwis_10133650')
ec_daily <- filter(dat, sitecode %in% ec_code)
ec_inst <- filter(dat_iv, sitecode == ec_code[1]) %>%
  pivot_wider(values_from = value, names_from = variable) %>%
  mutate(datetime = lubridate::with_tz(dateTimeUTC, tz = 'MST')) %>%
  select(-dateTimeUTC)

ggplo

ggplot(ec_daily, aes(Date, GPP, col = sitecode)) +
  geom_point()
  
ec_inst %>%
  select(DO_mgL, Depth_m) %>%
  xts(order.by = ec_inst$datetime) %>%
  dygraph() %>%
  dygraphs::dyRangeSelector()

plot(ec_daily$DO.amp, ec_daily$GPP)
plot(ec_daily$Date, ec_daily$GPP + ec_daily$ER)
plot(ec_daily$Date, ec_daily$depth_m, log = 'y')
plot(ec_daily$depth_m, ec_daily$GPP)

# Pecos river:
p_code <- c('nwis_08446500')
p_daily <- filter(dat, sitecode %in% p_code)
p_inst <- filter(dat_iv, sitecode == p_code) %>%
  pivot_wider(values_from = value, names_from = variable) %>%
  mutate(datetime = lubridate::with_tz(dateTimeUTC, tz = 'MST'))# %>%
  # select(-dateTimeUTC)
write_csv(p_inst, 'data/pecos_instantaneous.dat')
ggplot(p_daily, aes(Date, GPP, col = Name)) +
  geom_point() +
  geom_point(aes(y = ER))+
  ylab( 'Metabolism')
  

plot(p_daily$DO.amp, p_daily$GPP)
plot(p_daily$light_PAR, p_daily$GPP)
plot(p_daily$discharge_m3s, p_daily$GPP, log = 'x')

p_daily %>%
  select(GPP, discharge_m3s) %>%
  xts(order.by = p_daily$Date) %>%
  dygraph() %>%
  dyRangeSelector() 

# Grand River
gcode <- 'nwis_04119400'
gdaily <- filter(dat, sitecode == gcode)
plot(gdaily$Date, gdaily$GPP)
plot(gdaily$discharge_m3s, gdaily$GPP, log = 'xy')
plot(gdaily$light_PAR, gdaily$GPP, log = 'y')
plot(gdaily$Date, gdaily$discharge_m3s)

ggplot(dat, aes(Date, GPP)) +
  geom_point() +
  facet_wrap(~Name, scales = 'free')
by(Date) %>%
  summarize(depth = mean(Depth_m, na.rm = T)) %>%
  left_join(ec_daily, by = 'Date') %>%
  select(Date, depth, depth_m)

plot(depth$depth, depth$depth_m)
)
%>%
  group_by(date) %>%
  datorophyll

pcodes <- select(parameterCdFile, parameter_cd, 
                 parameter_nm, parameter_units) %>%
  as_tibble() %>%
  filter(parameter_cd %in% c('00095', '63680','99409','62361')) 
pcodes$name <- c('spc', 'turbidity', 'chl', 'tss')

sites <- substr(unique(dat$sitecode), 6, 13)

uv_raw <- readNWISuv(siteNumbers = sites, parameterCd = unique(pcodes$parameter_cd))
dv_raw <- readNWISdv(siteNumbers = sites, parameterCd = unique(pcodes$parameter_cd))

uv <- uv_raw %>% 
  as_tibble() %>%
  select(site_no, dateTime, ends_with('_00000')) %>%
  pivot_longer(cols = ends_with('00000'), names_to = 'parameter_cd', 
               values_to = 'value') %>%
  mutate(parameter_cd = substr(parameter_cd, nchar(parameter_cd) - 10, 
                               nchar(parameter_cd) - 6)) %>%
  # mutate(parameter_cd = stringr::str_match(parameter_cd, 'X_(.+)?_00000')[, 2],
  #        note = case_when(nchar(parameter_cd) == 5 ~ NA_character_,
  #                         TRUE ~ substr(parameter_cd, 1, 
  #                                       nchar(parameter_cd)-6)),
  #        parameter_cd = case_when(nchar(parameter_cd) == 5 ~ parameter_cd,
  #                                 TRUE ~ substr(parameter_cd, 
  #                                               nchar(parameter_cd)-4,
  #                                               nchar(parameter_cd)))) %>%
  left_join(pcodes, by = 'parameter_cd') %>%
  select(-parameter_nm) %>%
  group_by(site_no, dateTime, parameter_cd, name, parameter_units) %>%
  mutate(value = mean(value, na.rm = T)) %>%
  ungroup()

dv <- dv_raw %>% 
  as_tibble() %>%
  select(site_no, Date, ends_with('_00003')) %>%
  pivot_longer(cols = ends_with('00003'), names_to = 'parameter_cd', 
               values_to = 'value') %>%
  mutate(parameter_cd = stringr::str_match(parameter_cd, 'X_(.+)?_00003')[, 2]) %>%
  left_join(pcodes, by = 'parameter_cd') %>%
  select(-parameter_nm)

dat <- dv %>% 
  mutate(sitecode = paste('nwis', site_no, sep = '_')) %>%
  pivot_wider(id_cols = c('sitecode', 'Date'), names_from = name, values_from = value) %>%
  right_join(dat, by = c('sitecode', 'Date'))


SrwdecckntohRhtdt.

# nest steps:
# - get instantaneous dataframes loaded, plot DO over time to see if met is legit
# - look at chla + gpp relationship where available
# - plot transition point. 
# - read steve's paper. See if you can do that analysis
nwis <- read_csv('data/nwis_data_auto_rivers.csv')