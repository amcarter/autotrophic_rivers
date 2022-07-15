# Pair daily model output data with input sensor data for auto rivers 
# and with normalized residence time calculations (see loticlentic_synthesis)

# A Carter
# 9/26/2021

# Setup ####
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(lubridate)

# dataframe of daily data from all autotrophic rivers in streampulse database
dat <- read_csv('data/autotrophic_siteyears_daily_filled_data.csv')
in_dat <- read_csv('../loticlentic_synthesis/data/normalized_residence_times_and_kO2.csv',
                   guess_max = 1000000)

all <- in_dat %>%
  select(-GPP, -ER, -regionID, -ends_with('sd')) %>%
  rename(Date = date) %>%
  right_join(dat, by = c('sitecode', 'Date')) %>%
  mutate(NEP_C = GPP_C_filled + ER_C_filled)

write_csv(all, 'data/autotrophic_siteyears_daily_filled_2.csv')

ggplot(all, aes(kO2, iT_norm, col = NEP_C)) +
  geom_point()

ggplot(all, aes(kO2, NEP_C)) + 
  geom_point()
ggplot(all, aes(log(iTR_days), NEP_C)) + 
  geom_point()

# build a dataframe of inStantaneous data ####
sites <- unique(dat$sitecode)
setwd('../loticlentic_synthesis/data/powell_data_import/model_inputs/') 
files <- list.files()
auto <- data.frame()
cc <- c('regionID', 'siteID', 'dateTimeUTC', 'variable', 'value')
for(f in 1:length(files)){
  ff <- files[f]
  if(f == 1){ tmp <- read_csv(ff) } else {
    tmp <- read_csv(ff, col_names = F)
    colnames(tmp) <- cc 
    }
    tmp <- tmp %>%
      mutate(siteID = paste0('nwis_', substr(siteID, 6, nchar(siteID))),
             year = year(dateTimeUTC)) %>%
      filter(siteID %in% sites) 
    d <- filter(dat, sitecode %in% unique(tmp$siteID)) %>%
      mutate(year = year(date))
    tmp <- filter(tmp, year %in% unique(d$year))
  
    auto <- bind_rows(auto, tmp)
    if(nrow(auto)>10000000){
      write_csv(auto, paste0('C:/Users/Alice Carter/git/autotrophic_rivers/data/',
                             'inst/autotrophic_siteyears_instantaneous_',
                              f, '.csv'))
      auto <- data.frame()  
      }
}

setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
files <- dir('data/inst', full.names = TRUE) %>% sort()

files_n <- stringr::str_match(files, '[^0-9]([0-9]+).csv')[, 2] %>%
  as.numeric()
files_ord <- files[order(files_n)]
files_newnames <- sapply(1:length(files),
                         function(x) sub('[0-9]+', x, files_ord[x]))
file.rename(from = files_ord,
            to = files_newnames)

