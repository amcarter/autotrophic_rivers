# Trial spectral analysis on the five percent of most autotrophic and most
# heterotrophic rivers in the Powell center database 
# See extract_autotrophic_rivers.R for how this dataframe was generated 

# A Carter
# 2021-10

setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(lubridate)

# read in compiled auto and hetero streams:
dat <- read_csv('data/autohetero_quantiles_siteyears_daily_filled_data.csv')

#group into autotrophic and heterotrophic:
auto <- filter(dat, trophic == 'autotrophic')
hetero <- filter(dat, trophic == 'heterotrophic')

# loop through site years to calculate spectral densities
# extract __
for(s in unique(auto$sitecode)){
  site <- filter(auto, sitecode == s)
  for(y in unique(site$year)){
    sy <- filter(site, year == y)
    spectrum(site$GPP_C_filled)
    x.spec <- spectrum(sy$GPP_C_filled)
    spx <- x.spec$freq
    spy <- 2*x.spec$spec
    plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")
    