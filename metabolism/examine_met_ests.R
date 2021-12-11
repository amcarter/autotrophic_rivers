# Examine snake river metabolism estimates
library(tidyverse)
setwd("C:/Users/Alice Carter/git/autotrophic_rivers/")

# load data ####
snake <- readRDS('data/fit_snake.rds')
# extract DO fits
dat <- snake@data
plot(dat$DO.obs, dat$DO.mod)
par(mfrow = c(3,3), mar = c(3,3,1,1))
for(d in unique(raw$date)){
  # rr <- filter(raw, date == d)
  tt <- filter(dat, date == d)
  plot(tt$solar.time, tt$DO.obs, type = 'l', main = tt$date[1])
  # points(tt$solar.time, tt$DO.obs)
  points(tt$solar.time, tt$DO.mod)
}

raw <- read_csv('data/raw_snake.csv') %>%
  mutate(date = date(solar.time))

