# Test SAM stan models
library(tidyverse)
library(rstan)
setwd('C:/Users/Alice Carter/git/autotrophic_rivers')
source('src/SAM/generate_stan_SAM_models.R')

# read in data from autotrophic rivers
dat <- read_csv('data/selected_autotrophic_rivers_daily.csv')

# east canyon creek below I-80. Data from 2011 are questionable, so filtered out
east <- filter(dat, sitecode == 'nwis_10133650', year >= 2012)

# for now, fill in missing data with linear interpolation:
east <- east %>%
  mutate(GPP = zoo::na.approx(GPP, na.rm = F),             
         ER = zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP))

sam_east_dat <- list(R = -east$ER, P = east$GPP, N = length(east$GPP), 
                    nweight = 5, alpha = rep(1,5))
fit <- stan(file = 'src/SAM/stan/SAM_pi_5.stan', data = sam_east_dat, 
            warmup = 500, iter = 1000, 
            chains = 4, cores = 4)


fit2 <- stan(file = 'src/SAM/stan/SAM_oipi_5.stan', data = sam_east_dat,
             warmup = 500, iter = 1000,
             chains = 4, cores = 4)


setwd('C:/Users/Alice Carter/git/autotrophic_rivers')