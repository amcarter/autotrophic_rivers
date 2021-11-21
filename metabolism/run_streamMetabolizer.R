## Model Metabolism #
# this is still old code, not ready for running on snake river data

# update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE,
#                 repos=c("http://owi.usgs.gov/R", "https://cran.rstudio.com"))
# devtools::install_github("USGS-R/streamMetabolizer", ref="develop")

library(rstan)
library(tidyverse)
library(ggplot2)
library(streamMetabolizer)
library(lubridate)
library(dygraphs)

setwd("C:/Users/Alice Carter/git/autotrophic_rivers/")

## Read in Data ####
sites <- read_csv("siteData/NHCsite_metadata.csv") %>%
  slice(1:5,7)

# select variables for metabolism
read_metdata <- function(site){
  MP <- read_csv(paste0("metabolism/processed/", site, ".csv"), guess_max = 100000) %>%
    select(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge)
  return(MP)
}

NHC <- read_metdata("NHC")
PM <- read_metdata("PM")
CBP <- read_metdata("CBP")
WB <- read_metdata("WB")
WBP <- read_metdata("WBP")
# PWC <- read_metdata("PWC")
UNHC <- read_metdata("UNHC")


# # Visualize the data #####
# dat <- CBP
# 
# dat %>% unitted::v() %>%
#   mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
#   select(solar.time, starts_with('DO')) %>%
#   gather(type, DO.value, starts_with('DO')) %>%
#   mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
#   ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')
# 
# labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)',
#             light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
# dat %>% unitted::v() %>%
#   mutate(discharge = log(discharge)) %>%
#   select(solar.time, depth, temp.water, light, discharge) %>%
#   gather(type, value, depth, temp.water, light, discharge) %>%
#   mutate(
#     type=ordered(type, levels=c('depth','temp.water','light','discharge')),
#     units=ordered(labels[type], unname(labels))) %>%
#   ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')


## Set bayes specs #####
bayes_name <- mm_name(type='bayes', pool_K600="binned", 
                          err_obs_iid=TRUE, err_proc_iid = TRUE, 
                          ode_method = "trapezoid", deficit_src='DO_mod', 
                          engine='stan')


set_up_model <- function(dat, bayes_name, site, one = T){
  
  ## Set bayes specs
  bayes_specs <- specs(bayes_name)
  bayes_specs$keep_mcmc_data <- FALSE
  ## Based on range of log daily Q
  daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>% group_by(date) %>%
    summarize(discharge = mean(discharge, na.rm = T),
              depth = mean(depth, na.rm = T))
  Qrng <- range(log(daily$discharge), na.rm = T)
  delta = 2
  n = 6
  while(delta > 1){
    n = n + 1
    delta <- (Qrng[2]-Qrng[1])/n
  }  
  nodes <- seq(Qrng[1], Qrng[2], length = n)
  ## Based on Pete Raymond's data
  slope <- sites[sites$sitecode==site, ]$slope
  # from bob
  # from Joanna's paper
  if(one == T){
    lnK600 <- 4.77+0.55*log(slope)+(-0.52*(log(median(dat$depth, na.rm = T))))
    # lnK600 <- 6.59 + 0.72 * log(slope) - 0.065* log(median(daily$discharge, na.rm = T))
    bayes_specs$K600_lnQ_nodes_meanlog <- c(rep(lnK600, n))
  } else {
    # Kmeanlog <- 6.59 + 0.72 * log(slope) - 0.065* nodes
    daily$lnK600 = 4.77+0.55*log(slope)+(-0.52*(log(daily$depth)))
    a <- summary(lm(lnK600 ~ log(discharge), daily))$coefficients[,1]
    Kmeanlog <- a[1] + a[2] * nodes
    bayes_specs$K600_lnQ_nodes_meanlog <- Kmeanlog
  }
  bayes_specs$K600_lnQ_nodes_centers <- nodes
  bayes_specs$K600_lnQ_nodes_sdlog <- c(rep(0.7, 7))
  
  return(bayes_specs)
}

dir.create('data/metabolism/modeled/finalQ')
# Model Runs ####
# CBP 
bayes_specs <- set_up_model(CBP, bayes_name, "CBP", one = F)
fit <- metab(bayes_specs, CBP)
saveRDS(fit, "metabolism/modeled/finalQ/fit_cbp_uninformed_raymond_K.rds")
rm(fit)
gc()

# PM 
bayes_specs <- set_up_model(PM, bayes_name, "PM", one = F)
fit <- metab(bayes_specs, PM)
saveRDS(fit, "metabolism/modeled/finalQ/fit_pm_uninformed_raymond_K.rds")
rm(fit)
gc()

# WB 
bayes_specs <- set_up_model(WB, bayes_name, "WB", one = F)
fit <- metab(bayes_specs, WB)
saveRDS(fit, "metabolism/modeled/finalQ/fit_wb_uninformed_raymond_K.rds")
rm(fit)
gc()

# WBP 
bayes_specs <- set_up_model(WBP, bayes_name, "WBP", one = F)
fit <- metab(bayes_specs, WBP)
saveRDS(fit, "metabolism/modeled/finalQ/fit_wbp_uninformed_raymond_K.rds")
rm(fit)
gc()

# PWC 
bayes_specs <- set_up_model(PWC, bayes_name, "PWC", one = F)
fit <- metab(bayes_specs, PWC)
saveRDS(fit, "metabolism/modeled/finalQ/fit_pwc_uninformed_raymond_K.rds")
rm(fit)
gc()

# NHC 
for(i in 2017:2019){
  dd = ymd_hms(paste0(i,"-03-01 04:00:00"))
  n = 365*24*60*60
  if(i == 2019) { n = (365 + 25) *24*60*60}
  dat <- NHC %>% 
    filter(solar.time >= dd,
           solar.time <= dd + n)
  bayes_specs <- set_up_model(dat, bayes_name, "NHC", one = F)
  fit <- metab(bayes_specs, dat)
  saveRDS(fit, paste0("metabolism/modeled/finalQ/fit_nhc_",i,
                      "_uninformed_raymond_K.rds"))
  rm(fit)
  gc()
}
# UNHC 
for(i in 2017:2019){
  dd = ymd_hms(paste0((i),"-03-01 04:00:00"))
  n = 365*24*60*60
  if(i == 2019) { n = (365 + 25) *24*60*60}
  dat <- UNHC %>% 
    filter(solar.time >= dd,
           solar.time <= dd + n)
  bayes_specs <- set_up_model(dat, bayes_name, "UNHC", one = F)
  fit <- metab(bayes_specs, dat)
  saveRDS(fit, paste0("metabolism/modeled/finalQ/fit_unhc_",i,
                      "_uninformed_raymond_K.rds"))
  rm(fit)
  gc()
}

