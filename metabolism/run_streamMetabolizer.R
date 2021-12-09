## Model Metabolism #
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
dat <- read_csv("data/raw_snake.csv") %>%
  select(-regionID, -sitecode, -datetime_UTC)


# # Visualize the data #####
dat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)',
            light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
dat %>% unitted::v() %>%
  mutate(discharge = log(discharge)) %>%
  select(solar.time, depth, temp.water, light, discharge) %>%
  gather(type, value, depth, temp.water, light, discharge) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light','discharge')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')


## Set bayes specs #####
bayes_name <- mm_name(type='bayes', pool_K600="binned", 
                          err_obs_iid=TRUE, err_proc_iid = TRUE, 
                          ode_method = "trapezoid", deficit_src='DO_mod', 
                          engine='stan')
## Set bayes specs
bayes_specs <- specs(bayes_name)

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
bayes_specs$K600_lnQ_nodes_centers <- nodes
  # ## Based on Pete Raymond's data
  # slope <- sites[sites$sitecode==site, ]$slope
  # # from bob
  # # from Joanna's paper
  # if(one == T){
  #   lnK600 <- 4.77+0.55*log(slope)+(-0.52*(log(median(dat$depth, na.rm = T))))
  #   # lnK600 <- 6.59 + 0.72 * log(slope) - 0.065* log(median(daily$discharge, na.rm = T))
  #   bayes_specs$K600_lnQ_nodes_meanlog <- c(rep(lnK600, n))
  # } else {
  #   # Kmeanlog <- 6.59 + 0.72 * log(slope) - 0.065* nodes
  #   daily$lnK600 = 4.77+0.55*log(slope)+(-0.52*(log(daily$depth)))
  #   a <- summary(lm(lnK600 ~ log(discharge), daily))$coefficients[,1]
  #   Kmeanlog <- a[1] + a[2] * nodes
  #   bayes_specs$K600_lnQ_nodes_meanlog <- Kmeanlog
  # }
  # bayes_specs$K600_lnQ_nodes_sdlog <- c(rep(0.7, 7))
  # 

# Model Run ####
fit <- metab(bayes_specs, dat)
saveRDS(fit, "data/fit_snake.rds")
