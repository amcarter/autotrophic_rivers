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
dat <- read_csv("data/raw_snake.csv") 
sitecode <- paste0('nwis_', substr(dat$sitecode[1], 6, 13))
dat <- dat %>%
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
# fit <- metab(bayes_specs, dat)
# saveRDS(fit, "data/fit_snake.rds")
fit <- readRDS("data/fit_snake.rds")
ff <- streamMetabolizer::get_fit(fit)
ff$warnings
ff$errors
ff$daily %>%
  select(date, ends_with('daily_Rhat')) %>%
  pivot_longer(cols = ends_with('hat'), names_to = 'var', values_to = 'rhat') %>%
  ggplot(aes(x = rhat, y = var)) + 
  geom_boxplot()
# The Rhats all look fine.

met <- ff$daily %>%
  select(date, GPP = GPP_mean, GPP_sd, GPP_2.5pct, GPP_97.5pct, 
         GPP_n_eff, GPP_Rhat, ER = ER_mean, ER_sd, ER_2.5pct, ER_97.5pct,
         ER_n_eff, ER_Rhat, K600 = K600_daily_mean, K600_sd = K600_daily_sd,
         K600_2.5pct = K600_daily_2.5pct, K600_97.5pct = K600_daily_97.5pct,
         K600_n_eff = K600_daily_n_eff, K600_Rhat = K600_daily_Rhat, warnings,
         errors)

# Compare new metab estimates to the powell center estimates ####
pw_metab <- read_csv('../loticlentic_synthesis/data/powell_data_import/compiled_daily_model_results.csv') %>%
  filter(site_name == sitecode) 
  
pw_metab <- pw_metab %>%
  rename(GPP_pw = GPP, ER_pw = ER, K600_pw = K600)

comp <- left_join(met, pw_metab, by = 'date') %>%
  mutate(GPP.Rhat = case_when(GPP.Rhat < 1.05 ~ 0,
                              GPP.Rhat >= 1.05 ~ 1,
                              TRUE ~ NA_real_),
         ER.Rhat = case_when(ER.Rhat < 1.05 ~ 0,
                              ER.Rhat >= 1.05 ~ 1,
                              TRUE ~ NA_real_))
jpeg(filename = 'figures/snake_met_fit_comparison.jpeg', width = 7, height = 7,
     res = 300, units = 'in')
  ggplot(comp, aes(date, GPP)) +
    geom_line(col = 'forestgreen') +
    geom_line(aes(y = ER), col = 'black') +
    geom_point(aes(y = GPP_pw, col = GPP.Rhat)) +
    geom_point(aes(y = ER_pw, col = ER.Rhat))
dev.off()

ggplot(comp, aes(GPP, GPP_pw, col = GPP.Rhat)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(comp, aes(ER, ER_pw, col = ER.Rhat)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

# The new metab estimates have considerably better Rhats, but they do not 
# otherwise differ significantly from the powell center estimates

ggplot(comp, aes(log(discharge), (K600))) +
  geom_point() +
  geom_point(aes(y = (K600_pw)), col = 'steelblue')

jpeg(filename = 'figures/snake_K600_fit_comparison.jpeg', width = 5, height = 5,
     res = 300, units = 'in')
  par(mfrow = c(2,2), mar = c(4,4,1,0), oma = c(0,0,0,1))
  plot(comp$K600, comp$ER, main = 'new estimates', xlab = 'K600', ylab = 'ER', 
       xlim = c(0.35, 1.25), ylim = c(-35, 0))
  plot(comp$K600_pw, comp$ER_pw, main = 'powell estimates', xlab = 'K600',
       ylab = '', xlim = c(0.35, 1.25), ylim = c(-35, 0))
  plot(comp$discharge, comp$K600, log = 'xy', ylim = c(0.35, 1.25), ylab = 'K600',
       xlab = 'discharge')
  plot(comp$discharge, comp$K600_pw, log = 'xy', xlab = 'discharge', ylab = '', 
       ylim = c(0.35, 1.25))
dev.off()

met <- pw_metab %>% 
  select(-starts_with(c('GPP', 'ER', 'K600')), -resolution) %>%
  right_join(met, by = 'date')
  
write_csv(met, 'data/met_snake.csv')
