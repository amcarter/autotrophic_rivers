# Test SAM stan models
library(tidyverse)
library(rstan)
library(shinystan)
setwd('C:/Users/Alice Carter/git/autotrophic_rivers')

# run this if the models have been updated
# source('src/SAM/generate_stan_SAM_models.R')

# read in data from autotrophic rivers
dat <- read_csv('data/selected_autotrophic_rivers_daily.csv')


# initial model runs ####
# east canyon creek below I-80. Data from 2011 are questionable, so filtered out
east <- filter(dat, sitecode == 'nwis_10133650', year >= 2012)

# for now, fill in missing data with linear interpolation:
east <- east %>%
  mutate(GPP = zoo::na.approx(GPP, na.rm = F),             
         ER = zoo::na.approx(ER, na.rm = F),
         GPP.sd = (GPP.upper - GPP.lower)/4) %>%
  filter(!is.na(GPP))

# write_csv(east, 'data/east_metab.csv')
R <- -east$ER

# look for trends in error on GPP obs
plot(east$GPP, east$GPP.sd)#, log = 'xy') 
mu_obs = mean(east$GPP.sd, na.rm = T)
tau_obs = sd(east$GPP.sd, na.rm = T)


sam_east_dat <- list(R = -east$ER, P = east$GPP, N = length(east$GPP), 
                    nweight = 5, alpha = rep(1,5))
fit <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_east_dat,
            warmup = 500, iter = 1000,
            chains = 4, cores = 4)
# 
# plot(fit)
# launch_shinystan(fit)
# 
# saveRDS(fit, 'src/SAM/stan/fits/east_sam_fit.rds')

fit <- readRDS('src/SAM/stan/fits/east_sam_fit.rds')
R_hat <- summary(fit, pars = 'R_hat')$summary %>%
  as_tibble() %>%
  pull(mean)

# Do predictions correlate with data?
plot(R, R_hat)
abline(0,1)

# are we accounting for all of the autocorrelation in the respiration?
par(mfrow = c(1,2))
acf(R, main = 'Autocorrelation of ER')
acf(R - R_hat, main = 'Autocorrelation of ER residuals')

# No and it looks like the pattern in the ER residuals might be seasonal
par(mfrow = c(1,1))
plot(east$Date, R-R_hat, xlab = 'date', ylab = 'ER residuals')


# Autoregressive state space model fit

ar1_dat <- list(R = -east$ER, N = length(east$ER), mu_obs = mu_obs, 
                tau_obs = tau_obs)
# fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# 
# saveRDS(fit_ar1, 'src/SAM/stan/fits/east_ar1_fit.rds')
fit_ar1 <- readRDS('src/SAM/stan/fits/east_ar1_fit.rds')
east$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]
pairs(fit_ar1, pars = c('a0', 'a1', 'sigma_obs', 'sigma_proc'))
plot(fit_ar1, pars = c('a0', 'a1', 'sigma_obs', 'sigma_proc'))
print(fit_ar1, pars = c('a0', 'a1', 'sigma_obs', 'sigma_proc'))
# combined model fit

sam_ar1_dat <- list(R = -east$ER, P = east$GPP, N = length(east$ER),
                    nweight = 5, alpha = rep(1, 5), mu_obs = mu_obs,
                    tau_obs = tau_obs)
# fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
#                  warmup = 500, iter = 1000,
#                  chains = 4, cores = 4)
# saveRDS(fit_comb, 'src/SAM/stan/fits/east_combined_fit.rds')
fit_comb <- readRDS('src/SAM/stan/fits/east_combined_fit.rds')
east$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]
plot(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_obs', 'sigma_proc'))
print(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_obs', 'sigma_proc'))


# SAM fit
sam_dat <- list(R = -east$ER, P = east$GPP, N = length(east$ER),
                nweight = 5, alpha = rep(1, 5))
# fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# saveRDS(fit_sam, 'src/SAM/stan/fits/east_sam_fit.rds')
fit_sam <- readRDS('src/SAM/stan/fits/east_sam_fit.rds')
east$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]


p4<- ggplot(east, aes(-ER, R_hat_ar1)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'AR1 model')
p5<- ggplot(east, aes(-ER, R_hat_sam)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'SAM model')
p6<- ggplot(east, aes(-ER, R_hat_comb)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'Combined model')

p1 <- plot(fit_ar1, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,1)
p2 <- plot(fit_sam, pars = c('a0', 'a1', 'w')) + xlim(-0.2,1)
p3 <- plot(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,1)

jpeg('figures/east_model_comp.jpeg', width = 6, height = 6, units = 'in', res = 300)
  ggpubr::ggarrange(p4, p5, p6, p1, p2, p3, nrow = 2, ncol = 3)
dev.off()

# pecos river ####
pecos <- filter(dat, sitecode == 'nwis_08446500') %>%
  mutate(GPP <- zoo::na.approx(GPP, na.rm = F),
         ER <- zoo::na.approx(ER, na.rm = F),
         GPP.sd = (GPP.upper - GPP.lower)/4) %>%
  filter(!is.na(GPP), !is.na(ER))

mu_obs = mean(pecos$GPP.sd, na.rm = T)
tau_obs = sd(pecos$GPP.sd, na.rm = T)
# SAM fit
sam_dat <- list(R = -pecos$ER, P = pecos$GPP, N = length(pecos$ER),
                nweight = 5, alpha = rep(1, 5))
# fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# saveRDS(fit_sam, 'src/SAM/stan/fits/pecos_sam_fit.rds')
fit_sam <- readRDS('src/SAM/stan/fits/pecos_sam_fit.rds')
pecos$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]

# Autoregressive state space model fit 

ar1_dat <- list(R = -pecos$ER, N = length(pecos$ER), mu_obs = mu_obs,
                tau_obs = tau_obs)
# fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# saveRDS(fit_ar1, 'src/SAM/stan/fits/pecos_ar1_fit.rds')
fit_ar1 <- readRDS('src/SAM/stan/fits/pecos_ar1_fit.rds')
pecos$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]

# combined model fit 

sam_ar1_dat <- list(R = -pecos$ER, P = pecos$GPP, N = length(pecos$ER),
                    nweight = 5, alpha = rep(1, 5), mu_obs = mu_obs,
                    tau_obs = tau_obs)
# fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
#                  warmup = 500, iter = 1000,
#                  chains = 4, cores = 4)
# saveRDS(fit_comb, 'src/SAM/stan/fits/pecos_combined_fit.rds')
fit_comb <- readRDS('src/SAM/stan/fits/pecos_combined_fit.rds')
pecos$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]

p4<- ggplot(pecos, aes(-ER, R_hat_ar1)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'AR1 model')
p5<- ggplot(pecos, aes(-ER, R_hat_sam)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'SAM model')
p6<- ggplot(pecos, aes(-ER, R_hat_comb)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'Combined model')

p1 <- plot(fit_ar1, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,6)
p2 <- plot(fit_sam, pars = c('a0', 'a1', 'w')) + xlim(-0.2,6)
p3 <- plot(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,6)

jpeg('figures/pecos_model_comp.jpeg', width = 6, height = 6, units = 'in', res = 300)
  ggpubr::ggarrange(p4, p5, p6, p1, p2, p3, nrow = 2, ncol = 3)
dev.off()


# grand river ####
grand <- filter(dat, sitecode == 'nwis_04119400') %>%
  mutate(GPP <- zoo::na.approx(GPP, na.rm = F),
         ER <- zoo::na.approx(ER, na.rm = F),
         GPP.sd = (GPP.upper - GPP.lower)/4) %>%
  filter(!is.na(GPP), !is.na(ER))

mu_obs = mean(grand$GPP.sd, na.rm = T)
tau_obs = sd(grand$GPP.sd, na.rm = T)

plot(grand$GPP, grand$GPP.sd)
hist(grand$GPP.sd)
# SAM fit
sam_dat <- list(R = -grand$ER, P = grand$GPP, N = length(grand$ER),
                    nweight = 5, alpha = rep(1, 5))
# fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
#                  warmup = 500, iter = 1000,
#                  chains = 4, cores = 4)
# saveRDS(fit_sam, 'src/SAM/stan/fits/grand_sam_fit.rds')
fit_sam <- readRDS('src/SAM/stan/fits/grand_sam_fit.rds')
grand$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]

# Autoregressive state space model fit 

ar1_dat <- list(R = -grand$ER, N = length(grand$ER), tau_obs = tau_obs, 
                mu_obs = mu_obs)
# fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# saveRDS(fit_ar1, 'src/SAM/stan/fits/grand_ar1_fit.rds')
fit_ar1 <- readRDS('src/SAM/stan/fits/grand_ar1_fit.rds')
grand$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]

# combined model fit 

sam_ar1_dat <- list(R = -grand$ER, P = grand$GPP, N = length(grand$ER),
                    nweight = 5, alpha = rep(1, 5), mu_obs = mu_obs,
                    tau_obs = tau_obs)
# fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
#                  warmup = 500, iter = 1000,
#                  chains = 4, cores = 4)
# saveRDS(fit_comb, 'src/SAM/stan/fits/grand_combined_fit.rds')
fit_comb <- readRDS('src/SAM/stan/fits/grand_combined_fit.rds')
grand$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]

p4<- ggplot(grand, aes(-ER, R_hat_ar1)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'AR1 model')
p5<- ggplot(grand, aes(-ER, R_hat_sam)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'SAM model')
p6<- ggplot(grand, aes(-ER, R_hat_comb)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'Combined model')

p1 <- plot(fit_ar1, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,2.4)
p2 <- plot(fit_sam, pars = c('a0', 'a1', 'w')) + xlim(-0.2,2.4)
p3 <- plot(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,2.4)

jpeg('figures/grand_model_comp.jpeg', width = 6, height = 6, units = 'in', res = 300)
  ggpubr::ggarrange(p4, p5, p6, p1, p2, p3, nrow = 2, ncol = 3)
dev.off()

# snake river ####
snake <- filter(dat, sitecode == 'nwis_13173600') %>%
  mutate(GPP <- zoo::na.approx(GPP, na.rm = F),
         ER <- zoo::na.approx(ER, na.rm = F),
         GPP.sd = (GPP.upper - GPP.lower)/4) %>%
  filter(!is.na(GPP), !is.na(ER))

plot(snake$GPP, snake$GPP.sd)
mu_obs = mean(snake$GPP.sd, na.rm = T)
tau_obs = sd(snake$GPP.sd, na.rm = T)

# SAM fit
sam_dat <- list(R = -snake$ER, P = snake$GPP, N = length(snake$ER),
                nweight = 5, alpha = rep(1, 5))
# fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# saveRDS(fit_sam, 'src/SAM/stan/fits/snake_sam_fit.rds')
fit_sam <- readRDS( 'src/SAM/stan/fits/snake_sam_fit.rds')
snake$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]

# Autoregressive state space model fit 

ar1_dat <- list(R = -snake$ER, N = length(snake$ER), mu_obs = mu_obs, 
                tau_obs = tau_obs)
# fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
#                 warmup = 500, iter = 1000,
#                 chains = 4, cores = 4)
# saveRDS(fit_ar1, 'src/SAM/stan/fits/snake_ar1_fit.rds')
fit_ar1 <- readRDS('src/SAM/stan/fits/snake_ar1_fit.rds')
snake$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]

# combined model fit 
sam_ar1_dat <- list(R = -snake$ER, P = snake$GPP, N = length(snake$ER),
                    nweight = 5, alpha = rep(1, 5), mu_obs = mu_obs, 
                    tau_obs = tau_obs)
# fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
#                  warmup = 500, iter = 1000,
#                  chains = 4, cores = 4)
# saveRDS(fit_comb, 'src/SAM/stan/fits/snake_combined_fit.rds')
fit_comb <- readRDS('src/SAM/stan/fits/snake_combined_fit.rds')
snake$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]


jpeg('figures/snake_residual_acfs.jpeg', width = 6, height = 6, 
     units = 'in', res = 300)
  par(mfrow = c(2,2), mar = c(2,2,1,1))
  acf(snake$ER)
  mtext('ER', 3, -2)
  acf(snake$ER - snake$R_hat_ar1)
  mtext('ar1 residuals', 3, -2)
  acf(snake$ER - snake$R_hat_sam)
  mtext('SAM residuals', 3, -2)
  acf(snake$ER - snake$R_hat_comb)
  mtext('combined residuals', 3, -2)
dev.off()

p4<- ggplot(snake, aes(-ER, R_hat_ar1)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'AR1 model')
p5<- ggplot(snake, aes(-ER, R_hat_sam)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'SAM model')
p6<- ggplot(snake, aes(-ER, R_hat_comb)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'Combined model')

p1 <- plot(fit_ar1, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,3.2)
p2 <- plot(fit_sam, pars = c('a0', 'a1', 'w')) + xlim(-0.2,3.2)
p3 <- plot(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,3.2)

jpeg('figures/snake_model_comp.jpeg', width = 6, height = 6, units = 'in', res = 300)
  ggpubr::ggarrange(p4, p5, p6, p1, p2, p3, nrow = 2, ncol = 3)
dev.off()


jpeg('figures/sd_variance.jpeg', width = 6, height = 6, units = 'in', res = 300)
  par(mfrow = c(2,2), mar = c(2,2,2,1), oma = c(3,3,0,0))
  plot(east$GPP, east$GPP.sd, xlab = "", ylab = "", main = "East CC", log = 'xy')  
  plot(pecos$GPP, pecos$GPP.sd, xlab = "", ylab = "", main = "Pecos", log = 'xy')  
  plot(grand$GPP, grand$GPP.sd, xlab = "", ylab = "", main = "Grand", log = 'xy')  
  plot(snake$GPP, snake$GPP.sd/snake$GPP, xlab = '', ylab = '', main = "Snake")  
  par(new = T, mfrow = c(1,1))
  mtext('GPP (gO2/m2/d)', 1, 3.2)  
  mtext('GPP stdev (gO2/m2/d)', 2, 3.2)
dev.off()
  