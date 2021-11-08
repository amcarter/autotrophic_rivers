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
         ER = zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP))

write_csv(east, 'data/east_metab.csv')
R <- -east$ER

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

ar1_dat <- list(R = -east$ER, N = length(east$ER))
fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)

saveRDS(fit_ar1, 'src/SAM/stan/fits/east_ar1_fit.rds')

# combined model fit

sam_ar1_dat <- list(R = -east$ER, P = east$GPP, N = length(east$ER),
                    nweight = 5, alpha = rep(1, 5))
fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
saveRDS(fit_comb, 'src/SAM/stan/fits/east_combined_fit.rds')


# SAM fit
sam_dat <- list(R = -east$ER, P = east$GPP, N = length(east$ER),
                nweight = 5, alpha = rep(1, 5))
fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
east$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]
# saveRDS(fit_sam, 'src/SAM/stan/fits/east_sam_fit.rds')

# Autoregressive state space model fit 

ar1_dat <- list(R = -east$ER, N = length(east$ER))
fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
east$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]
# saveRDS(fit_ar1, 'src/SAM/stan/fits/east_ar1_fit.rds')

# combined model fit 

sam_ar1_dat <- list(R = -east$ER, P = east$GPP, N = length(east$ER),
                    nweight = 5, alpha = rep(1, 5))
fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
east$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]
# saveRDS(fit_comb, 'src/SAM/stan/fits/east_combined_fit.rds')

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
         ER <- zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP), !is.na(ER))

# SAM fit
sam_dat <- list(R = -pecos$ER, P = pecos$GPP, N = length(pecos$ER),
                nweight = 5, alpha = rep(1, 5))
fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
pecos$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]
# saveRDS(fit_sam, 'src/SAM/stan/fits/pecos_sam_fit.rds')

# Autoregressive state space model fit 

ar1_dat <- list(R = -pecos$ER, N = length(pecos$ER))
fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
pecos$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]
# saveRDS(fit_ar1, 'src/SAM/stan/fits/pecos_ar1_fit.rds')

# combined model fit 

sam_ar1_dat <- list(R = -pecos$ER, P = pecos$GPP, N = length(pecos$ER),
                    nweight = 5, alpha = rep(1, 5))
fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
pecos$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]
# saveRDS(fit_comb, 'src/SAM/stan/fits/pecos_combined_fit.rds')

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
         ER <- zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP), !is.na(ER))

# SAM fit
sam_dat <- list(R = -grand$ER, P = grand$GPP, N = length(grand$ER),
                    nweight = 5, alpha = rep(1, 5))
fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
grand$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]
# saveRDS(fit_sam, 'src/SAM/stan/fits/grand_sam_fit.rds')

# Autoregressive state space model fit 

ar1_dat <- list(R = -grand$ER, N = length(grand$ER))
fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
grand$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]
# saveRDS(fit_ar1, 'src/SAM/stan/fits/grand_ar1_fit.rds')

# combined model fit 

sam_ar1_dat <- list(R = -grand$ER, P = grand$GPP, N = length(grand$ER),
                    nweight = 5, alpha = rep(1, 5))
fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
grand$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]
# saveRDS(fit_comb, 'src/SAM/stan/fits/grand_combined_fit.rds')

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
         ER <- zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP), !is.na(ER))

# SAM fit
sam_dat <- list(R = -snake$ER, P = snake$GPP, N = length(snake$ER),
                nweight = 5, alpha = rep(1, 5))
fit_sam <- stan(file = 'src/SAM/stan/SAM.stan', data = sam_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
snake$R_hat_sam <- summary(fit_sam, pars = 'R_hat')$summary[,1]
# saveRDS(fit_sam, 'src/SAM/stan/fits/snake_sam_fit.rds')

# Autoregressive state space model fit 

ar1_dat <- list(R = -snake$ER, N = length(snake$ER))
fit_ar1 <- stan(file = 'src/SAM/stan/ar1_model.stan', data = ar1_dat,
                warmup = 500, iter = 1000,
                chains = 4, cores = 4)
snake$R_hat_ar1 <- summary(fit_ar1, pars = 'R_hat')$summary[,1]
# saveRDS(fit_ar1, 'src/SAM/stan/fits/snake_ar1_fit.rds')

# combined model fit 

sam_ar1_dat <- list(R = -snake$ER, P = snake$GPP, N = length(snake$ER),
                    nweight = 5, alpha = rep(1, 5))
fit_comb <- stan(file = 'src/SAM/stan/SAM_ar1.stan', data = sam_ar1_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
snake$R_hat_comb <- summary(fit_comb, pars = 'R_hat')$summary[,1]
# saveRDS(fit_comb, 'src/SAM/stan/fits/snake_combined_fit.rds')

p4<- ggplot(snake, aes(-ER, R_hat_ar1)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'AR1 model')
p5<- ggplot(snake, aes(-ER, R_hat_sam)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'SAM model')
p6<- ggplot(snake, aes(-ER, R_hat_comb)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x='measured ER', y = 'model predicted ER', title = 'Combined model')

p1 <- plot(fit_ar1, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,1)
p2 <- plot(fit_sam, pars = c('a0', 'a1', 'w')) + xlim(-0.2,1)
p3 <- plot(fit_comb, pars = c('a0', 'a1', 'a2', 'w', 'sigma_proc', 'sigma_obs')) + xlim(-0.2,1)

jpeg('figures/snake_model_comp.jpeg', width = 6, height = 6, units = 'in', res = 300)
  ggpubr::ggarrange(p4, p5, p6, p1, p2, p3, nrow = 2, ncol = 3)
dev.off()