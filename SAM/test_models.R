# Test SAM stan models with simulated data
library(tidyverse)
library(beepr)
library(rstan)
library(shinystan)
setwd('C:/Users/alice.carter/git/autotrophic_rivers')
source('src/SAM/simulate/simulation_functions.R')
source('src/stan_helpers.R')
rstan_options(javascript = FALSE)
# Model Versions ####
# Model elements:
# 1.detC = Latent Detrital Carbon
#   logpi = lognormal process error on latent carbon state
#   logC = log Carbon with normal process error
#   ss = shear stress based loss term (segatto et al 2021)
# 2.SAM = Stochastic Antecedent Model based on GPP
#   N = number of antecedent intervals included
#   logint = SAM calculated using log distributed antecedent intervals
#
# Other things that might get added:
# - log Carbon with normally distributed process error
# - a different model for carbon or detrital respiration
# - fit C0 rather than fixed
# Data from powell center rivers ####
dat <- read_csv('data/data_working/high_quality_daily_met.csv')

# Format sites for simulation

auto_sites <- filter(dat, trophic_stat == 'auto') %>%
    select(long_name, site_name, GPP) %>%
    group_by(long_name, site_name) %>%
    summarize(N = n(), GPP = mean(GPP, na.rm = T)) %>%
    ungroup() %>% filter(N>600) %>% arrange(-GPP, -N)
het_sites <- filter(dat, trophic_stat == 'hetero', width <50) %>%
    group_by(long_name, site_name) %>%
    summarize(N = n(), NEP = mean(GPP+ER, na.rm = T),
              width = mean(width, na.rm = T), LAI = mean(LAI, na.rm = T)) %>%
    select(long_name, site_name, N, NEP, width, LAI) %>%
    ungroup() %>% filter(N>400, LAI >2) %>%
    arrange(width,  -N)

# dat %>% filter(long_name %in% unique(auto_sites$long_name)) %>%
#     ggplot(aes(DOY, ER, col = log(width))) +
#     geom_point() + geom_point(aes(y = GPP)) +
#     facet_wrap(.~long_name, scales = 'free')
# dat %>% filter(long_name %in% unique(het_sites$long_name)) %>%
#     ggplot(aes(DOY, ER, col = LAI)) + geom_point() +
#     facet_wrap(.~long_name, scales = 'free')
# pick some sites
# Heterotrophic: Fanno Creek,OR; New River, TN
# Autotrophic: Grand River, MI; East Canyon Creek,Jeremy; Clackamas; Pecos, sheffield


# simulate data
prep_data <- function(dat){
    dd <- dat %>%
        arrange(date) %>%
        select(sitecode = site_name, long_name, date, DOY, GPP, ER,
               depth, discharge, Stream_PAR, LAI, temp_C = temp.water) %>%
        mutate(across(c(-sitecode, -long_name, -date, -DOY),
                      zoo::na.approx, na.rm = F),
            light = Stream_PAR/max(Stream_PAR, na.rm = T),
            Q = discharge/max(discharge, na.rm = T)) %>%
        filter(!is.na(ER))

    dd$litter <- calc_litter_from_LAI(dd)

    return(dd)
}

newR <- filter(dat, grepl('^NEW RIVER', long_name)) %>%
    prep_data() %>% filter(date < as.Date('2016-01-01'))

fanno <- filter(dat, grepl('^FANNO', long_name)) %>%
    prep_data()

grand <- filter(dat, grepl('^GRAND', long_name)) %>%
    prep_data() %>% filter(date >= as.Date('2016-01-01')) %>%slice(1:731)

east <- filter(dat, grepl('^EAST.*?RANCH', long_name)) %>%
    prep_data()

clack <- filter(dat, grepl('^CLACKA', long_name)) %>%
    prep_data() %>% filter(date >= as.Date('2010-01-01') &
                               date < as.Date('2014-01-01'))

pecos <- filter(dat, grepl('^Pecos.*?Sheff', long_name)) %>%
    prep_data() %>% filter(date >= as.Date('2014-01-01')&
                               date < as.Date('2016-01-01'))

# bind_rows(newR, fanno, grand, east, clack, pecos) %>%
#     ggplot(aes(date, GPP, col = sitecode)) +
#     geom_line() + geom_line(aes(y = ER)) +
#     facet_wrap(.~long_name, ncol = 2, scales = 'free')
# bind_rows(newR, fanno, grand, east, clack, pecos) %>%
#     ggplot(aes(date, LAI, col = sitecode)) +
#     geom_line() + geom_line(aes(y = litter)) +
#     facet_wrap(.~long_name, ncol = 2, scales = 'free')


# Detrital Carbon Models: ####
# Run simulations on the New River and Clackamas:
# Constants
E_a = 0.63            # activation energy for heterotrophic respiration
k_b = 8.6173 * 10^-5  # Boltzmann's constant in eV/K
ARf = 0.44            # the fraction of GPP respired by autotrophs

# Model: detC_logpi_1 ####
# define parameters
C0 = 100       # Initial organic C
K_20 = .01        # Heterotrophic respiration at 20 C
beta_s = 0.8   # Percent removal of organic carbon from max storm
sigma_proc = .02 # lognormally distributed, so this is a % error
sigma_obs = 0.08

simulate_detC_logpi <- function(ss){
    ss <- ss %>%
        mutate(C = 0.5,
               R = 0.5,
               R_obs = ER) %>%
        select(date, GPP, ER, R_obs, R, light, Q, temp_C, litter, C)

    ndays <- nrow(ss) # number of days
    ss$K = calc_rate_coef(ss$temp_C, K_20 = K_20)

    Chat = numeric()
    Chat[1] = C0
    ss$C[1] = exp(rnorm(1, log(C0), sigma_proc))
    ss$AR = -ARf * ss$GPP
    ss$R[1] = ss$AR[1] - ss$K[1] * ss$C[1]

    for(i in 2:ndays){
        Chat[i] = (ss$C[i-1] + ss$litter[i] + ss$R[i-1] + ss$GPP[i-1])*
            (1-beta_s*ss$Q[i])
        ss$C[i] = exp(rnorm(1, log(Chat[i]), sigma_proc))
        ss$R[i] = ss$AR[i] - ss$K[i]*ss$C[i]
    }
    # observation model:
    ss$R_obs = rnorm(ndays, ss$R, sigma_obs)

    return(ss)
}

# Simulate and run
sim_new <- simulate_detC_logpi(newR)
sim_new %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_new), R_obs = sim_new$R_obs, P = sim_new$GPP,
#                  C0 = 100, temp = sim_new$temp_C, litter = sim_new$litter,
#                  Q = sim_new$Q)
# mod <- stan('src/SAM/stan/detC_logpi_1.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             warmup = 500, iter = 1000)
# beep(sound = 8)
# saveRDS(mod, 'src/SAM/stan/fits/detC_logpi_1_new_sim.rds')

sim_clack <- simulate_detC_logpi(clack)
sim_clack %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_clack), R_obs = sim_clack$R_obs, P = sim_clack$GPP,
#                  C0 = 100, temp = sim_clack$temp_C, litter = sim_clack$litter,
#                  Q = sim_clack$Q)
# mod <- stan('src/SAM/stan/detC_logpi_1.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             warmup = 500, iter = 1000)
# beep(sound = 8)
# saveRDS(mod, 'src/SAM/stan/fits/detC_logpi_1_clack_sim.rds')

# evaluate model fit
mod <- readRDS('src/SAM/stan/fits/detC_logpi_1_new_sim.rds')
pars = c('beta_s', 'sigma_obs', 'sigma_proc', 'K_20')
print(mod, pars)
t <- traceplot(mod, ncol = 2, pars)
p <- plot_post_sim(mod, pars, vals = c( beta_s, sigma_obs, sigma_proc*10, K_20*100),
              xlim = c(0,1.5))
png('figures/simulation_fits/detC_logpi_sim_newriver.png', height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery det C (New River)',
                                                          size = 14))
dev.off()

pp <- calc_pp_ests(mod, newR, pars= pars, factors = c(1, 1, .1, .01),
                   sim_func = simulate_detC_logpi)
png('figures/simulation_fits/PPcheck_detC_logpi_newriver_sim.png')
    plot_pp_interval(sim_new$R_obs, pp)#, ylim = c(-30,0))
dev.off()


mod <- readRDS('src/SAM/stan/fits/detC_logpi_1_clack_sim.rds')
pars = c('beta_s', 'sigma_obs', 'sigma_proc', 'K_20')
print(mod, pars)
t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'sigma_obs', 'sigma_proc', 'K_20'))
p <- plot_post_sim(mod, pars = c('beta_s', 'sigma_obs', 'sigma_proc', 'K_20'),
              vals = c( beta_s, sigma_obs, sigma_proc*10, K_20*100), xlim = c(0, 1.2))
png('figures/simulation_fits/detC_logpi_sim_clackamas.png', height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery det C (New River)',
                                                          size = 14))
dev.off()

pp <- calc_pp_ests(mod, clack, pars= pars, factors = c(1, 1, .1, .01),
                   sim_func = simulate_detC_logpi)
png('figures/simulation_fits/PPcheck_detC_logpi_clack_sim.png')
    plot_pp_interval(sim_clack$R_obs, pp, xrng = c(1,365))
dev.off()


# Stochastic Antecedent Models: ####
# Model: SAM5 ####
# define parameters
beta_0 = 0.3   # antecedent P coefficient
beta_p = 0.4   # antecedent P coefficient
sigma_obs = 0.08
nweights <- 5
w <- c(0,0,.5,.5,0)
antdays = 5

# Constants
ARf = 0.44            # the fraction of GPP respired by autotrophs

simulate_SAM5 <- function(ss, ANT, antdays){
    ss <- ss %>%
        mutate(R = ER,
               R_obs = ER)

    ndays <- nrow(ss) # number of days
    ss$Pant = ANT %*% w
    ss$AR = -ARf * ss$GPP

    for(i in (antdays+1):ndays){
        ss$R[i] = beta_0 - beta_p * ss$Pant[i] + ss$AR[i]
    }

    ss$R_obs = rnorm(ndays, ss$R, sigma_obs)

    return(ss)
}

# Simulate and run
ANT <- as.matrix(calc_antecedent_drivers(newR$GPP, 5, 1))
sim_new <- simulate_SAM5(newR, ANT, antdays)

# stan_dat <- list(ndays = nrow(sim_new), nweights = 5, antdays = 5,
#                  R_obs = sim_new$R_obs, P = sim_new$GPP, ANT = ANT)
# mod <- stan('src/SAM/stan/SAM5.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             warmup = 500, iter = 1000)
# beep(sound = 8)
# saveRDS(mod, 'src/SAM/stan/fits/SAM5_new_sim.rds')

ANT <- as.matrix(calc_antecedent_drivers(clack$GPP, 5, 1))
sim_clack <- simulate_SAM5(clack, ANT, antdays )

# stan_dat <- list(ndays = nrow(sim_clack), nweights = 5, antdays = 5,
#                  R_obs = sim_clack$R_obs, P = sim_clack$GPP, ANT = ANT)
# mod <- stan('src/SAM/stan/SAM5.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             warmup = 500, iter = 1000)
# beep(sound = 8)
# saveRDS(mod, 'src/SAM/stan/fits/SAM5_clack_sim.rds')

# evaluate model fit
mod <- readRDS('src/SAM/stan/fits/SAM5_new_sim.rds')
pars = c('beta_0', 'beta_p', 'sigma_obs', 'w')
print(mod, pars)
t <- traceplot(mod, ncol = 2, pars)
p <- plot_post_sim(mod, pars, vals = c(beta_0, beta_p, sigma_obs, w),
                   xlim = c(0,1.5))
png('figures/simulation_fits/SAM5_sim_newriver.png', height = 300, width = 650)
plot <- ggpubr::ggarrange(p, t)
ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery det C (New River)',
                                                      size = 14))
dev.off()

pp <- calc_pp_ests(mod, newR, pars= pars,
                   sim_func = simulate_SAM5)
png('figures/simulation_fits/PPcheck_detC_logpi_newriver_sim.png')
plot_pp_interval(sim_new$R_obs, pp)#, ylim = c(-30,0))
dev.off()
# Model: SAM5_detC_logpi_1 ####
# define parameters
C0 = 100       # Initial organic C
K_20 = .01        # Heterotrophic respiration at 20 C
beta_s = 0.8   # Percent removal of organic carbon from max storm
beta_p = 0.3   # antecedent P coefficient
sigma_proc = .02 # lognormally distributed, so this is a % error
sigma_obs = 0.08
nweights <- 5
w <- c(0,0,.5,.5,0)

# Constants
E_a = 0.63            # activation energy for heterotrophic respiration
k_b = 8.6173 * 10^-5  # Boltzmann's constant in eV/K
ARf = 0.44            # the fraction of GPP respired by autotrophs

simulate_detC_SAM5 <- function(ss){
    ss <- ss %>%
        mutate(C = 0.5,
               R = 0.5,
               HR = 0.5, HR_alg = .5,
               R_obs = ER)

    ndays <- nrow(ss) # number of days
    ss$K = calc_rate_coef(ss$temp_C, K_20 = K_20)

    ss$Pant = ss$GPP

    for (i in (nweights+1):ndays){
        Pvec <- numeric(nweights)
        for(j in 1:nweights){
            Pvec[j] <- w[j]*ss$GPP[i-j]
        }
        ss$Pant[i]<-sum(Pvec)
    }


    Chat = numeric()
    Chat[1] = C0
    ss$C[1] = exp(rnorm(1, log(C0), sigma_proc))
    ss$HR[1] = - ss$K[1] * ss$C[1]

    for(i in 2:ndays){
        Chat[i] = (ss$C[i-1] + ss$litter[i] + ss$HR[i-1])*
            (1-beta_s*ss$Q[i])
        ss$C[i] = exp(rnorm(1, log(Chat[i]), sigma_proc))
        ss$HR[i] = -ss$K[i]*ss$C[i]
    }

    ss$AR = -ARf * ss$GPP
    ss$HR_alg = - beta_p * ss$Pant
    ss$R = ss$AR + ss$HR + ss$HR_alg

    # observation model:
    ss$R_obs = rnorm(ndays, ss$R, sigma_obs)

    return(ss)
}

# Simulate and run
sim_new <- simulate_detC_SAM5(newR)
sim_new %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_new), nweights = 5, R_obs = sim_new$R_obs,
#                  P = sim_new$GPP, C0 = 100, temp = sim_new$temp_C,
#                  litter = sim_new$litter, Q = sim_new$Q)
# mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             warmup = 500, iter = 1000)
# beep(sound = 8)
# saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_new_sim.rds')

sim_clack <- simulate_detC_SAM5(clack)
sim_clack %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_clack), nweights = 5, R_obs = sim_clack$R_obs,
#                  P = sim_clack$GPP, C0 = 100, temp = sim_clack$temp_C,
#                  litter = sim_clack$litter, Q = sim_clack$Q)
# mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             warmup = 500, iter = 1000)
# beep(sound = 8)
# saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_clack_sim.rds')

# look at model fit
mod <- readRDS('src/SAM/stan/fits/SAM5_detC_logpi_new_sim.rds')
pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs', 'sigma_proc',
         'w[1]', 'w[2]', 'w[3]', 'w[4]', 'w[5]')

print(mod, pars)
t <- traceplot(mod, ncol = 2, pars)
p <- plot_post_sim(mod, pars,
                   vals = c( beta_s, beta_p, K_20*100, sigma_obs, sigma_proc*10, w),
                   xlim = c(0,1.5))
png('figures/simulation_fits/SAM5_detC_logpi_sim_newriver.png', height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery SAM5 det C (New River)',
                                                      size = 14))
dev.off()

pp <- calc_pp_ests(mod, newR, pars= pars, factors = c(1, 1, .01, 1, .1, rep(1, 5)),
                   sim_func = simulate_detC_SAM5)
png('figures/simulation_fits/PPcheck_SAM5_detC_logpi_newriver_sim.png')
    plot_pp_interval(sim_new$R_obs, pp, xrng = c(1,365))
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_new,
                'figures/simulation_fits/PPcheck_SAM5_detC_logpi_newriver.png')

mod <- readRDS('src/SAM/stan/fits/SAM5_detC_logpi_clack_sim.rds')
pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs', 'sigma_proc',
         'w[1]', 'w[2]', 'w[3]', 'w[4]', 'w[5]')

print(mod, pars)
t <- traceplot(mod, ncol = 2, pars)
p <- plot_post_sim(mod, pars,
                   vals = c( beta_s, beta_p, K_20*100, sigma_obs, sigma_proc*10, w),
                   xlim = c(0,1.5))
png('figures/simulation_fits/SAM5_detC_logpi_sim_clack.png', height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery SAM5 det C (Clackamas)',
                                                          size = 14))
dev.off()

pp <- calc_pp_ests(mod, clack, pars= pars, factors = c(1, 1, .01, 1, .1, rep(1, 5)),
                   sim_func = simulate_detC_SAM5)
png('figures/simulation_fits/PPcheck_SAM5_detC_logpi_clack_sim.png')
    plot_pp_interval(sim_clack$R_obs, pp, xrng = c(1,365))
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_clack,
                'figures/simulation_fits/PPcheck_SAM5_detC_logpi_clackamas.png')

# Run on actual rivers ####
# Clackamas
sim_clack <- simulate_detC_SAM5(clack)

stan_dat <- list(ndays = nrow(sim_clack), nweights = 5, R_obs = sim_clack$ER,
                 P = sim_clack$GPP, C0 = 100, temp = sim_clack$temp_C,
                 litter = sim_clack$litter, Q = sim_clack$Q)
mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_clack.rds')

t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                       'sigma_proc',  'w'))
p <- plot(mod, pars = c('beta_s', 'beta_p', 'K_20',
                        'sigma_obs', 'sigma_proc', 'w'))

png('figures/simulation_fits/SAM5_detC_logpi_clack.png', height = 300, width = 650)
ggpubr::ggarrange(p, t)
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_clack,
                     'figures/simulation_fits/PPcheck_SAM5_detC_logpi_clackamas.png')


# Pecos
sim_pecos <- simulate_detC_SAM5(pecos)

stan_dat <- list(ndays = nrow(sim_pecos), nweights = 5, R_obs = sim_pecos$ER,
                 P = sim_pecos$GPP, C0 = 100, temp = sim_pecos$temp_C,
                 litter = sim_pecos$litter, Q = sim_pecos$Q)
mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_pecos.rds')

t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                       'sigma_proc',  'w'))
p <- plot(mod, pars = c('beta_s', 'beta_p', 'K_20',
                        'sigma_obs', 'sigma_proc', 'w'))

png('figures/simulation_fits/SAM5_detC_logpi_pecos.png', height = 300, width = 650)
ggpubr::ggarrange(p, t)
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_pecos,
                     'figures/simulation_fits/PPcheck_SAM5_detC_logpi_pecos.png')


# Grand
sim_grand <- simulate_detC_SAM5(grand)

stan_dat <- list(ndays = nrow(sim_grand), nweights = 5, R_obs = sim_grand$ER,
                 P = sim_grand$GPP, C0 = 100, temp = sim_grand$temp_C,
                 litter = sim_grand$litter, Q = sim_grand$Q)
mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_grand.rds')

t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                       'sigma_proc',  'w'))
p <- plot(mod, pars = c('beta_s', 'beta_p', 'K_20',
                        'sigma_obs', 'sigma_proc', 'w'))

png('figures/simulation_fits/SAM5_detC_logpi_grand.png', height = 300, width = 650)
ggpubr::ggarrange(p, t)
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_grand,
                     'figures/simulation_fits/PPcheck_SAM5_detC_logpi_grand.png')



# Fanno
sim_fanno <- simulate_detC_SAM5(fanno)

stan_dat <- list(ndays = nrow(sim_fanno), nweights = 5, R_obs = sim_fanno$ER,
                 P = sim_fanno$GPP, C0 = 100, temp = sim_fanno$temp_C,
                 litter = sim_fanno$litter, Q = sim_fanno$Q)
mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_fanno.rds')

t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                       'sigma_proc',  'w'))
p <- plot(mod, pars = c('beta_s', 'beta_p', 'K_20',
                        'sigma_obs', 'sigma_proc', 'w'))

png('figures/simulation_fits/SAM5_detC_logpi_fanno.png', height = 300, width = 650)
ggpubr::ggarrange(p, t)
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_fanno,
                     'figures/simulation_fits/PPcheck_SAM5_detC_logpi_fanno.png')


# East
sim_east <- simulate_detC_SAM5(east)

stan_dat <- list(ndays = nrow(sim_east), nweights = 5, R_obs = sim_east$ER,
                 P = sim_east$GPP, C0 = 100, temp = sim_east$temp_C,
                 litter = sim_east$litter, Q = sim_east$Q)
mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_east.rds')
mod <- readRDS('src/SAM/stan/fits/SAM5_detC_logpi_east.rds')

t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                       'sigma_proc',  'w'))
p <- plot(mod, pars = c('beta_s', 'beta_p', 'K_20',
                        'sigma_obs', 'sigma_proc', 'w'))

png('figures/simulation_fits/SAM5_detC_logpi_east.png', height = 300, width = 650)
ggpubr::ggarrange(p, t)
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_east,
                     'figures/simulation_fits/PPcheck_SAM5_detC_logpi_east.png')


# newR
sim_newR <- simulate_detC_SAM5(newR)

stan_dat <- list(ndays = nrow(sim_newR), nweights = 5, R_obs = sim_newR$ER,
                 P = sim_newR$GPP, C0 = 100, temp = sim_newR$temp_C,
                 litter = sim_newR$litter, Q = sim_newR$Q)
mod <- stan('src/SAM/stan/SAM5_detC_logpi.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
# beep(sound = 8)
saveRDS(mod, 'src/SAM/stan/fits/SAM5_detC_logpi_newR.rds')

t <- traceplot(mod, ncol = 2, pars = c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                       'sigma_proc',  'w'))
p <- plot(mod, pars = c('beta_s', 'beta_p', 'K_20',
                        'sigma_obs', 'sigma_proc', 'w'))

png('figures/simulation_fits/SAM5_detC_logpi_newR.png', height = 300, width = 650)
ggpubr::ggarrange(p, t)
dev.off()

plot_ER_ppcheck_SAM5(mod, sim_newR,
                     'figures/simulation_fits/PPcheck_SAM5_detC_logpi_newR.png')



test_params <- extract(mod, c('beta_s', 'beta_p', 'K_20', 'sigma_obs',
                                    'sigma_proc',  'w'))


