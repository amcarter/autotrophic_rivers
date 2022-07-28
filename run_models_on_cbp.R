# Test SAM stan models with simulated data
library(tidyverse)
library(beepr)
library(rstan)
library(shinystan)
source('src/models/simulate/simulation_functions.R')
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
dat <- read_csv('../nhc_50yl/data/met_and_drivers_for_model.csv')
mdat <- read_csv('../nhc_50yl/data/siteData/NHCsite_metadata.csv') %>%
    select(site = sitecode, slope = slope_wbx)
dat <- dat %>% left_join(mdat) %>%
    select(-PAR_surface)
LAI <- read_csv('../nhc_50yl/data/light/drivers/daily_modeled_light_all_sites.csv')
dat <- dat %>% left_join(LAI, by = c('site', 'date'))
# simulate data
prep_data <- function(dat){
    dd <- dat %>%
        arrange(date) %>%
        select(site, date, GPP, ER, ER.upper, ER.lower, slope, depth,
               discharge, temp_C = temp.water, PAR_surface, LAI) %>%
        complete(date = seq(min(date), max(date), by = 'day'))%>%
        mutate(litter = calc_litter_oneyear(500, LAI),
               across(c(-site, -date),
                      zoo::na.approx, na.rm = F),
               light = PAR_surface/max(PAR_surface, na.rm = T),
               Q = discharge/max(discharge, na.rm = T),
               tau_gm2 = calc_tau_gm2(depth, slope),
               tau = tau_gm2/max(tau_gm2, na.rm = T)) %>%
        filter(!is.na(ER))

    return(dd)
}

cbp <- filter(dat, site == "CBP") %>%
    prep_data()

nhc <- filter(dat, site == "NHC") %>%
    prep_data()
# data simulation ####
E_a = 0.63            # activation energy for heterotrophic respiration
k_b = 8.6173 * 10^-5  # Boltzmann's constant in eV/K
ARf = 0.44            # the fraction of GPP respired by autotrophs
# Detrital Carbon Models: ####

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
    ARf = 0.44
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
sim_cbp <- simulate_detC_logpi(cbp)
sim_cbp %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer(-date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)
pars = c('beta_s', 'sigma_obs', 'sigma_proc', 'K_20')

stan_dat <- list(ndays = nrow(sim_cbp), R_obs = sim_cbp$R_obs, P = sim_cbp$GPP,
                 C0 = 100, temp = sim_cbp$temp_C, litter = sim_cbp$litter,
                 Q = sim_cbp$Q)
mod <- stan('src/models/stan/detC_logpi_1.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/models/stan/fits/detC_logpi_1_cbp_sim.rds')

# mod <- readRDS('src/models/stan/fits/detC_logpi_1_cbp_sim.rds')
# pp <- calc_pp_ests(mod, cbp, pars= pars, factors = c(1, 1, .1, .01),
#                    sim_func = simulate_detC_logpi)
#
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/SAM/stan/fits/detC_logpi_1_cbp_sim_pp.rds')
m <- readRDS('src/models/stan/fits/detC_logpi_1_cbp_sim_pp.rds')
mod <- m$mod
pp <- m$post_preds
png('figures/nhc/PPcheck_detC_logpi_cbp_sim.png')
    plot_pp_interval(sim_cbp$R_obs, pp)#, ylim = c(-30,0))
dev.off()

print(mod, pars)
t <- traceplot(mod, ncol = 1, pars)
p <- plot_post_sim(mod, pars, vals = c(beta_s, sigma_obs, sigma_proc*10, K_20*100),
                   xlim = c(0,1.5))
png('figures/nhc/detC_logpi_sim_cbp.png',
    height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery det C log process error (CBP)',
                                                      size = 14))
dev.off()




# Model: detC_logpi_ss ####
# define parameters
C0 = 100         # Initial organic C
K_20 = .01       # Heterotrophic respiration at 20 C
beta_s = 0.8     # Percent removal of organic carbon from max storm
tau0 = 0.5       # lowest shear stress that causes a disturbance
sigma_proc = .02 # lognormally distributed, so this is a % error
sigma_obs = 0.08

simulate_detC_logpi_ss <- function(ss){
    ss <- ss %>%
        mutate(C = 0.5,
               R = 0.5,
               R_obs = ER) %>%
        select(date, GPP, ER, R_obs, R, light, Q, temp_C, litter, tau, C)

    ndays <- nrow(ss) # number of days
    ss$K = calc_rate_coef(ss$temp_C, K_20 = K_20)
    ss$disturb = calc_disturbance(ss$tau, tau0)

    Chat = numeric()
    Chat[1] = C0
    ss$C[1] = exp(rnorm(1, log(C0), sigma_proc))
    ARf = 0.44
    ss$AR = -ARf * ss$GPP
    ss$R[1] = ss$AR[1] - ss$K[1] * ss$C[1]

    for(i in 2:ndays){
        Chat[i] = ss$C[i-1]*(1-0.8 *ss$disturb[i]) +
            ss$litter[i] + ss$R[i-1] + ss$GPP[i-1]
        ss$C[i] = exp(rnorm(1, log(Chat[i]), sigma_proc))
        ss$R[i] = ss$AR[i] - ss$K[i]*ss$C[i]
    }
    # observation model:
    ss$R_obs = rnorm(ndays, ss$R, sigma_obs)

    return(ss)
}

# Simulate and run
sim_cbp <- simulate_detC_logpi_ss(cbp)
sim_cbp %>% select(date, GPP, disturb, R, ER,C, litter )%>%
    pivot_longer(-date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)
pars = c( 'tau0', 'sigma_obs', 'sigma_proc_scaled', 'K_20_scaled')

stan_dat <- list(ndays = nrow(sim_cbp), R_obs = sim_cbp$R_obs, P = sim_cbp$GPP,
                 C0 = 100, temp = sim_cbp$temp_C, litter = sim_cbp$litter,
                 tau = sim_cbp$tau)
mod <- stan('src/models/stan/detC_logpi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/models/stan/fits/detC_logpi_ss_cbp_sim.rds')

# # mod <- readRDS('src/models/stan/fits/detC_logpi_ss_cbp_sim.rds')
# pp <- calc_pp_ests(mod, cbp, pars= pars, factors = c(1, 1, .1, .01),
#                    sim_func = simulate_detC_logpi_ss)
#
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/models/stan/fits/detC_logpi_ss_cbp_sim_pp.rds')
m <- readRDS('src/models/stan/fits/detC_logpi_ss_cbp_sim_pp.rds')
mod <- m$mod
pp <- m$post_preds
png('figures/nhc/PPcheck_detC_logpi_ss_cbp_sim.png')
    plot_pp_interval(sim_cbp$R_obs, pp)#, ylim = c(-30,0))
dev.off()

print(mod, pars)
t <- traceplot(mod, ncol = 1, pars)
p <- plot_post_sim(mod, pars, vals = c(beta_s, sigma_obs, sigma_proc*10, K_20*100),
                   xlim = c(0,1.5))
png('figures/nhc/detC_logpi_ss_sim_cbp.png',
    height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery det C log pe, shear stress (CBP)',
                                                      size = 14))
dev.off()

cmod <- apply(data.frame(extract(mod, pars = 'C')), 2, median)
plot(sim_cbp$date, sim_cbp$C, type = 'l')
points(sim_cbp$date, cmod)
mean(sim_cbp$C)
mean(cmod)*2.89
glimpse(cmod)
fit$summary[2:(nrow(fit$summary)-2),1]

# Model: detC ####
# define parameters
C0 = 100       # Initial organic C
pars = list(
    K_20 = .01,        # Heterotrophic respiration at 20 C
    beta_s = 0.8,   # Percent removal of organic carbon from max storm
    sigma_proc = 0.02, # lognormally distributed, so this is a % error
    sigma_obs = 0.08

)
simulate_detC <- function(ss, pars){
    ss <- ss %>%
        mutate(logC = 0.5,
               R = 0.5,
               R_obs = ER) %>%
        select(date, GPP, ER, R_obs, R, light, Q, temp_C, litter, logC)

    ndays <- nrow(ss) # number of days
    ss$K = calc_rate_coef(ss$temp_C, K_20 = pars$K_20)

    Chat = numeric()
    Chat[1] = C0
    ss$logC[1] = rnorm(1, log(C0), pars$sigma_proc)
    ARf = 0.44
    ss$AR = -ARf * ss$GPP
    ss$R[1] = ss$AR[1] - ss$K[1] * exp(ss$logC[1])

    for(i in 2:ndays){
        Chat[i] = (exp(ss$logC[i-1]) + ss$litter[i] + ss$R[i-1] + ss$GPP[i-1])*
            (1-pars$beta_s*ss$Q[i])
        ss$logC[i] = rnorm(1, log(Chat[i]), pars$sigma_proc)
        ss$R[i] = ss$AR[i] - ss$K[i]*exp(ss$logC[i])
    }
    # observation model:
    ss$R_obs = rnorm(ndays, ss$R, sigma_obs)

    return(ss)
}

# Simulate and run
sim_cbp <- simulate_detC(cbp, pars)
sim_cbp %>% select(date, GPP, Q, R, ER, logC, litter )%>%
    pivot_longer(-date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_cbp), R_obs = sim_cbp$R_obs, P = sim_cbp$GPP,
                 C0 = 100, temp = sim_cbp$temp_C, litter = sim_cbp$litter,
                 Q = sim_cbp$Q)
mod <- stan('src/models/stan/detC.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/models/stan/fits/detC_cbp_sim.rds')

# mod <- readRDS('src/models/stan/fits/detC_cbp_sim.rds')
# pp <- calc_pp_ests(mod, cbp, pars= pars, factors = c(1, 1, .1, .01),
#                    sim_func = simulate_detC)
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/models/stan/fits/detC_sim_cbp_pp.rds')
m <- readRDS('src/models/stan/fits/detC_sim_cbp_pp.rds')
mod <- m$mod
pp <- m$post_preds

print(mod, names(pars))
t <- traceplot(mod, ncol = 1, names(pars))
p <- plot_post_sim(mod, names(pars),
                   vals = c( pars$beta_s, pars$sigma_obs,
                             pars$sigma_proc*10, pars$K_20*100),
                   xlim = c(0,1.5))
png('figures/nhc/detC_sim_cbp.png', height = 300, width = 650)
    plot <- ggpubr::ggarrange(p, t)
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob('Parameter recovery det logC (CBP)',
                                                      size = 14))
dev.off()

png('figures/nhc/PPcheck_detC_cbp_sim.png')
    plot_pp_interval(sim_cbp$R_obs, pp)#, ylim = c(-30,0))
dev.off()



# Run on actual data ####

stan_dat <- list(ndays = nrow(cbp), R_obs = cbp$ER, P = cbp$GPP,
                 C0 = 100, temp = cbp$temp_C, litter = cbp$litter,
                 Q = cbp$Q)
mod <- stan('src/models/stan/detC.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            warmup = 500, iter = 1000)
beep(sound = 8)
saveRDS(mod, 'src/models/stan/fits/detC_cbp.rds')

# mod <- readRDS('src/models/stan/fits/detC_cbp.rds')
pp <- calc_pp_ests(mod, cbp, pars= pars, factors = c(1, 1, .1, .01),
                   sim_func = simulate_detC)
saveRDS(list(mod = mod, post_preds = pp),
        'src/models/stan/fits/detC_cbp_pp.rds')
m <- readRDS('src/models/stan/fits/detC_cbp_pp.rds')
mod <- m$mod
pp <- m$post_preds

print(mod, names(pars))

png('figures/nhc/PPcheck_detC_cbp.png')
    plot_pp_interval(sim_cbp$ER, pp)#, ylim = c(-30,0))
dev.off()

a <- traceplot(mod, ncol = 1, (pars))
b <- plot(mod, pars = (pars))
c <- pairs(mod, pars = (pars))

png('figures/simulation_fits/rivers/detC_logpi_ss_cbp.png', height = 800, width = 400)
    ggpubr::ggarrange(a, b, ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
fit_pars <- fit_pars[-5]
cbp_fit <- simulate_detC(cbp, fit_pars)
plot(cbp_fit$date, cbp_fit$R_obs, type = 'l')
lines(cbp_fit$date, cbp_fit$R_obs, col = 2)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_clack_fit2.png',
    width = 700, height = 500)
    plot_ER_breakdown(cbp_fit, fit_pars)#, xrng = c(31, 1445))
dev.off()

clackant <- as.data.frame(ANT)
clack_fit$ant1 = clackant[,1]
clack_fit$ant2 = clackant[,2]
clack_fit$ant3 = clackant[,3]
clack_fit$ant4 = clackant[,4]
clack_fit$LAI = clack$LAI
clack_fit$discharge = clack$discharge
dd <- slice(clack_fit, 31: 1445)
png('figures/GPP.png', width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
plot(dd$date, dd$GPP, type = 'l', col = 'forestgreen', xlab = 'Date',
     ylab = expression(paste("GPP (g"~O[2]~"m"^"-2"~" d"^"-1"*")")), lwd = 1.5)
dev.off()
png('figures/antGPP.png',  width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
plot(dd$date, dd$GPP, type = 'n', col = 'forestgreen', xlab = 'Date',
     ylab = "antecedent GPP", lwd = 1.5)
lines(dd$date, dd$ant1, col = alpha('goldenrod'))
lines(dd$date, dd$ant2, col = alpha('forestgreen', .5))
lines(dd$date, dd$ant3, col = alpha('forestgreen'))
lines(dd$date, dd$ant4)
legend('top', c(expression(paste(Delta,'t:')), '1-2 days', '3-6 days',
                '7-14 days', '15-30 days'), cex = 1.1,
       col = c(NA, alpha('goldenrod'), alpha('forestgreen', 0.5),
               'forestgreen', 1), lty = c(NA, 1, 1, 1, 1), bty = 'n', ncol = 5)
dev.off()
png('figures/litter.png',  width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
plot(dd$date, dd$LAI, type = 'l', xlab = 'Date',
     ylab = 'LAI', lwd = 1.5, ylim = c(0.9, 3.4))
par(new = T)
plot(dd$date, dd$litter, pch = 19, col = '#F4B570',
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend('topleft', c('LAI', 'litterfall'), cex = 1.1,
       col = c(1, '#F4B570'), pch = c(NA, 19), lty = c(1, NA),
       bty = 'n', ncol = 2, inset = .05)

dev.off()
png('figures/Q.png', width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
plot(dd$date, dd$discharge, type = 'l', col = 'steelblue', xlab = 'Date',
     ylab = expression(paste("Discharge (m"^"3"~" s"^"-1"*")")), lwd = 1.5)
dev.off()
png('figures/temp.png',  width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
plot(dd$date, dd$temp_C, type = 'l',xlab = 'Date',
     ylab = expression(paste("Temperature (",degree,"C)")), lwd = 1.5)
dev.off()
png('figures/C.png',  width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
plot(dd$date, dd$C, type = 'l',xlab = 'Date',
     ylab = expression(paste("Carbon (g m"^"-2"*")")), lwd = 1.5)
legend('topleft', 'Latent State: \nCarbon Storage', inset = .05, bty = 'n')
dev.off()
png('figures/tau.png',  width = 650, height = 250)
par(mar = c(5, 5, 0, 1))
# plot(dd$date, dd$discharge, type = 'n', col = 'steelblue', xlab = 'Date',
#      ylab = expression(paste("Discharge (m"^"3"~" s"^"-1"*")")), lwd = 1.5)
# polygon(c(dd$date, rev(dd$date)), c(dd$discharge, rep(0, nrow(dd))),
#         col = 'lightblue', border = NA)
# par(new = T)
plot(dd$date, dd$tau, type = 'l',xlab = 'Date',
     ylab = expression(paste(tau, ' (shear stress)')),
     lwd = 1.5)
# ylab = expression(paste(tau, " (shear stress, g m"^"-2"*")")))
abline(h = 0.7, lty = 2, col = 'brown3')
legend('topright', c(expression(tau), expression(tau[0])),
       col = c(1, 'brown3'), fill = c(NA, NA), border = NA,
       lty = c(1,2), bty = 'n', ncol = 3, cex = 1.1)

dev.off()

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_clack_PPcheck.png')
plot_pp_interval(sim_clack$ER, pp, xrng = c(31, nrow(clack)))#, ylim = c(-3,0.5))
dev.off()
#
# plot_ER_ppcheck_SAM5(mod, sim_clack,
#                 'figures/simulation_fits/rivers/PPcheck_SAM5_detC_logpi_clack.png')

# Pecos
ANT <- as.matrix(calc_antecedent_drivers(pecos$GPP, 4, 2))
sim_pecos <- simulate_SAMint_detC_logpi_ss(pecos, pars, ANT, antdays, 10)
sim_pecos %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_pecos), nweights = nweights, antdays = antdays,
#                  R_obs = sim_pecos$ER, P = sim_pecos$GPP, C0 = 10,
#                  temp = sim_pecos$temp_C, litter = sim_pecos$litter,
#                  tau = sim_pecos$tau, ANT = ANT)
# mod <- stan('src/SAM/stan/SAMint_detC_logpi_ss.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             iter = 2000, control = list(max_treedepth = 12))
# saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_ss_pecos.rds')
# beep(sound = 8)
# # mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_pecos.rds')
# pp <- calc_pp_ests_SAM(mod, pecos, pars = pars, ANT = ANT, antdays = antdays, C0 = 10,
#                        sim_func = simulate_SAMint_detC_logpi_ss)
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/SAM/stan/fits/SAMint_detC_logpi_ss_pecos_pp.rds')
# m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_pecos_pp.rds')
mod <- m$mod
pp <- m$post_preds

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_pecos_PPcheck.png')
plot_pp_interval(sim_pecos$ER, pp, xrng = c(31, nrow(pecos)))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_pecos.png', height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
pecos_fit <- simulate_SAMint_detC_logpi_ss(pecos, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_pecos_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(pecos_fit, fit_pars)
dev.off()

# Grand
ANT <- as.matrix(calc_antecedent_drivers(grand$GPP, 4, 2))
sim_grand <- simulate_SAMint_detC_logpi_ss(grand, pars, ANT, antdays, C0 = 10)
sim_grand %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)
#
# stan_dat <- list(ndays = nrow(sim_grand), nweights = nweights, antdays = antdays,
#                  R_obs = sim_grand$ER, P = sim_grand$GPP, C0 = 10,
#                  temp = sim_grand$temp_C, litter = sim_grand$litter,
#                  tau = sim_grand$tau, ANT = ANT)
# mod <- stan('src/SAM/stan/SAMint_detC_logpi_ss.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             iter = 2000, control = list(max_treedepth = 12))
# saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_ss_grand.rds')
# beep(sound = 8)
#
# # mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_grand.rds')
# pp <- calc_pp_ests_SAM(mod, grand, pars = pars, ANT = ANT, antdays = antdays, C0 = 10,
#                        sim_func = simulate_SAMint_detC_logpi_ss)
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/SAM/stan/fits/SAMint_detC_logpi_ss_grand_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_grand_pp.rds')
mod <- m$mod
pp <- m$post_preds

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_grand_PPcheck.png')
plot_pp_interval(sim_grand$ER, pp, xrng = c(31, nrow(grand)))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_grand.png', height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()


fit_pars <- get_pars(mod, names(pars), as_list = T)
grand_fit <- simulate_SAMint_detC_logpi_ss(grand, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_grand_fit.png',
    width = 800, height = 600)
plot_ER_breakdown(grand_fit, fit_pars)
dev.off()

# Fanno
ANT <- as.matrix(calc_antecedent_drivers(fanno$GPP, 4, 2))
sim_fanno <- simulate_SAMint_detC_logpi_ss(fanno, pars, ANT, antdays)
sim_fanno %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)
#
# stan_dat <- list(ndays = nrow(sim_fanno), nweights = nweights, antdays = antdays,
#                  R_obs = sim_fanno$ER, P = sim_fanno$GPP, C0 = 100,
#                  temp = sim_fanno$temp_C, litter = sim_fanno$litter,
#                  tau = sim_fanno$tau, ANT = ANT)
# mod <- stan('src/SAM/stan/SAMint_detC_logpi_ss.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             iter = 2000, control = list(max_treedepth = 12))
# saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_ss_fanno.rds')
# beep(sound = 8)
# mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_fanno.rds')
# pp <- calc_pp_ests_SAM(mod, fanno, pars = pars, ANT = ANT, antdays = antdays,
#                        sim_func = simulate_SAMint_detC_logpi_ss)
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/SAM/stan/fits/SAMint_detC_logpi_ss_fanno_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_fanno_pp.rds')
mod <- m$mod
pp <- m$post_preds

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_fanno_PPcheck.png')
plot_pp_interval(sim_fanno$R_obs, pp, xrng = c(365, 730))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_fanno.png', height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
fanno_fit <- simulate_SAMint_detC_logpi_ss(fanno, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_fanno_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(fanno_fit, fit_pars)
dev.off()

# East
ANT <- as.matrix(calc_antecedent_drivers(east$GPP, 4, 2))
sim_east <- simulate_SAMint_detC_logpi_ss(east, pars, ANT, antdays)
sim_east %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_east), nweights = nweights, antdays = antdays,
#                  R_obs = sim_east$ER, P = sim_east$GPP, C0 = 100,
#                  temp = sim_east$temp_C, litter = sim_east$litter,
#                  tau = sim_east$tau, ANT = ANT)
# mod <- stan('src/SAM/stan/SAMint_detC_logpi_ss.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             iter = 2000, control = list(max_treedepth = 12))
# saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_ss_east.rds')
# beep(sound = 8)
# mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_east.rds')
# pp <- calc_pp_ests_SAM(mod, east, pars = pars, ANT = ANT, antdays = antdays,
#                        sim_func = simulate_SAMint_detC_logpi_ss)
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/SAM/stan/fits/SAMint_detC_logpi_ss_east_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_east_pp.rds')
mod <- m$mod
pp <- m$post_preds

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_east_PPcheck.png')
plot_pp_interval(sim_east$ER, pp, xrng = c(732, 1462))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_east.png', height = 800, width = 400)
ggpubr::ggarrange(a, b,ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
east_fit <- simulate_SAMint_detC_logpi_ss(east, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_east_fit.png',
    width = 800, height = 600)
plot_ER_breakdown(east_fit, fit_pars, xrng = c(732, 1462))
dev.off()

# newR
ANT <- as.matrix(calc_antecedent_drivers(newR$GPP, 4, 2))
sim_new <- simulate_SAMint_detC_logpi_ss(newR, pars, ANT, antdays)
sim_new %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

# stan_dat <- list(ndays = nrow(sim_new), nweights = nweights, antdays = antdays,
#                  R_obs = sim_new$ER, P = sim_new$GPP, C0 = 100,
#                  temp = sim_new$temp_C, litter = sim_new$litter,
#                  tau = sim_new$tau, ANT = ANT)
# mod <- stan('src/SAM/stan/SAMint_detC_logpi_ss.stan',
#             data = stan_dat,
#             chains = 4,  cores = 4,
#             iter = 2000, control = list(max_treedepth = 12))
# saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_ss_new.rds')
# beep(sound = 8)
# mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_new.rds')
# pp <- calc_pp_ests_SAM(mod, newR, pars = pars, ANT = ANT, antdays = antdays,
#                        sim_func = simulate_SAMint_detC_logpi_ss)
# saveRDS(list(mod = mod, post_preds = pp),
#         'src/SAM/stan/fits/SAMint_detC_logpi_ss_new_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_new_pp.rds')
mod <- m$mod
pp <- m$post_preds

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_new_PPcheck.png')
plot_pp_interval(sim_new$ER, pp, xrng = c(31, nrow(newR)))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_new.png', height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()


fit_pars <- get_pars(mod, names(pars), as_list = T)
new_fit <- simulate_SAMint_detC_logpi_ss(newR, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_ss_new_fit.png',
    width = 800, height = 600)
plot_ER_breakdown(new_fit, fit_pars)
dev.off()


# Run with fixed Obs error####
#clackamas
ANT <- as.matrix(calc_antecedent_drivers(clack$GPP, 4, 2))
obs_error = median(-clack$ER.lower + clack$ER.upper)/4
obs_error_sd = sd(-clack$ER.lower + clack$ER.upper)/4
sim_clack <- simulate_SAMint_detC_logpi_ss(clack, pars, ANT, antdays)
sim_clack %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_clack), nweights = nweights, antdays = antdays,
                 R_obs = sim_clack$ER, P = sim_clack$GPP, C0 = 100,
                 temp = sim_clack$temp_C, litter = sim_clack$litter,
                 tau = sim_clack$tau, ANT = ANT, sigma_obs_prior = obs_error,
                 sigma_obs_sd = obs_error_sd)
mod <- stan('src/SAM/stan/SAMint_detC_logpi_fixedoi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            iter = 2000, control = list(max_treedepth = 12))
saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_clack.rds')
beep(sound = 8)

mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_clack.rds')
pp <- calc_pp_ests_SAM(mod, clack, pars = pars, ANT = ANT, antdays = antdays,
                       sim_func = simulate_SAMint_detC_logpi_ss)
saveRDS(list(mod = mod, post_preds = pp),
        'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_clack_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_clack_pp.rds')
mod <- m$mod
pp <- m$post_preds

a <- traceplot(mod, ncol = 2, pars = names(pars))#c('K_20', 'tau0', 'beta_p', 'sigma_proc','w') )
b <- plot(mod,  pars = names(pars))#c('K_20', 'tau0', 'beta_p', 'sigma_proc','w') )
c <- pairs(mod, pars = names(pars))#c('K_20', 'tau0', 'beta_p', 'sigma_proc','lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_clack.png',
    height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_clack_PPcheck.png')
plot_pp_interval(sim_clack$R_obs, pp, xrng = c(365, 730))#, ylim = c(-3,0.5))
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
clack_fit <- simulate_SAMint_detC_logpi_ss(clack, fit_pars, ANT, antdays)

plot_ER_breakdown(clack_fit, pars)

#
# plot_ER_ppcheck_SAM5(mod, sim_clack,
#                 'figures/simulation_fits/rivers/PPcheck_SAM5_detC_logpi_clack.png')

# Pecos
ANT <- as.matrix(calc_antecedent_drivers(pecos$GPP, 4, 2))
sim_pecos <- simulate_SAMint_detC_logpi_ss(pecos, pars, ANT, antdays)
obs_error = median(-pecos$ER.lower + pecos$ER.upper)/4
obs_error_sd = sd(-pecos$ER.lower + pecos$ER.upper)/4

# sim_pecos %>% select(date, GPP, disturb, R, ER,C, litter )%>%
#     pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
#     ggplot(aes(date, value, col = var)) +
#     geom_line()+
#     facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_pecos), nweights = nweights, antdays = antdays,
                 R_obs = sim_pecos$ER, P = sim_pecos$GPP, C0 = 100,
                 temp = sim_pecos$temp_C, litter = sim_pecos$litter,
                 tau = sim_pecos$tau, ANT = ANT, sigma_obs_prior = obs_error,
                 sigma_obs_sd = obs_error_sd)
mod <- stan('src/SAM/stan/SAMint_detC_logpi_fixedoi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            iter = 1000, control = list(max_treedepth = 12))
saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_pecos.rds')
beep(sound = 8)
pp <- calc_pp_ests_SAM(mod, pecos, pars = pars, ANT = ANT, antdays = antdays,
                       sim_func = simulate_SAMint_detC_logpi_ss)
saveRDS(list(mod = mod, post_preds = pp),
        'src/SAM/stan/fits/SAMint_detC_logpi_ss_pecos_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_ss_pecos_pp.rds')
mod <- m$mod
pp <- m$post_preds
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_pecos_PPcheck.png')
plot_pp_interval(sim_pecos$ER, pp)#, xrng = c(365, 730))#, ylim = c(-3,0.5))
dev.off()

a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_pecos.png',
    height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
pecos_fit <- simulate_SAMint_detC_logpi_ss(pecos, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_pecos_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(pecos_fit, fit_pars)
dev.off()

# Grand
ANT <- as.matrix(calc_antecedent_drivers(grand$GPP, 4, 2))
sim_grand <- simulate_SAMint_detC_logpi_ss(grand, pars, ANT, antdays)
obs_error = median(-grand$ER.lower + grand$ER.upper)/4
obs_error_sd = sd(-grand$ER.lower + grand$ER.upper)/4

sim_grand %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_grand), nweights = nweights, antdays = antdays,
                 R_obs = sim_grand$ER, P = sim_grand$GPP, C0 = 100,
                 temp = sim_grand$temp_C, litter = sim_grand$litter,
                 tau = sim_grand$tau, ANT = ANT, sigma_obs_prior = obs_error,
                 sigma_obs_sd = obs_error_sd)
mod <- stan('src/SAM/stan/SAMint_detC_logpi_fixedoi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            iter = 2000, control = list(max_treedepth = 12))
saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_grand.rds')
beep(sound = 8)

pp <- calc_pp_ests_SAM(mod, grand, pars = pars, ANT = ANT, antdays = antdays,
                       sim_func = simulate_SAMint_detC_logpi_ss)
saveRDS(list(mod = mod, post_preds = pp),
        'src/SAM/stan/fits/SAMint_detC_logpi_fixedio_ss_grand_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedio_ss_grand_pp.rds')
mod <- m$mod
pp <- m$post_preds
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_grand_PPcheck.png')
plot_pp_interval(sim_grand$ER, pp)#, xrng = c(365, 730))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_grand.png',
    height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
grand_fit <- simulate_SAMint_detC_logpi_ss(grand, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_grand_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(grand_fit, fit_pars)
dev.off()

# Fanno
ANT <- as.matrix(calc_antecedent_drivers(fanno$GPP, 4, 2))
sim_fanno <- simulate_SAMint_detC_logpi_ss(fanno, pars, ANT, antdays)
obs_error = median(-fanno$ER.lower + fanno$ER.upper)/4
obs_error_sd = sd(-fanno$ER.lower + fanno$ER.upper)/4

sim_fanno %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_fanno), nweights = nweights, antdays = antdays,
                 R_obs = sim_fanno$ER, P = sim_fanno$GPP, C0 = 100,
                 temp = sim_fanno$temp_C, litter = sim_fanno$litter,
                 tau = sim_fanno$tau, ANT = ANT, sigma_obs_prior = obs_error,
                 sigma_obs_sd = obs_error_sd)
mod <- stan('src/SAM/stan/SAMint_detC_logpi_fixedoi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            iter = 1000, control = list(max_treedepth = 12))
saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_fanno.rds')
mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_fanno.rds')
beep(sound = 8)

pp <- calc_pp_ests_SAM(mod, fanno, pars = pars, ANT = ANT, antdays = antdays,
                       sim_func = simulate_SAMint_detC_logpi_ss)
saveRDS(list(mod = mod, post_preds = pp),
        'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_fanno_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_fanno_pp.rds')
mod <- m$mod
pp <- m$post_preds
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_fanno_PPcheck.png')
plot_pp_interval(sim_fanno$ER, pp)#, xrng = c(365, 730))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_fanno.png',
    height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
fanno_fit <- simulate_SAMint_detC_logpi_ss(fanno, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_fanno_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(fanno_fit, fit_pars)
dev.off()

# East
ANT <- as.matrix(calc_antecedent_drivers(east$GPP, 4, 2))
sim_east <- simulate_SAMint_detC_logpi_ss(east, pars, ANT, antdays)
obs_error = median(-east$ER.lower + east$ER.upper)/4
obs_error_sd = sd(-east$ER.lower + east$ER.upper)/4

sim_east %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_east), nweights = nweights, antdays = antdays,
                 R_obs = sim_east$ER, P = sim_east$GPP, C0 = 100,
                 temp = sim_east$temp_C, litter = sim_east$litter,
                 tau = sim_east$tau, ANT = ANT, sigma_obs_prior = obs_error,
                 sigma_obs_sd = obs_error_sd)
mod <- stan('src/SAM/stan/SAMint_detC_logpi_fixedoi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            iter = 2000, control = list(max_treedepth = 12))
saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_east.rds')
beep(sound = 8)
pp <- calc_pp_ests_SAM(mod, east, pars = pars, ANT = ANT, antdays = antdays,
                       sim_func = simulate_SAMint_detC_logpi_ss)
saveRDS(list(mod = mod, post_preds = pp),
        'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_east_pp.rds')
m <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_east_pp.rds')
mod <- m$mod
pp <- m$post_preds

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_east_PPcheck.png')
plot_pp_interval(sim_east$ER, pp)#, xrng = c(365, 730))#, ylim = c(-3,0.5))
dev.off()
a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_east.png',
    height = 800, width = 400)
ggpubr::ggarrange(a, b,ncol = 1)
dev.off()

fit_pars <- get_pars(mod, names(pars), as_list = T)
east_fit <- simulate_SAMint_detC_logpi_ss(east, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_east_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(east_fit, fit_pars)
dev.off()

# newR
ANT <- as.matrix(calc_antecedent_drivers(newR$GPP, 4, 2))
sim_new <- simulate_SAMint_detC_logpi_ss(newR, pars, ANT, antdays)
obs_error = median(-newR$ER.lower + newR$ER.upper)/4
obs_error_sd = sd(-newR$ER.lower + newR$ER.upper)/4
sim_new %>% select(date, GPP, Q, R, ER,C, litter )%>%
    pivot_longer( -date, names_to = 'var', values_to = 'value') %>%
    ggplot(aes(date, value, col = var)) +
    geom_line()+
    facet_wrap(.~var, scales = 'free_y', ncol = 2)

stan_dat <- list(ndays = nrow(sim_new), nweights = nweights, antdays = antdays,
                 R_obs = sim_new$ER, P = sim_new$GPP, C0 = 100,
                 temp = sim_new$temp_C, litter = sim_new$litter,
                 tau = sim_new$tau, ANT = ANT, sigma_obs_prior = obs_error,
                 sigma_obs_sd = obs_error_sd)
mod <- stan('src/SAM/stan/SAMint_detC_logpi_fixedoi_ss.stan',
            data = stan_dat,
            chains = 4,  cores = 4,
            iter = 2000, control = list(max_treedepth = 12))
saveRDS(mod, 'src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_new.rds')
mod <- readRDS('src/SAM/stan/fits/SAMint_detC_logpi_fixedoi_ss_new.rds')
beep(sound = 8)

a <- traceplot(mod, ncol = 2, names(pars))
b <- plot(mod, pars = names(pars))
# c <- pairs(mod, pars = c('K_20', 'tau0', 'beta_p', 'sigma_proc',
#                          'sigma_obs', 'lp__'))

png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_new.png',
    height = 800, width = 400)
ggpubr::ggarrange(a, b, ncol = 1)
dev.off()


fit_pars <- get_pars(mod, names(pars), as_list = T)
new_fit <- simulate_SAMint_detC_logpi_ss(newR, fit_pars, ANT, antdays)
png('figures/simulation_fits/rivers/SAMint_detC_logpi_fixedoi_ss_new_fit.png',
    width = 600, height = 400)
plot_ER_breakdown(new_fit, fit_pars)
dev.off()
