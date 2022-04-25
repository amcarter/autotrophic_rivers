# Stan Helper functions
#A Carter


# plotting ####
plot_post_sim <- function(fit, pars, vals, title = NULL, xlim = NULL){
  dd <- data.frame(x = vals, y = length(vals):1)
  p <- rstan::plot(fit, show_density = T, fill_color = 'grey',
                   pars = pars) +
    geom_point(data = dd, aes(x = x, y = y), size = 3,
               shape = 17, col = 'brown3')+
    xlim(xlim)+
    ggtitle(title)
  return(p)
}

# post pred checks
plot_ER_ppcheck <- function(mod, ss, figname){
    fit <- summary(mod)
    ndays <- nrow(ss)
    # extract parameters
    K20_mean <- fit$summary[4, 1]
    K20_sd <- fit$summary[4, 3]
    beta_s_mean <- fit$summary[1,1]
    beta_s_sd <- fit$summary[1,3]
    sigma_obs_mean = fit$summary[3,1]
    sigma_obs_sd <- fit$summary[3,3]
    sigma_proc_mean <- fit$summary[2,1]
    sigma_proc_sd <- fit$summary[2,3]

    # predict data
    ss$K_pred = calc_rate_coef(ss$temp_C, K_20 = K20_mean/100)
    ss$K_pred.lower = calc_rate_coef(ss$temp_C, K_20 = (K20_mean- K20_sd)/100)
    ss$K_pred.upper = calc_rate_coef(ss$temp_C, K_20 = (K20_mean+ K20_sd)/100)


    Chat = Chat.lower = Chat.upper = numeric()
    Chat[1] = Chat.lower[1] = Chat.upper[1] = C0
    ss$C_pred <- ss$C_pred.lower <- ss$C_pred.upper <- NA_real_
    ss$R_pred <- ss$R_pred.lower <- ss$R_pred.upper <- NA_real_
    ss$C_pred[1] = ss$C_pred.lower = ss$C_pred.upper =
        exp(rnorm(1, log(C0), sigma_proc_mean/10))
    ss$R_pred[1] = ss$AR[1] - ss$K_pred[1] * ss$C_pred[1]
    ss$R_pred.lower[1] = ss$AR[1] - ss$K_pred.lower[1] * ss$C_pred[1]
    ss$R_pred.upper[1] = ss$AR[1] - ss$K_pred.upper[1] * ss$C_pred[1]

    for(i in 2:ndays){
        Chat[i] = (Chat[i-1] + ss$litter[i] + ss$R_pred[i-1] + ss$GPP[i-1])*
            (1-beta_s_mean*ss$Q[i])
        Chat.lower[i] = (ss$C_pred.lower[i-1] + ss$litter[i] + ss$R_pred.lower[i-1] +
                             ss$GPP[i-1])*(1-(beta_s_mean + beta_s_sd)*ss$Q[i])
        Chat.upper[i] = (ss$C_pred.upper[i-1] + ss$litter[i] + ss$R_pred.upper[i-1] +
                             ss$GPP[i-1])*(1-(beta_s_mean- beta_s_sd)*ss$Q[i])
        ss$C_pred[i] = exp(rnorm(1, log(Chat[i]), sigma_proc_mean/10))
        ss$C_pred.lower[i] = exp(rnorm(1, log(Chat.lower[i]), sigma_proc_mean/10))
        ss$C_pred.upper[i] = exp(rnorm(1, log(Chat.upper[i]), sigma_proc_mean/10))
        ss$R_pred[i] = ss$AR[i] - ss$K_pred[i]*ss$C_pred[i]
        ss$R_pred.upper[i] = ss$AR[i] - ss$K_pred.upper[i]*ss$C_pred.upper[i]
        ss$R_pred.lower[i] = ss$AR[i] - ss$K_pred.lower[i]*ss$C_pred.lower[i]
    }

    # observation model:
    ss$R_mod = rnorm(ndays, ss$R_pred, sigma_obs_mean)
    ss$R_mod.upper = rnorm(ndays, ss$R_pred.upper, sigma_obs_mean)
    ss$R_mod.lower = rnorm(ndays, ss$R_pred.lower, sigma_obs_mean)

    ss$C_modeled <- fit$summary[5:(ndays+4),1]

    ests <- ss %>%
        select(date, C_modeled, ER_modeled = R_mod,
               ER_measured = ER, GPP_measured = GPP) %>%
        mutate(GPP_modeled = NA_real_, C_measured = NA_real_)
    ests %>%
        pivot_longer(ends_with(c("_measured", "_modeled")), names_to = c('variable', 'est'),
                     names_sep = '_', values_to = 'value') %>%
        ggplot(aes(date, value, col = est)) +
        geom_line() +
        facet_wrap(.~variable, scales = 'free', ncol = 1)

    png(figname)
        par(mfrow = c(2, 1), mar = c(2, 4, 1, 2))
        plot(ss$date, ss$ER, type = 'l')
        lines(ss$date, ss$R_mod, col = 'steelblue')
        polygon(c(ss$date, rev(ss$date)), c(ss$R_mod.lower, rev(ss$R_mod.upper)),
                col = alpha('steelblue',0.5), border = NA)
        legend('bottomright', c('ER', 'modeled ER'), col = c(1, 'steelblue'),
               lty = c(1,1), bty = 'n')

        plot(ss$date, ss$R_mod, type = 'l')
        lines(ss$date, ss$AR, col = 'forestgreen')
        legend('bottomright', c('total R', 'auto R'), col = c(1, 'forestgreen'),
               lty = c(1,1), bty = 'n')
    dev.off()
}
plot_ER_ppcheck_SAM5 <- function(mod, ss, figname){
    fit <- summary(mod)
    ndays <- nrow(ss)
    # extract parameters
    K20_mean <- fit$summary[ndays +10, 1]
    K20_sd <- fit$summary[ndays+10, 3]
    beta_s_mean <- fit$summary[1,1]
    beta_p_mean <- fit$summary[2,1]
    w_mean <- fit$summary[3:7,1]
    beta_s_sd <- fit$summary[1,3]
    sigma_obs_mean = fit$summary[9,1]
    sigma_obs_sd <- fit$summary[3,3]
    sigma_proc_mean <- fit$summary[8,1]
    sigma_proc_sd <- fit$summary[2,3]

    # predict data
    ss$K_pred = calc_rate_coef(ss$temp_C, K_20 = K20_mean/100)
    ss$K_pred.lower = calc_rate_coef(ss$temp_C, K_20 = (K20_mean- K20_sd)/100)
    ss$K_pred.upper = calc_rate_coef(ss$temp_C, K_20 = (K20_mean+ K20_sd)/100)


    Chat = Chat.lower = Chat.upper = numeric()
    Chat[1] = Chat.lower[1] = Chat.upper[1] = C0
    ss$C_pred <- ss$C_pred.lower <- ss$C_pred.upper <- NA_real_
    ss$R_pred <- ss$R_pred.lower <- ss$R_pred.upper <- NA_real_
    ss$HR_det <- ss$HR_alg <- NA_real_
    ss$C_pred[1] = ss$C_pred.lower = ss$C_pred.upper =
        exp(rnorm(1, log(C0), sigma_proc_mean/10))
    ss$HR_det[1] =  - ss$K_pred[1] * ss$C_pred[1]
    ss$R_pred.upper[1] = ss$AR[1] - ss$K_pred.upper[1] * ss$C_pred[1]

    for(i in 2:ndays){
        Chat[i] = (Chat[i-1] + ss$litter[i] + ss$HR_det[i-1])*
            (1-beta_s_mean*ss$Q[i])
        Chat.lower[i] = (ss$C_pred.lower[i-1] + ss$litter[i] + ss$R_pred.lower[i-1] +
                             ss$GPP[i-1])*(1-(beta_s_mean + beta_s_sd)*ss$Q[i])
        Chat.upper[i] = (ss$C_pred.upper[i-1] + ss$litter[i] + ss$R_pred.upper[i-1] +
                             ss$GPP[i-1])*(1-(beta_s_mean- beta_s_sd)*ss$Q[i])
        ss$C_pred[i] = exp(rnorm(1, log(Chat[i]), sigma_proc_mean/10))
        ss$C_pred.lower[i] = exp(rnorm(1, log(Chat.lower[i]), sigma_proc_mean/10))
        ss$C_pred.upper[i] = exp(rnorm(1, log(Chat.upper[i]), sigma_proc_mean/10))
        ss$HR_det[i] = - ss$K_pred[i]*ss$C_pred[i]
        ss$R_pred.upper[i] = ss$AR[i] - ss$K_pred.upper[i]*ss$C_pred.upper[i]
        ss$R_pred.lower[i] = ss$AR[i] - ss$K_pred.lower[i]*ss$C_pred.lower[i]
    }

    ss$AR = -ARf * ss$GPP
    ss$HR_alg = -beta_p_mean * ss$Pant
    ss$R_pred = ss$AR + ss$HR_det + ss$HR_alg

    # observation model:
    ss$R_mod = rnorm(ndays, ss$R_pred, sigma_obs_mean)
    ss$R_mod.upper = rnorm(ndays, ss$R_pred.upper, sigma_obs_mean)
    ss$R_mod.lower = rnorm(ndays, ss$R_pred.lower, sigma_obs_mean)

    ss$C_modeled <- fit$summary[5:(ndays+4),1]

    ests <- ss %>%
        select(date, C_modeled, ER_modeled = R_mod,
               ER_measured = ER, GPP_measured = GPP) %>%
        mutate(GPP_modeled = NA_real_, C_measured = NA_real_)
    ests %>%
        pivot_longer(ends_with(c("_measured", "_modeled")), names_to = c('variable', 'est'),
                     names_sep = '_', values_to = 'value') %>%
        ggplot(aes(date, value, col = est)) +
        geom_line() +
        facet_wrap(.~variable, scales = 'free', ncol = 1)

    png(figname)
        par(mfrow = c(2, 1), mar = c(2, 4, 1, 2))
        plot(ss$date, ss$ER, type = 'l')
        lines(ss$date, ss$R_mod, col = 'steelblue')
        legend('bottomright', c('ER', 'modeled ER'), col = c(1, 'steelblue'),
               lty = c(1,1), bty = 'n', ncol = 2)

        plot(ss$date, ss$R_mod, type = 'l')
        polygon(c(ss$date, rev(ss$date)), c(ss$AR, rep(0, ndays)),
                col = alpha('goldenrod', .5), border = NA)
        polygon(c(ss$date, rev(ss$date)), c(ss$AR, rev(ss$AR + ss$HR_alg)),
                col = alpha('forestgreen', .5), border = NA)
        polygon(c(ss$date, rev(ss$date)), c(ss$AR + ss$HR_alg, rev(ss$R_mod)),
                col = alpha('sienna', .5), border = NA)
        legend('bottomright',
               c('autotrophic R', 'hetero algal R', 'hetero detrital R'),
                fill = c('goldenrod', 'forestgreen', 'sienna'), border = NA,
                bty = 'n', ncol = 3)
    dev.off()
}

# distributions ####
est_beta_params <- function(mu, sigma) {
  alpha <- ((1 - mu) / sigma^2 - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

