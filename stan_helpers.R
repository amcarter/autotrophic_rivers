# Stan Helper functions
#A Carter


# plotting ####
plot_post_sim <- function(fit, pars, vals, xlim = NULL){
  dd <- data.frame(x = vals, y = length(vals):1)
  p <- rstan::plot(fit, show_density = T, fill_color = 'grey',
                   pars = pars) +
    geom_point(data = dd, aes(x = x, y = y), size = 3,
               shape = 17, col = 'brown3')+
    xlim(xlim)
  return(p)
}

# distributions ####
est_beta_params <- function(mu, sigma) {
  alpha <- ((1 - mu) / sigma^2 - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

