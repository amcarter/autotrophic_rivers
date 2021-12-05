# Stan Helper functions
#A Carter


# plotting ####
plot_post_sim <- function(fit, pars, vals){
  dd <- data.frame(x = vals, y = length(vals):1)
  p <- rstan::plot(fit, show_density = T, fill_color = 'grey',
                   pars = pars) + 
    geom_point(data = dd, aes(x = x, y = y), size = 3, 
               shape = 17, col = 'brown3')
  return(p)
}

