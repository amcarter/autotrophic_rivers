# adapted from bob's Stan_report_2.rmd

library(rstan)
library(shinystan)
setwd("C:/Users/Alice Carter/git/autotrophic_rivers")
source('src/SAM/generate_stan_SAM_models.R')
# Basic SAM implementation ####
#   ER = a0 + a1*f(GPP[t-4:t]) + epsilon 

# Fake Data: Assign weight with 50% each on lags day 1 and 2

# set parameters
a0 = 0
a1 = 0.5
w <- c(0.5,0.5,0,0,0)
sigma = 0.2

# generate data
P <- numeric(100)
R <- numeric(100)
P[1] <- 10
for (i in 2:100)
  P[i] <- 1+0.9* P[i-1]+rnorm(1,0,0.85)

Pant <- numeric(100)
Pant[1:5]<-P[1:5]

for (i in 6:100){
  Pvec<-numeric(5)
  for(j in 1:5){
    Pvec[j]<-w[j]*P[i-(j-1)]
  }
  Pant[i]<-sum(Pvec)
  
}

R <- a0 + a1 * P[1]
for (i in 2:100){
  R[i] <- a0 + a1*Pant[i] + rnorm(1, 0, sigma)
}

plot(P,R)

## Run SAM model
sim_dat <- list(R = R, P = P, N = length(P), nweight = 5, alpha = rep(1, 5))
fit_fake <- stan(file = 'src/SAM/stan/SAM.stan', 
                 data = sim_dat, 
                 warmup = 500, iter = 1000, 
                 chains = 4, cores = 4)
print(fit_fake, pars=c("a0", "a1", "w", 'sigma'))
plot(fit_fake, pars = c("a0", "a1", "w", 'sigma'))
# launch_shinystan(fit_fake)
saveRDS(fit_fake, 'src/SAM/stan/fits/fit_fake.rds')


## Predict data
post <- summary(fit_fake, pars = c('a0', 'a1', 'w', 'sigma'))$summary
R_hat <- summary(fit_fake, pars = 'R_hat')$summary %>%
  as_tibble() %>%
  pull(mean)

# Do predictions correlate with data?
plot(R, R_hat)
abline(0,1)

# are we accounting for all of the autocorrelation in the respiration?
par(mfrow = c(1,2))
acf(R, main = 'Autocorrelation of ER')
acf(R - R_hat, main = 'Autocorrelation of ER residuals')

  
# SAM model for different previous intervals of GPP ####
## Fake Data: Assign weight with 50% each on lags day 1 and 2

# set parameters
a0 = 0
a1 = 0.5
w <- c(0.4, 0.3, 0.3, 0)
sigma = 0.2

# generate data
P <- numeric(365)
R <- numeric(365)
P[1] <- 10
for (i in 2:365)
  P[i] <- 1+0.9* P[i-1]+rnorm(1,0,0.85)

Pant <- numeric(365)
Pant[1:30]<-P[1:30]

for (i in 31:365){
  Pweek = sum(P[(i-7):(i-1)])/7
  Pmonth = sum(P[(i-30):(i-1)])/30
  Pant[i] = w[1]*P[i] + w[2]*P[i-1] + w[3]*Pweek + w[4]*Pmonth
}


R <- a0 + a1 * P[1]
for (i in 2:365){
  R[i] <- a0 + a1 * Pant[i] + rnorm(1, 0, sigma)
}

plot(P,R)

pacf(R)

## How did we do?
sim_dat <- list(R = R, P = P, N = length(P), nweight = 4, alpha = rep(1, 4))
fit_fake2 <- stan(file = 'src/SAM/stan/SAM_intervals.stan', 
                  data = sim_dat, 
                  warmup = 500, iter = 1000, 
                  chains = 4, cores = 4)
print(fit_fake2, pars=c("a0", "a1", "w", 'sigma'))
plot(fit_fake2, pars = c("a0", "a1", "w", 'sigma'))

saveRDS(fit_fake2, 'src/SAM/stan/fits/fit_fake2.rds')

## Predict data
p_est <- summary(fit_fake2, pars = c('a0', 'a1', 'w', 'sigma'), 
                 probs = c(0.025, 0.975))$summary

a0 <- p_est[1,1]
a1 <- p_est[2,1]
w <- p_est[3:6,1]

Pant_pred = numeric(365)
Pant_pred[1:30] <- P[1:30]

for (i in 31:365){
  Pweek = sum(P[(i-7):(i-1)])/7
  Pmonth = sum(P[(i-30):(i-1)])/30
  Pant_pred[i] = w[1]*P[i] + w[2]*P[i-1] + w[3]*Pweek + w[4]*Pmonth
}


r_pred <- numeric(365)
for(i in 1:365){
  r_pred[i] = a0 + a1 * Pant_pred[i]}

# Do predictions correlate with data?
plot(R, r_pred)
abline(0,1)

# are we accounting for all of the autocorrelation in the respiration?
par(mfrow = c(2,1))
acf(R, main = 'Autocorrelation of ER')
acf(R - r_pred, main = 'Autocorrelation of ER residuals')


# in both cases we have fully accounted for autocorrelation in respiration
# Next: Figure out how to include posterior predictive estimates in the model code

