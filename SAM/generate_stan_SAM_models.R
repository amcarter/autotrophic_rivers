# Stochastic Antecedent Model based on Ogle et al 2015
# base model is from Bob Hall
# 10/25/2021
setwd("C:/Users/Alice Carter/git/autotrophic_rivers")
source('src/stan_helpers.R')
library(rstan)
library(shinystan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TURE)
#
# Basic SAM code ####
# ER as a function of GPP on days t-4:t

sink("src/SAM/stan/SAM.stan")

cat("
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real a0;
      real <lower =0> a1;
      simplex [nweight] w; 
      real <lower=0> sigma;
    }
    
    transformed parameters{
      vector  [N] Pant;
      Pant[1:5] = P[1:5];

      for (i in 6:N){
        vector  [nweight] Pvec;
        for(j in 1:nweight){ 
          Pvec[j]=w[j]*P[i-(j-1)];
        }
        Pant[i]=sum(Pvec);
      }
    }
    
    model {
      for (i in 6:N){
        R[i] ~ normal(a0 + a1*Pant[i], sigma); // likelihood
      }
      

    a0 ~ normal(0,5); //priors
    a1 ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma ~ normal(0,1) T[0,];
    }
    
    generated quantities{
      vector [N] R_hat;
      R_hat[1:5] = R[1:5];
      for(i in 6:N){
        R_hat[i] = normal_rng(a0 + a1 * Pant[i], sigma);
      }
      
    }" ,fill=TRUE)

sink()

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
saveRDS(list(fit = fit_fake, dat = sim_dat), 'src/SAM/stan/fits/fit_fake.rds')


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



# Basic autoregressive model ####

sink('src/SAM/stan/ar1_model.stan')
cat("
  data {
    int <lower = 0> N; // Sample size
    vector[N] R_obs;   // time series of observations
  }
  
  parameters {
    real a0; // Intercept
    real <lower = 0, upper = 1> phi;//ar1 coeff
    real <lower = 0> sigma_proc;    //process error
    real <lower = 0> sigma_obs;     //observation error
    vector <lower = 0> [N] R;       //state vector
  }
  
  model {
    // initial state
    R[1] ~ normal(R_obs[1], sigma_proc);
    
    // process model
    R[2:N] ~ normal(a0 + phi * R[1:(N-1)], sigma_proc);
    
    //observation model
    R_obs ~ normal(R, sigma_obs);
   
    //priors
    a0 ~ normal(0,1);
    phi ~ beta(1,1);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(0,1) T[0,];
  }
  
  generated quantities {
    vector [N] R_hat;
    vector [N] R_tilde;
    R_hat[1] = R_obs[1];
    R_tilde[1] = normal_rng(R_obs[1], sigma_obs);
    for(i in 2:N){
      R_hat[i] = normal_rng(a0 + phi * R_hat[i-1], sigma_proc);
      R_tilde[i] = normal_rng(R_hat[i], sigma_obs);  
    }
  }", fill = T)

sink()

# simulate data linear model ####

a0 <- 2
phi <- 0.8
sigma_proc <- 1
sigma_obs <- 0.3

R <- numeric(100)
R_obs <- numeric(100)
R[1] <- 10
R_obs[1] =  rnorm(1, R[1], sigma_obs)
for (i in 2:100){
  R[i] <- rnorm(1, a0 + phi * R[i-1], sigma_proc)
  R_obs[i] = rnorm(1, R[i], sigma_obs)
}


# Plot simulated data####
plot(R, type = 'l', xlab = 'time')
points(R_obs)
legend('topright', c('state', 'observations'), lty = c(1,0), pch = c(NA, 1), 
       ncol = 2, bty = 'n')

sim_dat <- list(R_obs = R_obs, N = length(R_obs))#, mu_obs = sigma_obs)
fit <- stan(file = 'src/SAM/stan/ar1_model.stan', data = sim_dat,
            warmup = 500, iter = 1000, 
            chains = 4, cores = 4)


# saveRDS(list(fit = fit, dat = sim_dat), 'src/SAM/stan/fits/simulated_ar1_fit.rds')
traceplot(fit,pars = c('a0', 'phi', 'sigma_proc', 'sigma_obs'))
# print(fit, pars = c('a0', 'phi', 'sigma_proc', 'sigma_obs'))
plot_post_sim(fit, pars = c('a0', 'phi', 'sigma_proc', 'sigma_obs'),
              vals = c(a0, phi, sigma_proc, sigma_obs))
pairs(fit, pars = c('a0', 'phi', 'sigma_proc', 'sigma_obs', 'lp__'))

y_rep <- as.matrix(fit, pars = 'R_tilde')
bayesplot::ppc_dens_overlay(R_obs, y_rep[1:200,])


launch_shinystan(fit)
# simulate data log model####

a0 <- 2
a1 <- 0.8
sigma_proc <- 0.05
sigma_obs <- 0.05

R <- numeric(100)
R_obs <- numeric(100)
R[1] <- 10
R_obs[1] =  rlnorm(1, log(R[1]), sigma_obs)
for (i in 2:100){
  R[i] <- rlnorm(1, log(a0 + a1 * R[i-1]), sigma_proc)
  R_obs[i] = rlnorm(1, log(R[i]), sigma_obs)
}

# Combined SAM and autoregressive model ####


sink("src/SAM/stan/SAM_ar1_lognormal.stan")

cat("
    data {
      int <lower = 0> N; // Sample size
      vector[N] R_obs;
      vector [N] P;
      real mu_obs;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
     real a0; // Intercept
     real <lower =0> a1;
     real < lower = 0, upper = 1> phi; 
     real < lower = 0 > sigma_obs;
     real < lower = 0 > sigma_proc;
     simplex [nweight] w; 
     vector [N] R;
    }
    
    transformed parameters{
      vector  [N] Pant;
      Pant[1:5] = P[1:5];
  
      for (i in 6:N){
        vector  [nweight] Pvec;
        for(j in 1:nweight){ 
          Pvec[j]=w[j]*P[i-(j-1)];
        }
        Pant[i]=sum(Pvec);
      }
    }    
    
    model {
      R[1:5] ~ lognormal(log(R_obs[1:5]), 0.05);
      for(t in 6:N){
        R[t] ~ lognormal(log(a0 + a1 * Pant[t] + phi * R[t-1]), sigma_proc);
        R_obs[t] ~ lognormal(log(R[t]), sigma_obs);
      }
      
      a0 ~ normal(0,1);
      a1 ~ normal(0,1) T[0,];
      phi ~ uniform(0,1);
      w ~ dirichlet(alpha);
      sigma_proc ~ normal(0,0.2) T[0,];
      sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
    }

  generated quantities {
    vector [N] R_hat;
    vector [N] R_tilde;
    vector [N] Pant_hat;
    
    Pant_hat[1:5] = P[1:5];
    for(i in 6:N){
      vector [nweight] Pvec_hat;
      for(j in 1:nweight){
        Pvec_hat[j] = w[j]*P[i-(j-1)];
      }
      Pant_hat[i] = sum(Pvec_hat);
    }
    
    R_hat[1:5] = R[1:5];
    for(i in 6:N){
      R_hat[i] = lognormal_rng(log(a0 + a1 * Pant_hat[i] + phi * R_hat[i-1]), sigma_proc);
    }
    for(i in 1:N){
      R_tilde[i] = lognormal_rng(log(R_hat[i]), sigma_obs);
    }    
  }", fill = T)
sink()

# Fake Data: Assign weight with 50% each on lags day 1 and 2

# set parameters
a0 = 0
a1 = 0.3
phi = 0.7
w <- c(0.5,0.5,0,0,0)
sigma_proc = 0.05
sigma_obs = 0.05

# generate data
P <- numeric(300)
R <- numeric(300)
P[1] <- 10
for (i in 2:300)
  P[i] <- rlnorm(1, log(1 + 0.9 * P[i-1]), 0.05)

Pant <- numeric(300)
Pant[1:5]<-P[1:5]

for (i in 6:300){
  Pvec<-numeric(5)
  for(j in 1:5){
    Pvec[j]<-w[j]*P[i-(j-1)]
  }
  Pant[i]<-sum(Pvec)
  
}

R <- a0 + a1 * P[1] + phi * 5
R_obs = rlnorm(1, log(R), sigma_obs)
for (i in 2:300){
  R[i] = rlnorm(1, log(a0 + a1 * Pant[i] + phi * R[i-1]), sigma_proc)
  R_obs[i] = rlnorm(1, log(R[i]), sigma_obs)
}
plot(R, type = 'l', ylim = c(6,15))
points(R_obs)
lines(P, lty = 2)

## Run SAM model
sim_dat <- list(R_obs = R_obs, P = P, N = length(P), mu_obs = 0.05, 
                nweight = 5, alpha = rep(1, 5))
fit_fake <- stan(file = 'src/SAM/stan/SAM_ar1_lognormal.stan',
                 data = sim_dat,
                 warmup = 500, iter = 1000,
                 chains = 1, cores = 4)

saveRDS(list(fit = fit_fake, dat = sim_dat), 'src/SAM/stan/fits/simulated_sam_ar1.rds')

print(fit_fake, pars=c("a0", "a1", "a2", "w", 'sigma_proc', 'sigma_obs'))

# plot estimates vs true values
p <- rstan::plot(fit, show_density = T, fill_color = 'grey',
                 pars = c('a0', 'a1', 'sigma_obs', 'sigma_proc'))
dd <- data.frame(x = c(a0, a1, sigma_obs, sigma_proc), y = c(4,3,2,1))
p + geom_point(data = dd, aes(x = x, y = y), size = 3, shape = 17, col = 'brown3')

pairs(fit, pars = c('a0', 'a1', 'sigma_obs', 'sigma_proc'))


fit_fake <- readRDS('src/SAM/stan/fits/simulated_sam_ar1.rds')
post <- extract(fit_fake)
plot(post$sigma_obs, post$sigma_proc)
# SAM model for different previous intervals of GPP ####
#   ER_t = f(GPP_t, GPP_t-1, GPP_week, GPP_month)

sink("src/SAM/stan/SAM_intervals.stan")

cat("
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real a0;
      real <lower =0> a1;
      simplex [nweight] w; 
      real <lower=0> sigma;
    }
    
    transformed parameters{
      vector  [N] Pant;
      Pant[1:30] = P[1:30];

      for (i in 31:N){
        real Pweek;
        real Pmonth;
        Pweek = sum(P[(i-7):(i-1)])/7;
        Pmonth = sum(P[(i-30):(i-1)])/30;
        
        Pant[i] = w[1]*P[i] + w[2]*P[i-1] + w[3]*Pweek + w[4]*Pmonth;
        
      }
    }
    
    model {
      for (i in 31:N){
        R[i] ~ normal(a0 + a1 * Pant[i], sigma); // likelihood
      }
      
    a0 ~ normal(0,5); //priors
    a1 ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    
    }" ,fill=TRUE)

sink()
