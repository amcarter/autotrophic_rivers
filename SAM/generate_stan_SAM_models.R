# Stochastic Antecedent Model based on Ogle et al 2015
# base model is from Bob Hall
# 10/25/2021
setwd("C:/Users/Alice Carter/git/autotrophic_rivers")
library(rstan)
library(shinystan)
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



# basic autoregressive model ####

sink('src/SAM/stan/ar1_model.stan')
cat("
  data {
    int <lower = 0> N; // Sample size
    vector[N] logR;
    real mu_obs;
  }
  
  parameters {
    real a0; // Intercept
    real <lower = 0, upper = 1> a1; 
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
    vector <lower = 0> [N] mu;
  }
  
  model {
    mu[1] ~ normal(exp(logR[1]), 0.1);
    mu[2:N] ~ normal(a0 + a1 * mu[1:(N-1)], sigma_proc);
    logR ~ normal(log(mu), sigma_obs);
   
   
    a0 ~ normal(0,1);
    a1 ~ uniform(0,1);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
  }
  
  generated quantities {
    vector [N] R_hat;
    vector [N] R_tilde;
    R_hat[1] = exp(logR[1]);
    R_tilde[1] = logR[1];
    for(i in 2:N){
      R_hat[i] = normal_rng(a0 + a1 * R_hat[i-1], sigma_proc);
      R_tilde[i] = normal_rng(log(R_hat[i]), sigma_obs);  
    }
  }", fill = T)


sink()

# simulate data
a0 = 1
a1 = 0.8
sigma_proc = 0.2
sigma_obs = 0.1

mu <- numeric(100)
mu[1] <- 5
R = log(mu[1]) + rnorm(1, 0, sigma_obs)
for (i in 2:100){
  mu[i] <- a0 + a1 * mu[i-1] + rnorm(1, 0, sigma_proc)
  R[i] = log(mu[i]) + rnorm(1, 0, 0.2)
}

plot(mu, type = 'l')
points(exp(R))

plot(mu, abs(exp(R) - mu))

sim_dat <- list(logR = R, N = length(R), mu_obs = sigma_obs)
fit <- stan(file = 'src/SAM/stan/ar1_model.stan', data = sim_dat,
            warmup = 500, iter = 5000, 
            chains = 4, cores = 4)


saveRDS(list(fit = fit, dat = sim_dat), 'src/SAM/stan/fits/simulated_ar1_fit.rds')
print(fit, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs'))
plot(fit, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs'))

pairs(fit, pars = c('a0', 'a1', 'sigma_proc', 'sigma_obs'))

launch_shinystan(fit)

# Combined SAM and autoregressive model ####


sink("src/SAM/stan/SAM_ar1.stan")

cat("

  data {
   int <lower = 0> N;
   vector[N] R;
   vector [N] P;
   int <lower = 0> nweight;
   vector [nweight] alpha;
   real mu_obs;
   real tau_obs;
  
  }
  
  parameters {
   real a0; // Intercept
   real < lower = 0, upper = 1> a1; 
   real <lower =0> a2;
   real < lower = 0 > sigma_obs;
   real < lower = 0 > sigma_proc;
   simplex [nweight] w; 
   vector [N] mu;
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
   mu[1:5] ~ normal(R[1:5], 0.01);
   for (i in 6:N){
    mu[i] ~ normal(a0 + a1 * mu[i-1] + a2*Pant[i], sigma_proc);
    R[i] ~ normal(mu[i], sigma_obs);
   }
  
   a0 ~ normal(0,5); //priors
   a1 ~ uniform(0,1);
   a2 ~ normal(0,1) T[0,];
   w ~ dirichlet(alpha);
   sigma_obs ~ normal(mu_obs,tau_obs) T[0,];
   sigma_proc ~ normal(0,1) T[0,];
  }
  
  generated quantities {
    vector [N] R_hat;
    R_hat[1:5] = R[1:5];
    for(i in 6:N){
      R_hat[i] = normal_rng(a0 + a1 * R_hat[i-1] + a2*Pant[i], sigma_proc);
    }
    
  }", fill = T)
sink()

# Fake Data: Assign weight with 50% each on lags day 1 and 2

# set parameters
a0 = 0
a1 = 0.7
a2 = 0.3
w <- c(0.5,0.5,0,0,0)
sigma_proc = 0.5
sigma_obs = 0.2

# generate data
P <- numeric(300)
R <- numeric(300)
P[1] <- 10
for (i in 2:300)
  P[i] <- 1+0.9* P[i-1]+rnorm(1,0,0.85)

Pant <- numeric(300)
Pant[1:5]<-P[1:5]

for (i in 6:300){
  Pvec<-numeric(5)
  for(j in 1:5){
    Pvec[j]<-w[j]*P[i-(j-1)]
  }
  Pant[i]<-sum(Pvec)
  
}

mu <- a0 + a1 * 5 + a2 *P[1]
R = rnorm(1, mu, sigma_obs)
for (i in 2:300){
  mu[i] = rnorm(1, a0 + a1*mu[i-1] + a2*Pant[i], sigma_proc)
  R[i] = rnorm(1, mu[i], sigma_obs)
}
plot(mu, type = 'l', ylim = c(2,12))
points(R)
lines(P, lty = 2)

## Run SAM model
sim_dat <- list(R = R, P = P, N = length(P), nweight = 5, alpha = rep(1, 5))
fit_fake <- stan(file = 'src/SAM/stan/SAM_ar1.stan',
                 data = sim_dat,
                 warmup = 500, iter = 1000,
                 chains = 4, cores = 4)
print(fit_fake, pars=c("a0", "a1", "a2", "w", 'sigma_proc', 'sigma_obs'))
plot(fit_fake)

saveRDS(list(fit = fit_fake, dat = sim_dat), 'src/SAM/stan/fits/simulated_sam_ar1.rds')

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
