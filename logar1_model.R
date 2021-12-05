library(rstan)
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
sink('src/SAM/stan/ar1_model_logR.stan')
cat("
  data {
    int <lower = 0> N;  // Sample size
    vector[N] logR_obs; // log ER observed
    real mu_obs;        // prior mean for obs error
  }
  
  parameters {
    real b0; // Intercept
    real <lower = 0, upper = 1> phi; 
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
    vector <lower = 0> [N] logR; // State variable for ER
  }
  
  model {
    logR[1] ~ normal(logR_obs[1], sigma_obs);
    logR[2:N] ~ normal(b0 + phi * logR[1:(N-1)], sigma_proc);
    logR_obs ~ normal(logR, sigma_obs);
   
    // priors
    b0 ~ normal(0,1);
    phi ~ beta(1,1);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
  }
  
  generated quantities {
    vector [N] R_hat;   // estimated underlying state
    vector [N] R_tilde; // estimated log ER obs
    R_hat[1] = logR[1];
    R_tilde[1] = logR[1];
    for(i in 2:N){
      R_hat[i] = normal_rng(b0 + phi * R_hat[i-1], sigma_proc);
      R_tilde[i] = normal_rng(R_hat[i], sigma_obs);  
    }
  }", fill = T)


sink()

# simulate data
b0 = 1      # intercept
phi = 0.38  # ar1 coefficient
sigma_proc = 0.1 # ~10% error
sigma_obs = 0.05  

logR <- numeric(100)  # ER state variable
logR[1] <- log(5)

# Observed ER, on a log scale
logR_obs = rnorm(1, logR[1], sigma_obs)
for (i in 2:100){
  logR[i] <- b0 + phi * logR[i-1] + rnorm(1, 0, sigma_proc)
  logR_obs[i] = rnorm(1, logR[i], sigma_obs)
}

plot(exp(logR), type = 'l')
points(exp(logR_obs))

sim_dat <- list(logR_obs = logR_obs, N = length(logR_obs), mu_obs = sigma_obs)
fit <- stan(file = 'src/SAM/stan/ar1_model_logR.stan', data = sim_dat,
            warmup = 500, iter = 5000, 
            chains = 4, cores = 4)

saveRDS(fit, 'src/ar1/log_ar1_sim_fit.rds')
