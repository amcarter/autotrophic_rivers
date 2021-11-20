
  data {
    int <lower = 0> N; // Sample size
    vector[N] logR;    // log ER observed
    real mu_obs;       // prior mean for obs error
  }
  
  parameters {
    real a0; // Intercept
    real <lower = 0, upper = 1> a1; 
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
    vector <lower = 0> [N] mu; // State variable for ER
  }
  
  model {
    mu[1] ~ normal(exp(logR[1]), 0.1);
    mu[2:N] ~ normal(a0 + a1 * mu[1:(N-1)], sigma_proc);
    logR ~ normal(log(mu), sigma_obs);
   
    // priors
    a0 ~ normal(0,1);
    a1 ~ uniform(0,1);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
  }
  
  generated quantities {
    vector [N] R_hat;   // estimated underlying state
    vector [N] R_tilde; // estimated log ER obs
    R_hat[1] = exp(logR[1]);
    R_tilde[1] = logR[1];
    for(i in 2:N){
      R_hat[i] = normal_rng(a0 + a1 * R_hat[i-1], sigma_proc);
      R_tilde[i] = normal_rng(log(R_hat[i]), sigma_obs);  
    }
  }
