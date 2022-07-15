
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

  }
