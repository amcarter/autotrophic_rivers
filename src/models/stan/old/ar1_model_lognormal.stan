
   data {
    int <lower = 0> N; // Sample size
    vector[N] R_obs;
    real mu_obs;
  }
  
  parameters {
    real a0; // Intercept
    real <lower = 0, upper = 1> a1; //phi
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
    vector <lower = 0> [N] R;
  
  }
  
  model {
    R[1] ~ lognormal(log(R_obs[1]), 0.05);
    for(t in 2:N){
      R[t] ~ lognormal(log(a0 + a1 * R[t-1]), sigma_proc);
      R_obs[t] ~ lognormal(log(R[t]), sigma_obs);
    }
   
    a0 ~ normal(0,1);
    a1 ~ uniform(0,1);
    sigma_proc ~ normal(0,0.2) T[0,];
    sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
  }
  
  generated quantities {
    
  }
