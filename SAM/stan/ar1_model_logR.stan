
  data {
    int <lower = 0> N; // Sample size
    vector[N] logR_obs;
    real mu_obs;
  }
  
  parameters {
    real a0; // Intercept
    real <lower = 0, upper = 1> a1; 
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
    vector <lower = 0> [N] logR;
    vector <lower = 0> [N] procerr;
    vector <lower = 0> [N] obserr;
  }
  
  model {
    logR[1] ~ normal(logR_obs[1], 0.1);
    logR_obs[1] = logR[1] * obserr[1];
    procerr[1] ~ lognormal(0, sigma_proc);
    obserr[1] ~ lognormal(0, sigma_obs);
    for(t in 2:N){
      logR[t] = a0 + a1 * logR[t-1] * procerr[t];
      procerr[t] ~ lognormal(0, sigma_proc);
      logR_obs[t] = logR[t] * obserr[t];
      obserr[t] ~ lognormal(0, sigma_obs);
    }
   
   
    a0 ~ normal(0,1);
    a1 ~ uniform(0,1);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
  }
  
  generated quantities {
    
    }
  }
