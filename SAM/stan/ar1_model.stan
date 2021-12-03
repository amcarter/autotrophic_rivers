
  data {
    int <lower = 0> N; // Sample size
    vector[N] R_obs;   // time series of observations
  }
  
  parameters {
    real a0; // Intercept
    real <lower = 0, upper = 1> phi; //ar1 coeff
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
  }
