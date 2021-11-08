
  data {
   int <lower = 0> N; // Sample size
   vector[N] R;
  }
  
  parameters {
   real a0; // Intercept
   real <lower = 0, upper = 1> a1; 
   real <lower = 0> sigma_proc;
   real <lower = 0> sigma_obs;
   vector [N] mu;
  }
  
  model {
   mu[1] ~ normal(R[1], 0.01);
   mu[2:N] ~ normal(a0 + a1 * mu[1:(N-1)], sigma_proc);
   R ~ normal(mu, sigma_obs);
   
   
   a0 ~ normal(0,1);
   a1 ~ uniform(0,1);
   sigma_proc ~ normal(0,1) T[0,];
   sigma_obs ~ normal(0,1) T[0,];
  }
  
  generated quantities {
    vector [N] R_hat;
    vector [N] R_tilde;
    R_hat[1] = R[1];
    R_tilde[1] = R[1];
    for(i in 2:N){
        R_hat[i] = normal_rng(a0 + a1 * R_hat[i-1], sigma_proc);
        R_tilde[i] = normal_rng(R_hat[i], sigma_obs);  
    }
  }
