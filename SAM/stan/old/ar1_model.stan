
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
  }
function () 
.Internal(getwd())
<bytecode: 0x0000017c4f986148>
<environment: namespace:base>
[1] "C:/Users/Alice Carter/git/autotrophic_rivers"
[1] "C:/Users/Alice Carter/git/autotrophic_rivers"
[1] "C:/Users/Alice Carter/git/autotrophic_rivers"
[1] "C:/Users/Alice Carter/git/autotrophic_rivers"

