
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
  }
