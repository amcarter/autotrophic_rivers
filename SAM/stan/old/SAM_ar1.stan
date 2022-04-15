

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
    
  }
