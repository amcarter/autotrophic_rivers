

  data {
    int <lower = 0> N;
    vector [N] P;
    vector [N] R;
    int <lower = 0> nweight;
    vector [nweight] alpha;
  }
  
  parameters {
    vector [N] mu; // latent, unobserved R values
    real <lower = 0> a;
    real b;
    simplex [nweight] w;
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
  }
  
  transformed parameters {
    vector [N] Pant;
    Pant[1:nweight] = P[1:nweight];
    
    for(i in (nweight + 1):N){
      vector [nweight] Pvec;
      for(j in 1:nweight){
        Pvec[j] = w[j] * P[i-(j-1)];
      }
      Pant[i] = sum(Pvec);
    }
  }
  
  model {
    for(i in (nweight + 1):N){
      mu[i] ~ normal(b + a * Pant[i], sigma_proc); //likelihood
      R[i] ~ normal(mu[i], sigma_obs);
    }
    
    // priors
    b ~ normal(0,5); 
    a ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(0,1) T[0,];
  }
  
  generated quantities{
  
  }
