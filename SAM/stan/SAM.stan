
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real <lower =0> a;
      real b;
      simplex [nweight] w; 
      real <lower=0> sigma_proc;
      real <lower=0> sigma_obs;
      vector [N] R_mod;
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
      //R_mod[1:5] = R[1:5];
      for (i in 6:N){
        R_mod[i] ~ normal(b+a*Pant[i], sigma_proc); // likelihood
      }
      
      for (i in 1:N){
        R[i] ~ normal(R_mod[i], sigma_obs);
      }

    b ~ normal(0,5); //priors
    a ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    
    }
