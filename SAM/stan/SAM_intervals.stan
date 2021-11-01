
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
      Pant[1:30] = P[1:30];

      for (i in 31:N){
        real Pweek;
        real Pmonth;
        Pweek = sum(P[(i-7):(i-1)])/7;
        Pmonth = sum(P[(i-30):(i-1)])/30;
        
        Pant[i] = w[1]*P[i] + w[2]*P[i-1] + w[3]*Pweek + w[4]*Pmonth;
        
      }
    }
    
    model {
      for (i in 31:N){
        R_mod[i] ~ normal(b+a*Pant[i], sigma_proc); // likelihood
      }
      
      for (i in 1:N){
        R[i] ~ normal(R_mod[i], sigma_obs);
      }

    b ~ normal(0,5); //priors
    a ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma_obs ~ normal(0,1) T[0,];
    sigma_proc ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    
    }
