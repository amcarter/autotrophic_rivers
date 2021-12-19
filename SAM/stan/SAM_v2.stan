
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      vector [N] temp;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real a0;  // intercept
      real <lower=0, upper=1> ARf; // autotrophic resp fraction
      real <lower=0> a1;
      real <lower=0> a2;
      real <lower=0> Ea;
      simplex [nweight] w; 
      real <lower=0> sigma;
    }
    
    transformed parameters{
      vector [N] Pant;
      vector [N] temp_K;
      temp_K = (temp + 273.15)/100;
      Pant[1:nweight] = P[1:nweight];

      for (i in (nweight + 1):N){
        vector  [nweight] Pvec;
        for(j in 1:nweight){ 
          Pvec[j]=w[j]*P[i-j];
        }
        Pant[i]=sum(Pvec);
      }
    }
    
    model {
      for (i in (nweight+1):N){
        R[i] ~ normal(a0 + ARf*P[i] + a1*Pant[i] + a2*exp(-Ea/temp_K[i]), sigma); // likelihood
      }
      
    //priors
    a0 ~ normal(0,5); 
    a1 ~ normal(0,1);
    a2 ~ normal(0,1);
    Ea ~ normal(0,1);
    ARf ~ beta(10.4, 13.2); //mean = 0.44, sd = 0.1
    w ~ dirichlet(alpha);
    sigma ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    
      
    }
