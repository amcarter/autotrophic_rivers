
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      vector [N] temp;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real <lower=0, upper=1> ARf; // autotrophic resp fraction
      real <lower=0> a1;
      real <lower=0> a2;
      //real <lower=0> Ea; excluded for now because I can't get the magnitudes right
      simplex [nweight] w; 
      real <lower=0> sigma;
    }
    
    transformed parameters{
      vector [N] Pant;
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
        R[i] ~ normal(ARf*P[i] + a1*Pant[i] + a2*exp(-1/temp[i]), sigma); // likelihood
      }
      
    //priors
    a1 ~ beta(1,1);
    a2 ~ normal(0,5);
    //Ea ~ normal(0,1);
    ARf ~ beta(10.4, 13.2); //mean = 0.44, sd = 0.1
    w ~ dirichlet(alpha);
    sigma ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    }
    
    
