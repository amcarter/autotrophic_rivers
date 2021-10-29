
    
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
    real <lower=0> sigma;
  }
    
  transformed parameters{
    vector  [N] Pant;
    Pant[1:5]= P[1:5];

    for (i in 6:N){
      vector  [nweight] Pvec;
      for(j in 1:nweight){
        Pvec[j]=w[j]*P[i-(j-1)];
      }
      Pant[i]=sum(Pvec);
    }

  }
    
  model {
    for ( i in 6:N){
      R[i] ~ normal(b+a*Pant[i], sigma); // likelihood
    }

    b~normal(0,5); //priors
    a~normal(0,1);
    w~dirichlet(alpha);
  }

  generated quantities{

  }
