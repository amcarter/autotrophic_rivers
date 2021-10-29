# Stochastic Antecedent Model based on Ogle et al 2015
# base model is from Bob Hall
# 10/25/2021
setwd('C:/Users/Alice Carter/git/autotrophic_rivers')
sink("src/SAM/stan/SAM_pi_5.stan")

cat("
    
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

  }", fill=TRUE)

sink()


# Model including observation error:

sink('src/SAM/stan/SAM_oipi_5.stan')

cat("

  data {
    int <lower = 0> N;
    vector [N] P;
    vector [N] R;
    int <lower = 0> nweight;
    vector [nweight] alpha;
  }
  
  parameters {
    real <lower = 0> a;
    real b;
    simplex [nweight] w;
    real <lower = 0> sigma_proc;
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
      R[i] ~ normal(b + a * Pant[i], sigma_proc); //likelihood
    }
    
    b ~ normal(0,5); // priors
    a ~ normal(0,1);
    w ~ dirichlet(alpha);
  }
  
  generated quantities{
  
  }", fill = TRUE)

sink()

