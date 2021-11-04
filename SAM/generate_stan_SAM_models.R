# Stochastic Antecedent Model based on Ogle et al 2015
# base model is from Bob Hall
# 10/25/2021
setwd("C:/Users/Alice Carter/git/autotrophic_rivers")

# Basic SAM code ####
# ER as a function of GPP on days t-4:t

sink("src/SAM/stan/SAM.stan")

cat("
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real a0;
      real <lower =0> a1;
      simplex [nweight] w; 
      real <lower=0> sigma;
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
      for (i in 6:N){
        R[i] ~ normal(a0 + a1*Pant[i], sigma); // likelihood
      }
      

    a0 ~ normal(0,5); //priors
    a1 ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma ~ normal(0,1) T[0,];
    }
    
    generated quantities{
      vector [N] R_hat;
      R_hat[1:5] = R[1:5];
      for(i in 6:N){
        R_hat[i] = normal_rng(a0 + a1 * Pant[i], sigma);
      }
      
    }" ,fill=TRUE)

sink()

# SAM model for different previous intervals of GPP ####
#   ER_t = f(GPP_t, GPP_t-1, GPP_week, GPP_month)

sink("src/SAM/stan/SAM_intervals.stan")

cat("
    data {
      int <lower = 0> N;
      vector [N] P;
      vector [N] R;
      int <lower = 0> nweight;
      vector [nweight] alpha;
    }
    
    parameters {
      real a0;
      real <lower =0> a1;
      simplex [nweight] w; 
      real <lower=0> sigma;
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
        R[i] ~ normal(a0 + a1 * Pant[i], sigma); // likelihood
      }
      
    a0 ~ normal(0,5); //priors
    a1 ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    
    }" ,fill=TRUE)

sink()
