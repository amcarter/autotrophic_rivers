//
// Stan program
// Respiration model with:
// - Latent state for detrital carbon
//      - log process error
//      - fixed C0
//      - linear scour term
//
// A Carter
// Created: 4/2/22
// Updated: 4/10/22
//

data {
  int<lower=0> ndays;           // length of time series
  int<lower=0> nweights;        // number of antecedent weights
  vector[ndays] R_obs;          // observed respiration (gC/m2/d)
  vector<lower=0>[ndays] P;     // observed GPP (gC/m2/d)
}

transformed data{
  //constants:
  real AR_f; // fraction of autotrophic respiration
  int<lower=0> istart[nweights] = {1,3,7,15}; // starting day of each past interval
  int<lower=0> iend[nweights] = {2,6,14,30}; // starting day of each past interval
  vector[ndays] AR; // autotrophic respiration
  vector[nweights] w_prior;

  AR_f = 0.44;

  AR = -AR_f * P;
  for(i in 1:nweights){
    w_prior[i] = 1.0/nweights;
  }
}

parameters {
  real intercept;
  real<lower=0,upper=1> beta_p;
  simplex[nweights] w;          //weights on intervals of past GPP
  real<lower=0> sigma_obs;
}

transformed parameters{
  vector<lower=0>[ndays] Pant; //antecedent productivity
  vector<lower=0>[ndays+30] Ppre; //GPP vector used to calc ante productivity

  for(i in 1:30){
    Ppre[i] = P[1];
  }
  Ppre[31:(ndays+30)] = P;

  for(i in 31:(ndays)){
    vector[nweights] Pvec;
    for(j in 1:nweights){
      vector[iend[j]-istart[j]+1] pp;
      for(k in istart[j]:iend[j]){
        pp[k-istart[j]+1] = Ppre[i-k];
      }
      Pvec[j] = w[j]*mean(pp);
    }
    Pant[i] = sum(Pvec);
  }
}

model {
  vector[ndays] R;  // actual respiration (unobserved)

  for(i in 31:ndays){
    R[i] = intercept + AR[i] + beta_p * Pant[i];
  }

  // Observation Model:
  R_obs ~ normal(R, sigma_obs);

  // Priors:
  intercept ~ normal(0,1);
  beta_p ~ uniform(0, 1);
  w ~ dirichlet(w_prior);
  sigma_obs ~ normal(0.08, 0.02);
}
