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
    int<lower=0> antdays;         // total number of antecedent days being considered
    vector[ndays] R_obs;          // observed respiration (gC/m2/d)
    vector<lower=0>[ndays] P;     // observed GPP (gC/m2/d)
    matrix<lower=0>[ndays,nweights] ANT;// antecedent conditions
}

transformed data{
    //constants:
    real AR_f; // fraction of autotrophic respiration
    vector[ndays] AR; // autotrophic respiration
    vector[nweights] w_prior;

    AR_f = 0.44;

    AR = -AR_f * P;
    for(i in 1:nweights){
        w_prior[i] = 1;
    }

}

parameters {
    real beta_0;
    real<lower=0> beta_p;
    simplex[nweights] w;          //weights on intervals of past GPP
    real<lower=0> sigma_obs;
}

transformed parameters{
    vector<lower=0>[ndays] Pant; //yesterday's productivity
    Pant = ANT * w;
}

model {
  vector[ndays] R;  // actual respiration (unobserved)

  for(i in (antdays+1):ndays){
    R[i] = beta_0 + AR[i] + beta_p * Pant[i];
    R_obs[i] ~ normal(R[i], sigma_obs);
  }


  // Priors:
  beta_0 ~ normal(0, 1);
  beta_p ~ normal(0, 5);
  w ~ dirichlet(w_prior);
  sigma_obs ~ normal(0.08, 0.02);

 }
