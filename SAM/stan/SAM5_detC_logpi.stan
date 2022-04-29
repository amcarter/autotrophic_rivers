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
    int<lower=0> antdays;         // total number of antecedent days
    vector[ndays] R_obs;          // observed respiration (gC/m2/d)
    vector<lower=0>[ndays] P;     // observed GPP (gC/m2/d)
    vector[ndays] temp;           // water temperature (C)
    vector<lower=0>[ndays] litter;// terrestrial litter input (gC/m2/d)
    vector<lower=0>[ndays] Q;     // discharge
    real<lower=0> C0;             // initial carbon storage
    matrix<lower=0>[ndays,nweights] ANT; //antecedent conditions matrix
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
    real<lower=0,upper=1> beta_s;
    real<lower=0,upper=1> beta_p;
    simplex[nweights] w;          //weights on intervals of past GPP
    real<lower=0> sigma_proc;
    real<lower=0> sigma_obs;
    vector<lower=0>[ndays] C; // organic carbon
    real<lower=0> K_20; // decay rate at 20C (1/day)
}

transformed parameters{
    vector[ndays] K;  // decay rate of organic carbon
    vector<lower=0>[ndays] Pant; //yesterday's productivity
    real E_a;  // activation energy for carbon decay (eV)
    real k_b;  // boltzmann's constant in eV/K
    E_a = 0.63;
    k_b = 8.6173 * 10^(-5);

    for(i in 1:ndays){
        K[i] = K_20/100 * exp(-E_a/k_b * (1.0/(temp[i]+273) - 0.003412969));
    }

    Pant = ANT * w;
}

model {
    vector[ndays] R;  // actual respiration (unobserved)
    vector[ndays] HR;  // detrital carbon respiration (unobserved)
    vector[ndays] Chat; //Carbon before process error

    // Process model

    //initialize
    for(i in 1:antdays){
        Chat[i] = C0;
        C[i] ~ lognormal(log(Chat[i]), sigma_proc/10);
        HR[i] = -K[i] * C[i];
    }
    // iterate through timesteps:
    for(i in (antdays+1):ndays){
        Chat[i] = (C[i-1] + litter[i] + HR[i-1])*(1 - beta_s*Q[i]);
        C[i] ~ lognormal(log(Chat[i]), sigma_proc/10);
        HR[i] = - K[i] * C[i];
    }

    // Observation Model:
    R = HR + AR - beta_p * Pant;
    R_obs ~ normal(R, sigma_obs);

    // Priors:
    beta_s ~ beta(3, 1.25);
    beta_p ~ uniform(0, 1);
    w ~ dirichlet(w_prior);
    K_20 ~ normal(1, 0.1);
    sigma_obs ~ normal(0.08, 0.02);
    sigma_proc ~ normal(0,0.1);
    C ~ lognormal(4.5, 0.75);
}
