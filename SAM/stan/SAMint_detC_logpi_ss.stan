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
    vector<lower=0>[ndays] tau;   // normalized shear stress = f(depth, width, slope)
    real<lower=0> C0;             // initial carbon storage
    matrix<lower=0>[ndays,nweights] ANT; //antecedent conditions matrix
}

transformed data{
    // autotrophic respiration
    real AR_f = 0.44; // fraction of autotrophic respiration
    vector[ndays] AR = -AR_f * P; // autotrophic respiration
    // shear stress parameters
    real tau_max = max(tau);

    vector[nweights] w_prior;
    for(i in 1:nweights){
        w_prior[i] = 1;
    }

}

parameters {
    real<lower=0> beta_p;       // coefficient on HR_algal
    real<lower=0> K_20_scaled;  // 100x decay rate at 20C (1/day)
    vector<lower=0>[ndays] C;   // organic carbon
    real<lower=0,upper=1> tau0; // lowest shear stress that causes disturbance
    simplex[nweights] w;        // weights on intervals of past GPP
    real<lower=0> sigma_proc_scaled;
    real<lower=0> sigma_obs;
}

transformed parameters{
    //rescaled parameters:
    real K_20 = K_20_scaled/100;
    real sigma_proc = sigma_proc_scaled/10;
    //constants:
    real E_a = 0.63;  // activation energy for carbon decay (eV)
    real k_b = 8.6173 * 10^(-5);  // boltzmann's constant in eV/K
    vector[ndays] K;  // decay rate of organic carbon
    vector[ndays] d2;
    vector<lower=0>[ndays] Pant; //yesterday's productivity

    Pant = ANT * w;

    K = K_20 * exp(-E_a/k_b * (1.0 ./ (temp+273) - 0.003412969));
    for(i in 1:ndays){
        if(tau[i] < tau0){
            d2[i] = 0;
        } else {
            d2[i] = 1/(tau_max - tau0);
        }
    }
}

model {
    vector[ndays] R;  // actual respiration (unobserved)
    vector[ndays] HR;  // detrital carbon respiration (unobserved)
    vector[ndays] Chat; //Carbon before process error
    vector[ndays] disturb;//disturbance from shear stress

    for(i in 1:ndays){//assumes 80% loss from largest storm
        disturb[i] = 1 - 0.8 * (((tau[i] - tau0) * d2[i]) ^ 2);
    }

    // Process model

    //initialize
    for(i in 1:antdays){
        Chat[i] = C0;
        C[i] ~ lognormal(log(Chat[i]), sigma_proc);
        HR[i] = -K[i] * C[i];
    }

    // iterate through timesteps:
    for(i in (antdays+1):ndays){
      Chat[i] = C[i-1]*disturb[i] + litter[i] + HR[i-1];
      C[i] ~ lognormal(log(Chat[i]), sigma_proc);
      HR[i] = -K[i] * C[i];
    }

    // Observation Model:
    R = HR + AR - beta_p * Pant;
    R_obs ~ normal(R, sigma_obs);

    // Priors:
    beta_p ~ normal(0, 1);
    tau0 ~ normal(0.5, 0.3);
    w ~ dirichlet(w_prior);
    sigma_obs ~ normal(0.08, 0.02);
    sigma_proc_scaled ~ normal(0,0.1);
    K_20_scaled ~ normal(1, 0.1);
    C ~ lognormal(4.5, 0.75);
}
