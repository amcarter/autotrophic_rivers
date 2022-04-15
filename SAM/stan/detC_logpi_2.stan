//
// Stan program
// Respiration model with:
// - Latent state for detrital carbon
//      - log process error
//      - fixed C0
//
// A Carter
// Created: 4/2/22
// Updated: 4/10/22
//

data {
  int<lower=0> ndays;           // length of time series
  vector[ndays] R_obs;          // observed respiration (gC/m2/d)
  vector<lower=0>[ndays] P;     // observed GPP (gC/m2/d)
  vector[ndays] temp;           // water temperature (C)
  vector<lower=0>[ndays] tau;   // Shear stress (mg/m2) = f(depth, width, slope)
  vector<lower=0>[ndays] litter;// terrestrial litter input (gC/m2/d)
  real<lower=0> C0;             // initial carbon storage
}

transformed data {
 real tau_max;
 real tau_prior;
 real K_20; // decay rate at 20C (1/day)
 real E_a;  // activation energy for carbon decay (eV)
 real k_b;  // boltzmann's constant in eV/K
 vector[ndays] K;  // decay rate of organic carbon

 K_20 = 0.01;
 E_a = 0.63;
 k_b = 8.6173 * 10^(-5);
 K = K_20 * exp(-E_a/k_b * (1.0/(temp + 273) - 0.003412969));

 tau_max = max(tau);
 tau_prior = sort_asc(tau)[(ndays*3)/4];

}

parameters {
  real<lower=0,upper=1> beta_s;
  real<lower=0> sigma_proc;
  real<lower=0> sigma_obs;
  vector<lower=0>[ndays] C; // organic carbon
  real<lower=0> tau0;       // lowest shear stress that causes disturbance
}

transformed parameters{
 vector<lower=0,upper=1>[ndays] ss; //disturbance from shear stress

 for(i in 1:ndays){
     if(tau[1] < tau0)
        ss[i] = 0;
     else
        ss[i] = ((tau[i] - tau0)/(tau_max - tau0))^2;
 }
}


model {
  //constants:
  real AR_f; // fraction of autotrophic respiration
  vector[ndays] R;  // actual respiration (unobserved)
  vector[ndays] Chat; //Carbon before process error

  AR_f = 0.44;

  // Process model

  //initialize
  Chat[1] = C0;
  C[1] ~ lognormal(log(Chat[1]), sigma_proc);
  R[1] = R_obs[1];


  // iterate through timesteps:
  for(i in 2:ndays){
    Chat[i] =  C[i-1] * (1 - beta_s * ss[i]) + litter[i] + R[i-1] + P[i-1];
    C[i] ~ lognormal(log(Chat[i]), sigma_proc);
    R[i] = -AR_f * P[i] - K[i] * C[i];
  }

  // Observation Model:
  R_obs ~ normal(R, sigma_obs);

  // Priors:
  beta_s ~ uniform(0, 1);
  tau0 ~ normal(tau_prior, 1);
  sigma_obs ~ normal(0.08, 0.02);
  sigma_proc ~ normal(0,0.05);
  C ~ lognormal(4.5, 0.75);
}
