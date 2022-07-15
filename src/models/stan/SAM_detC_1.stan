//
// Stan program
// Respiration model with:
// - Stochastic antecedent GPP
// - Latent state for detrital carbon
//
// A Carter
// Created: 4/2/22
// Updated:
//

data {
  int<lower=0> ndays;           // length of time series
  int<lower=0> nweights;        // number of antecedent intervals
  vector[ndays] R_obs;          // observed respiration (gC/m2/d)
  vector<lower=0>[ndays] P;     // observed GPP (gC/m2/d)
  vector[ndays] temp;           // water temperature (C)
  vector<lower=0>[ndays] tau;   // Shear stress (mg/m2) = f(depth, width, slope)
  vector<lower=0>[ndays] litter;// terrestrial litter input (gC/m2/d)
  real C0;
}

transformed data {
  vector<lower=0>[ndays] temp_K;// temperature (K)
  vector<upper=0>[ndays] AR;    // Autotrophic respiration
  vector[nweights] w_prior;     // simplex prior
  real<lower=0> tau_max;        // max shear stress
  real<lower=0> tau_prior;      // prior value on tau0 - set using a quantile on tau

  //constants:
  real AR_f; // fraction of autotrophic respiration
  real E_a;  // activation energy for carbon decay (eV)
  real k_b;  // boltzmann's constant in eV/K
  real K_20;          // decay rate at 20C (1/day)
  vector<lower=0>[ndays] K;       // decay rate per day

  AR_f = 0.44;
  E_a = 0.63;
  k_b = 8.6173 * 10^(-5);
  K_20 = 0.01;

  temp_K = temp + 273;
  AR = -AR_f * P;
  for(i in 1:nweights){
    w_prior[i] = 1.0/nweights;
  }
  tau_max = max(tau);
  tau_prior = sort_asc(tau)[ndays/2];
  K = K_20 * exp((E_a/k_b)*(1.0./temp_K - 1.0/293.0));
}

parameters {
  // real<lower=0> C0;            // initial detrital carbon in the system (g/m2)
  vector<lower=0>[ndays] C;    // detrital carbon at each timestep (g/m2)
  real<lower=0,upper=1> beta_s;// percent carbon lost at max shear stress
  real<lower=0> tau0;          // minimum shear stress for bed disturbance (mg/m2)
  real<lower=0> beta_p;        // coefficient on antecedent primary productivity
  simplex[nweights] w;         // weights on intervals of antecedent GPP
  real<lower=0> sigma_proc;    // process error
  real<lower=0> sigma_obs;     // observation error
}

transformed parameters {
  vector<lower=0>[ndays] Pant;    // antecedent productivity
  vector<lower=0>[ndays] ss_ratio;// squared ratio of shear stress in loss term

  Pant[1:nweights] = P[1:nweights];

  for(i in (nweights + 1):ndays){
    vector[nweights] Pvec;
    for(j in 1:nweights){
      Pvec[j] = w[j]*P[i-j];
    }
    Pant[i] = sum(Pvec);
  }

  for(i in 1:ndays){
    if(tau[i] < tau0){
      ss_ratio[i] = 0;
    } else {
      ss_ratio[i] = ((tau[i] - tau0)/(tau_max - tau0))^2;
    }
  }


}

model {
  vector[ndays] R;    // actual respiration (gC/m2/d)
  vector[ndays] HR_d; // respiration of detrital carbon
  vector[ndays] HR_p; // respiration of algal material
  vector[ndays] C_hat;// Carbon before process error

  // process model
  C_hat[1] = C0;
  C[1] ~ lognormal(log(C_hat[1]), sigma_proc);
  HR_d[1] = -K[1] * C[1];
  HR_p = -beta_p * Pant;
  print("C0 = ", C0, "C_hat = ", C_hat[1], " beta_s = ", beta_s, " ss = ", ss_ratio[1], " K20 = ", K_20, K[1]);

  for(i in 2:ndays){
    C_hat[i] = (C[i-1] + HR_d[i-1]) * (1 - beta_s * ss_ratio[i]) + litter[i];
    C[i] ~ lognormal(log(C_hat[i]), sigma_proc);
    HR_d[i] = -K[i] * C[i];
  }

  // observation model
  R = AR + HR_d + HR_p;
  R_obs ~ normal(R, sigma_obs);

  // priors
  // C0 ~ normal(100, 5);
  C[2:ndays] ~ lognormal(log(C[1:(ndays-1)]), 0.75);
  K_20 ~ normal(1, 0.1)T[0,];
  beta_s ~ beta(2.7, 1);
  tau0 ~ normal(tau_prior, 1);
  beta_p ~ normal(.5, .5)T[0,];
  w ~ dirichlet(w_prior);
  sigma_proc ~ normal(0, 0.02)T[0,];
  sigma_obs ~ normal(0.08, 0.02)T[0,];
}
