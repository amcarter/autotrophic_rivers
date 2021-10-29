// Stan simple autoregressive model

data {
 int <lower = 0> N; // Sample size
 vector[N] GPP_obs; // observed values, y
 vector[N] light_PAR; // driver 
 real GPP_obs_1;
}

parameters {

 real < lower = 0, upper = 1> phi; 
 real beta0; // Intercept
 real < lower = 0 > beta1; // constrained to positive values 
 real < lower = 0 > sigma;
 real < lower = 0 > omega;
 vector[N] GPP_mod;
}

model {
 
 GPP_mod[1] ~ normal(GPP_obs_1, 1);
 for(i in 2:N){
  GPP_mod[i] ~ normal(phi * GPP_mod[i-1] + beta0 + beta1 * light_PAR[i], sigma);
 }
 GPP_obs ~ normal(GPP_mod, omega);
 
 //priors
 phi ~ uniform(0,1);
 beta0 ~ normal(0,10);
 beta1 ~ lognormal(0,1);
 sigma ~ exponential(1);
 omega ~ exponential(1);
 
}

generated quantities {
} 
