# adapted from bob's Stan_report_2.rmd

library(rstan)
setwd("C:/Users/Alice Carter/git/autotrophic_rivers")

# Stan model

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
      real <lower =0> a;
      real b;
      simplex [nweight] w; 
      real <lower=0> sigma_proc;
      real <lower=0> sigma_obs;
      vector [N] R_mod;
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
      //R_mod[1:5] = R[1:5];
      for (i in 6:N){
        R_mod[i] ~ normal(b+a*Pant[i], sigma_proc); // likelihood
      }
      
      for (i in 1:N){
        R[i] ~ normal(R_mod[i], sigma_obs);
      }

    b ~ normal(0,5); //priors
    a ~ normal(0,1);
    w ~ dirichlet(alpha);
    sigma_proc ~ normal(0,1) T[0,];
    sigma_obs ~ normal(0,1) T[0,];
    }
    
    generated quantities{
    
    }" ,fill=TRUE)

sink()

##Fake Data: Assign weight with 50% each on lags day 1 and 2

# set parameters
a = 0.5
b = 0
w <- c(0.5,0.5,0,0,0)
sigma_proc = 0.2
sigma_obs = 0.2

# generate data
P <- numeric(100)
R <- numeric(100)
P[1] <- 10
for (i in 2:100)
  P[i] <- 1+0.9* P[i-1]+rnorm(1,0,0.85)

Pant <- numeric(100)
Pant[1:5]<-P[1:5]

for (i in 6:100){
  Pvec<-numeric(5)
  for(j in 1:5){
    Pvec[j]<-w[j]*P[i-(j-1)]
  }
  Pant[i]<-sum(Pvec)
  
}

R_mod <- 0.5 * P[1]
for (i in 2:100){
  R_mod[i] <- b + a*Pant[i] + rnorm(1, 0, sigma_proc)
}
for (i in 1:100){
  R[i] <- rnorm(1, R_mod[i], sigma_obs)
}

plot(P,R)

## How did we do?
sim_dat <- list(R = R, P = P, N = length(P), nweight = 5, alpha = rep(1, 5))
fit_fake <- stan(file = 'src/SAM/stan/SAM.stan', 
                 data = sim_dat, 
                 warmup = 500, iter = 1000, 
                 chains = 4, cores = 4)
print(fit_fake, pars=c("a", "b", "w", 'sigma_proc', 'sigma_obs'))
plot(fit_fake, pars = c("a", "b", "w", 'sigma_proc', 'sigma_obs'))

## Next:
# 1. Find example state space model code:
#     Another of Bob's projects, Kiona's code, stream metabolizer
# 2. Figure out why observation error is trading off with intercept
# 3. Make markdown file that documents what i've tried and what results I got
#     What is not working with state space, how to set first five latent states
# 4. clean up code and push to git

