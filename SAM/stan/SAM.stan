
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
      
    }
[[1]]
Stan model 'ar1_model' does not contain samples.

[[2]]
Stan model 'ar1_model' does not contain samples.

[[3]]
Stan model 'ar1_model' does not contain samples.

[[4]]
Stan model 'ar1_model' does not contain samples.

Stan model 'ar1_model' does not contain samples.
Stan model 'ar1_model' does not contain samples.
Stan model 'ar1_model' does not contain samples.

SAMPLING FOR MODEL 'ar1_model' NOW (CHAIN 1).
Chain 1: 
Chain 1: Gradient evaluation took 0 seconds
Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
Chain 1: Adjust your expectations accordingly!
Chain 1: 
Chain 1: 
Chain 1: Iteration:    1 / 5000 [  0%]  (Warmup)
Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc512d1741_ar1_model' at line 35)

Chain 1: Iteration:  500 / 5000 [ 10%]  (Warmup)
Chain 1: Iteration:  501 / 5000 [ 10%]  (Sampling)
Chain 1: Iteration: 1000 / 5000 [ 20%]  (Sampling)
Chain 1: Iteration: 1500 / 5000 [ 30%]  (Sampling)
Chain 1: Iteration: 2000 / 5000 [ 40%]  (Sampling)
Chain 1: Iteration: 2500 / 5000 [ 50%]  (Sampling)
Chain 1: Iteration: 3000 / 5000 [ 60%]  (Sampling)
Chain 1: Iteration: 3500 / 5000 [ 70%]  (Sampling)
Chain 1: Iteration: 4000 / 5000 [ 80%]  (Sampling)
Chain 1: Iteration: 4500 / 5000 [ 90%]  (Sampling)
Chain 1: Iteration: 5000 / 5000 [100%]  (Sampling)
Chain 1: 
Chain 1:  Elapsed Time: 1.286 seconds (Warm-up)
Chain 1:                15.309 seconds (Sampling)
Chain 1:                16.595 seconds (Total)
Chain 1: 

SAMPLING FOR MODEL 'ar1_model' NOW (CHAIN 1).
Chain 1: 
Chain 1: Gradient evaluation took 0 seconds
Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
Chain 1: Adjust your expectations accordingly!
Chain 1: 
Chain 1: 
Chain 1: Iteration:    1 / 5000 [  0%]  (Warmup)
Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Iteration:  500 / 5000 [ 10%]  (Warmup)
Chain 1: Iteration:  501 / 5000 [ 10%]  (Sampling)
Chain 1: Iteration: 1000 / 5000 [ 20%]  (Sampling)
Chain 1: Iteration: 1500 / 5000 [ 30%]  (Sampling)
Chain 1: Iteration: 2000 / 5000 [ 40%]  (Sampling)
Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Iteration: 2500 / 5000 [ 50%]  (Sampling)
Chain 1: Iteration: 3000 / 5000 [ 60%]  (Sampling)
Chain 1: Iteration: 3500 / 5000 [ 70%]  (Sampling)
Chain 1: Iteration: 4000 / 5000 [ 80%]  (Sampling)
Chain 1: Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)

Chain 1: Iteration: 4500 / 5000 [ 90%]  (Sampling)
Chain 1: Iteration: 5000 / 5000 [100%]  (Sampling)
Chain 1: 
Chain 1:  Elapsed Time: 1.152 seconds (Warm-up)
Chain 1:                9.202 seconds (Sampling)
Chain 1:                10.354 seconds (Total)
Chain 1: 
 [1] "Chain 1: Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:"                        
 [2] "Chain 1: Exception: normal_lpdf: Location parameter[1] is inf, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 18)"              
 [3] "Chain 1: If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,"
 [4] "Chain 1: but if this warning occurs often then your model may be either severely ill-conditioned or misspecified."                              
 [5] "Chain 1: "                                                                                                                                      
 [6] "Chain 1: Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:"                        
 [7] "Chain 1: Exception: normal_lpdf: Scale parameter is 0, but must be > 0!  (in 'model37dc144e2924_ar1_model' at line 19)"                         
 [8] "Chain 1: If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,"
 [9] "Chain 1: but if this warning occurs often then your model may be either severely ill-conditioned or misspecified."                              
[10] "Chain 1: "                                                                                                                                      
[11] "Error in sampler$call_sampler(args_list[[i]]) : "                                                                                               
[12] "  Exception: normal_rng: Location parameter is nan, but must be finite!  (in 'model37dc144e2924_ar1_model' at line 35)"                         
[13] "In addition: Warning message:"                                                                                                                  
[14] "In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :"                                                                      
[15] "  'C:/Rtools/mingw_/bin/g++' not found"                                                                                                         
[1] "error occurred during calling the sampler; sampling not done"
Stan model 'ar1_model' does not contain samples.
Stan model 'ar1_model' does not contain samples.
Stan model 'ar1_model' does not contain samples.
Stan model 'ar1_model' does not contain samples.
Stan model 'ar1_model' does not contain samples.
[[1]]
Stan model 'ar1_model' does not contain samples.

[[2]]
Stan model 'ar1_model' does not contain samples.

[[3]]
Stan model 'ar1_model' does not contain samples.

[[4]]
Stan model 'ar1_model' does not contain samples.

