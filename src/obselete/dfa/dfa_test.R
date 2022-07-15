# Dynamic factor analysis on GPP values in powell center
# This is a trial run using the greta package in R

# A Carter
# 30 Sept 2021

# setup #### 
# This script is vollowing the greta vignette at 
# https://mdscheuerell.github.io/gretaDFA/
# install.packages('greta')
# install.packages('tensorflow')
# install.packages('reticulate')

library(greta)
library(tensorflow)
library(MASS)
library(viridis)
library(corrplot)
installed.packages()[c("greta", "tensorflow"), "Version"] # v 0.3.1, v 2.6.0
tensorflow::tf_version() # 2

#these are needed only for plotting, and they can be heinous to install
# install.packages("igraph")
# install.packages("DiagrammeR")
library(igraph)
library(DiagrammeR)


# functions ####
zd <- function(M, sigma) {
  zd <- zeros(M)
  for(i in 1:M) {
    zd[i] <- ld(i, M = M, sigma = sigma, dim = 1)
  }
  return(zd)
}

# Make the Z matrix
# explicitly define the distribution needed to make the diagonal of the Z matrix
library (R6)
distrib <- .internals$nodes$constructors$distrib
distribution_node <- .internals$nodes$node_classes$distribution_node
check_dims <- .internals$utils$checks$check_dims
as.greta_array <- .internals$greta_arrays$as.greta_array
is_scalar <- .internals$utils$misc$is_scalar
fl <- .internals$utils$misc$fl

# function to call 
ld <- function (i, M, sigma, dim = 1) {
  distrib("ld", i, M, sigma, dim)
}

# Leung & Drton distn for prior of diag(Z)
ld_distribution <- R6Class (
  "ld_distribution",
  inherit = distribution_node,
  public = list(
    
    initialize = function (i, M, sigma, dim) {
      
      # check if (i, M, dim) are in counting set
      # (i)
      if (length(i) > 1 ||
          i <= 0 ||
          !is.finite(i) ||
          i != floor(i)) {
        
        stop ("i must be a scalar positive integer, but was: ",
              capture.output(dput(i)),
              call. = FALSE)
        
      }
      
      # (M)
      if (length(M) > 1 ||
          M <= 0 ||
          !is.finite(M) ||
          M != floor(M)) {
        
        stop ("M must be a scalar positive integer, but was: ",
              capture.output(dput(M)),
              call. = FALSE)
        
      }
      
      # (dim)
      if (length(dim) > 1 ||
          dim <= 0 ||
          !is.finite(dim) ||
          dim != floor(dim)) {
        
        stop ("dim must be a scalar positive integer, but was: ",
              capture.output(dput(dim)),
              call. = FALSE)
        
      }
      
      # check if i > M
      if (i > M) {
        
        stop ("i can be no larger than M",
              call. = FALSE)
        
      }
      
      # check if sigma is scalar
      # (sigma)
      if (length(sigma) > 1) {
        
        stop ("sigma must be a scalar positive integer, but was: ",
              capture.output(dput(sigma)),
              call. = FALSE)
        
      }
      
      i <- as.greta_array(i)
      M <- as.greta_array(M)
      sigma <- as.greta_array(sigma)
      
      self$bounds <- c(0, Inf)
      super$initialize("ld", dim, truncation = c(0, Inf))
      self$add_parameter(i, "i")
      self$add_parameter(M, "M")
      self$add_parameter(sigma, "sigma")
      
    },
    
    tf_distrib = function (parameters, dag) {
      
      i <- parameters$i
      M <- parameters$M
      sigma <- parameters$sigma
      
      # log pdf(x | i, M, sigma)
      log_prob = function (x) {
        (M - i) * tf$log(x) - x ^ fl(2) / (fl(2) * sigma)
      }
      
      list(log_prob = log_prob, cdf = NULL, log_cdf = NULL)
      
    },
    
    # no CDF for discrete distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL
    
  )
)

# simulate data ####
NN <- 30 # number of time series
TT <- 30 # number of time points
MM <- 3  # number of latent trends we will try to fit

# latent factors
set.seed(123)
## MM x TT matrix of innovations
ww <- matrix(rnorm(MM*TT, 0, 1), MM, TT) 
ww[,1] <- rnorm(MM, 0, sqrt(5))
## MM x TT matrix of scaled latent trends
xx <- t(scale(apply(ww,1,cumsum)))

plot(seq(1,30), xx[1,], type = 'b', ylim = c(-2,2), xlab = 'time', ylab = 'x_t',
     main = 'cumulative proc error')
lines(seq(1,30), xx[2,], type = 'b', col = 'goldenrod')
lines(seq(1,30), xx[3,], type = 'b', col = 'purple')

ZZ <- matrix(runif(NN*MM, -1, 1), NN, MM)
diag(ZZ) <- rev(sort(abs(diag(ZZ))))
ZZ[upper.tri(ZZ)] <- 0
ZZ <- round(ZZ, 2)

## obs var
obs_var <- 0.2^2
## obs errors
ee <- t(mvrnorm(TT, matrix(0,NN,1), diag(obs_var,NN,NN)))
## NN x TT matrix of observed data
yy <- ZZ %*% xx + ee
plot(1:30, yy[1,], type = 'l', ylim = c(-3,3))
for(i in 2:30){ lines(1:30, yy[i,])}


corrplot(cor(yy, yy), method="ellipse", type="lower",
         tl.col = "black", tl.srt = 0, tl.cex = 0.8, tl.offset = 0.7,
         cl.cex = 0.8, cl.offset = 0.9, cl.ratio = 0.2)
# Priors ####
## empty loadings matrix
ZZ_est <- zeros(NN,MM)
## define sigma
sigma_est = normal(0, 2, truncation = c(0, Inf))
## diagonal
idx_d <- row(ZZ_est) == col(ZZ_est)
ZZ_est_raw_d = zd(MM, sigma_est)
ZZ_est[idx_d] <- ZZ_est_raw_d
## sub-diagonal
idx_s <- lower.tri(ZZ_est)
ZZ_est_raw_s = normal(0, sigma_est, dim = sum(idx_s), truncation = c(-1,1))
ZZ_est[idx_s] <- ZZ_est_raw_s


head(calculate(ZZ_est, list(ZZ_est_raw_d = rep(1, sum(idx_d)),
                            ZZ_est_raw_s = rep(2, sum(idx_s)))))

# covariance matrix
## diagonal of R
RR_est = inverse_gamma(alpha = 1, beta = 3, 1, truncation=c(0, 1))

## inital factor
X0 = normal(mean = 0, sd = sqrt(10), dim = c(MM, 1))

## factors for t = 2-TT
XX= normal(mean = 0, sd = 1, dim = c(MM, TT))
## cumsum over proc errs for the random walk
xx_est <- t(apply(cbind(X0,XX), 1, "cumsum"))[,-1]

## scale data
yy_z <- t(scale(t(yy), scale = FALSE))
## vectorize data
yy_vec <- as_data(matrix(yy_z, 1, NN*TT))
## vectorize mean
Zx_vec <- t(c(ZZ_est %*% xx_est))
## define likelihood
distribution(yy_vec) = normal(Zx_vec, RR_est)

mod_defn <- model(xx_est, ZZ_est, RR_est, sigma_est)
mod_fit <- mcmc(mod_defn, verbose = FALSE,
                sampler = hmc(Lmin = 1, Lmax = 30, epsilon = 0.001, diag_sd = 1),
                warmup = 2000, n_samples = 5000, thin = 10, chains = 1,
                initial_values = initials(RR_est = runif(1, 0.1, 1),
                                          sigma_est = runif(1, 0.1, 1)
                                          # XX = matrix(rnorm(MM*TT), MM, TT),
                                          # ZZ_est_raw_s = matrix(rnorm(sum(idx_s),
                                          #                             0,
                                          #                             0.1),
                                          #                       ncol = 1)
                )
)
# mod_fit <- extra_samples(mod_fit, n_samples = 5000, thin = 50, verbose = FALSE)# sampling

# Assessing convergence ####
# plot the chains:
par(mfcol=c(NN,MM), mai=rep(0,4), omi=rep(0,4))
cnt <- 1
for(j in 1:MM) {
  for(i in 1:NN) {
    coda::traceplot(mod_fit[[1]][,paste0("ZZ_est[",i,",",j,"]")],
                    ylab="", xlab="", yaxt="n", xaxt="n")
    cnt <- cnt + 1
  }
}

# use the Heidelberger and Welch's convergence diagnostic
## number of estimated params/states
n_par <- dim(mod_fit[[1]])[2]
## H&W diag: pass (1) or fail (NA)
halfwidth_test <- coda::heidel.diag(mod_fit[[1]])[,4]
## proportion passing
round(sum(halfwidth_test, na.rm = TRUE) / n_par, 2)

# Examine the fits ####
mod_smry <- summary(mod_fit)
par_means <- mod_smry$statistics[,"Mean"]

# loadings matrix
## fitted Z
ZZ_fit <- par_means[grepl("ZZ",rownames(mod_smry$statistics))]
ZZ_fit <- matrix(ZZ_fit, NN, MM, byrow=FALSE)
## fitted Z
head(round(ZZ_fit, 2))

# correlation between actual Z and estimated Z
round(cor(as.vector(ZZ), as.vector(ZZ_fit)), 2)
par(mai=rep(0.1, 4), omi=rep(0.1, 4))
corrplot(cor(ZZ, ZZ_fit), method="ellipse", type="lower",
         tl.col = "black", tl.srt = 0, tl.cex = 0.8, tl.offset = 0.7,
         cl.cex = 0.8, cl.offset = 0.9, cl.ratio = 0.2)

# Factors
## fitted factors
xx_fit <- par_means[grepl("xx",rownames(mod_smry$statistics))]
xx_fit <- matrix(xx_fit, MM, TT, byrow=FALSE)

par(mai=rep(0.1, 4), omi=rep(0.1, 4))
corrplot(cor(t(xx), t(xx_fit)), method="ellipse", type="lower",
         tl.col = "black", tl.srt = 0, tl.cex = 0.8, tl.offset = 0.7,
         cl.cex = 0.8, cl.offset = 0.9, cl.ratio = 0.2)

round(par_means[grepl("RR",rownames(mod_smry$statistics))], 3)


## fitted values
yy_fit <- ZZ_fit %*% xx_fit
## corrrelation
cor_yy <- matrix(diag(cor(t(yy_z), t(yy_fit))), NN/5, 5)
## plots
#correlation between fitted and observed datasets
par(mai=rep(0.1, 4), omi=rep(0.1, 4))
corrplot(cor_yy, method="ellipse",
         tl.pos = "n",
         cl.cex = 0.8, cl.offset = 0.9,
         cl.ratio = 0.2) #cl.lim = c(0, 1))

sigma_fit <- par_means["sigma_est"]
round(sigma_fit, 2)


## rotation matrix
HH_inv <- varimax(ZZ_fit)$rotmat
## rotated Z
ZZ_rot <- ZZ_fit %*% HH_inv
round(ZZ_rot, 2)
