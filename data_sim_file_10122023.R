# install.packages("rstan")
# install.packages("R2jags")
# install.packages("rjags")
# install.packages("blavaan")
# install.packages("mvtnorm")

# library(rjags) # issue with jags instalation
# library(R2jags)
library(blavaan) # standard package
library(rstan)
library(mvtnorm) # mv normal data
source("functions_for_blavaan.R")
packageVersion("blavaan")
packageVersion("rjags")
packageVersion("rstan")

# BNFI, BTLI, BCFI, BMc, BΓadj, BΓˆ, BRMSEA //// PPPχ2
###########################################
# generate real data
###########################################
set.seed(1024)
N <- 250
phi <- diag(2)
eta <- rmvnorm(N,sigma=phi)
epsilon <- rmvnorm(N,sigma=diag(6))
lambda <- matrix(c(1,1,1,0,0,0,
                   0,0,0,1,1,1),6,2,byrow=F) # misspecify lambda -> 1 where a zero is (one of them)
x <- eta%*%t(lambda)+epsilon #observed data
###########################################
# generate fake data
###########################################
lambda_fake <- matrix(c(1,1,1,0,0,0,
                   1,0,0,1,1,1),6,2,byrow=F)
x_fake <- eta%*%t(lambda_fake)+epsilon#observed data
###########################################
# data analysis with real data
###########################################
# 1st model 
model <- '
eta1 =~ V1+V2+V3
eta2 =~ V4+V5+V6
'
# run the bayesian model
res1 <- bsem(model, data=as.data.frame(x),#target="jags",     ##################COnvergence rate , rhat stat 1, ess above 400
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1", std.lv = TRUE) # this is wrong
# extract model fit
res1 # summary(res1)
summary(res1)
'
  Statistic                                 MargLogLik         PPP
  Value                                             NA       0.540
'
blavFitIndices(res1)
'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.017        0.998        0.994        0.997 
'
#+++++++++++++++++++++++++++++++++
# https://stats.stackexchange.com/questions/532019/how-to-get-rmsea-cfi-of-blavaan-object-in-r
fit1_REAL_INDICIES <- bcfa(model, data = as.data.frame(x),
             n.chains = 4, burnin = 1000, sample = 1000)
null.model <- c(paste0("V", 1:6, " ~~ V", 1:6), paste0("V", 1:6, " ~ 1"))
fit0_REAL_INDICIES <- bcfa(null.model, data = as.data.frame(x),
             n.chains = 4, burnin = 1000, sample = 1000)

blavFitIndices(fit1_REAL_INDICIES, baseline.model = fit0_REAL_INDICIES)
'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc         BCFI         BTLI         BNFI 
       0.000        0.963        1.065        0.942        0.917        1.063        0.954 
Warning message:
'
#+++++++++++++++++++++++++++++++++
# 2nd model: jags
data.jags <- list(N=N,           #sample size
                  y=x,         #data
                  mu.0=c(0,0),   #zero-vector
                  phi.0=phi) #diagonal matrix

params <- c("ly","sigma.eps","sigma.eta")
model_real_jags <- jags.parallel(data=data.jags, 
                                 parameters.to.save=params,
                                 n.iter=2000, n.chains=4,n.thin=1,n.burnin = 1000,
                                 model.file="model1_cfa.txt")
model_real_jags

# 3rd model: stan
stan_data <- list(
  N=N,
  x=x
)

#######################################STAN CODE###########################
stan_model_code <- '
data {
  int<lower=0> N;
  matrix[N, 6] x; 
  // vector[4] ly;
}

parameters {
  matrix[N, 2]  eta;
  vector<lower=0>[6] epsilon; 
  cov_matrix[2] sigma; 
  vector[2] mu_0; 
  vector[4] ly;
}

transformed parameters {
  matrix[N, 6] mu_x;
  
  for (i in 1:N) {
    mu_x[i, 1] = eta[i, 1];
    mu_x[i, 2] = ly[1]*eta[i, 1];
    mu_x[i, 3] = ly[2]*eta[i, 1];

    mu_x[i, 4] = eta[i, 2];
    mu_x[i, 5] = ly[3]*eta[i, 2];
    mu_x[i, 6] = ly[4]*eta[i, 2];
  }
}

model {
  for (i in 1:N) {
    for (j in 1:6) {
      x[i,j] ~ normal(mu_x[i,j], epsilon[j]);
    }
    eta[i, 1:2] ~ multi_normal(mu_0,sigma);
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = 0;
    for (j in 1:6) {
      log_lik[i] += normal_lpdf(x[i,j] | mu_x[i,j], epsilon[j]);
    }
  }
}
'



#######################################STAN CODE###########################

stan_data <- list(
  N=N,
  x=x
  # ly=c(1,1, 1, 1)
)


# find.package('rstan')
model_real_stan <- stan(model_code = stan_model_code, iter = 1000, verbose = FALSE, data=stan_data )
posterior_samples <- extract(model_real_stan)
posterior_samples
posterior_samples$log_lik
posterior_samples$lp

posterior_samples$eta
posterior_samples$sigma
posterior_samples$epsilon
posterior_samples$mu_0
posterior_samples$mu_x
posterior_samples$ly

# BayesRelFit(obs, rep, nvar=6, pD, N=N, null_model=FALSE)
BayesRelFit(obs=1, rep=1, nvar=6, pD=1, N=N, null_model=FALSE)

# 
# ff <- BayesRelFit(obs=chisq.obs, rep=chisq.boot, pD=temp1$pd)
# 
# pd <- loo(blavaan:::case_lls(lavjags, origlavmodel, lavpartable, 
#                              lavsamplestats, lavoptions, lavcache, 
#                              origlavdata))$p_loo
# 
# 
# chisq.obs <- -2*(samplls[i, j, 1] -
#                    samplls[i, j, 2])
# 
# chisq.boot <- 2*diff(blavaan:::get_ll(lavmodel = lavmodel,
#                                       lavsamplestats = lavsamplestats,
#                                       lavdata = lavdata,
#                                       measure = measure))



# fit_indices <- blavFitIndices(model_real_stan)
# print(fit_indices)

# https://discourse.mc-stan.org/t/covariance-matrix-not-symmetric-for-lkj-prior/20932
# https://discourse.mc-stan.org/t/covariance-matrix-not-symmetric/29824/2
# https://mc-stan.org/loo/reference/extract_log_lik.html

###########################################
# data analysis with fake data
###########################################
# run the bayesian model
res2 <- bsem(model, data=as.data.frame(x_fake),# target="jags",
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1",std.lv = TRUE) # this is wrong
res2 #summary(res2), blavInspect(res2, 'neff'), blavInspect(res2, 'rhat')
'
  Statistic                                 MargLogLik         PPP
  Value                                             NA       0.000
'
blavFitIndices(res2)

'
      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.208        0.907        0.726        0.857 
'

#+++++++++++++++++++++++++++++++++
fit1_FAKE_INDICIES <- bcfa(model, data = as.data.frame(x_fake), 
                           n.chains = 4, burnin = 1000, sample = 1000)
null.model <- c(paste0("V", 1:6, " ~~ V", 1:6), paste0("V", 1:6, " ~ 1"))
fit0_FAKE_INDICIES <- bcfa(null.model, data = as.data.frame(x_fake),
                           n.chains = 4, burnin = 1000, sample = 1000)

blavFitIndices(fit1_FAKE_INDICIES, baseline.model = fit0_FAKE_INDICIES)

'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc         BCFI         BTLI         BNFI 
       0.210        0.907        0.721        0.857        0.820        0.760        0.808 
Warning message:
'
#+++++++++++++++++++++++++++++++++

# 2nd model: jags
data.jags <- list(N=N,           #sample size
                  y=x_fake,         #data
                  mu.0=c(0,0),   #zero-vector
                  phi.0=phi) #diagonal matrix

params <- c("ly","sigma.eps","sigma.eta")
model_real_jags <- jags.parallel(data=data.jags, 
                                 parameters.to.save=params,
                                 n.iter=2000, n.chains=4,n.thin=1,n.burnin = 1000,
                                 model.file="model1_cfa.txt")
model_real_jags
# 3rd model: stan



'
Questions to ask:
1) How to make BRMSEA  non NaN?
2) Am I doing it correctly so far?
3) # fitmeasures(res2) is not working anymore? 
'

##############################################
### Bayesian RMSEA from Rens
### MGV: modified to also calculate gammahat, adjusted gammahat
### TDJ: added McDonald's centrality index

### Bayesian RMSEA from Rens
### MGV: modified to also calculate gammahat, adjusted gammahat
### TDJ: added McDonald's centrality index
### if a null model information provided, it calculates CFI, TLI, NFI
BayesRelFit <-function(obs, rep, nvar, pD, N, ms = TRUE, Min1 = TRUE, Ngr = 1,
                       null_model=TRUE, obs_null=NULL, rep_null=NULL, pD_null=NULL){
  
  # # Compute number of parameters
  if(ms) p <- (((nvar * (nvar + 1)) / 2) + nvar)
  if(!ms) p <- (((nvar * (nvar + 1))/ 2) + 0)
  p <- p * Ngr
  # # Substract parameters and estimated parameters
  dif.ppD <- p - pD
  nonc <- ( ( obs-pD ) - dif.ppD )
  # # Correct if numerator is smaller than zero
  nonc[nonc < 0] <- 0
  # # Compute BRMSEA (with or without the -1 correction)
  if(Min1)
    BRMSEA <- sqrt(nonc / (dif.ppD * (N -1)))*sqrt(Ngr)
  if(!Min1) BRMSEA <- sqrt(nonc / (dif.ppD * N ))*sqrt(Ngr)
  
  ## compute GammaHat and adjusted GammaHat
  gammahat <- nvar / ( nvar+2* ((nonc)/(N-1))  )
  adjgammahat <- 1 - (((Ngr * nvar * (nvar + 1))/2)/dif.ppD) * (1 - gammahat)
  
  ## compute McDonald's centrality index
  Mc <- exp(-.5 * nonc/(N-1) )
  
  ## calculate fit when null model is provided
  if (null_model) {
    dif.ppD_null <- p - pD_null
    nonc_null <- ( ( obs_null-pD_null ) - dif.ppD_null )
    
    cfi <- (nonc_null - nonc)/nonc_null 
    tli <- ((( obs_null-pD_null )/dif.ppD_null) - (( obs-pD )/dif.ppD)) / (((( obs_null-pD_null )/dif.ppD_null))-1) 
    nfi <- (( obs_null-pD_null ) - ( obs-pD )) / ( obs_null-pD_null )
    
    out <- cbind(BRMSEA=BRMSEA, BGammaHat=gammahat, adjBGammaHat=adjgammahat, BMc = Mc, BCFI=cfi, BTLI=tli, BNFI=nfi)
  } else {
    out <- cbind(BRMSEA=BRMSEA, BGammaHat=gammahat, adjBGammaHat=adjgammahat, BMc = Mc)
  }
  
  return(out)
}

