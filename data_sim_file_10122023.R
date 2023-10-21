# install.packages("rstan")
# install.packages("R2jags")
install.packages("rjags")
# install.packages("blavaan")
# install.packages("mvtnorm")

library(rjags) # issue with jags instalation
library(R2jags)
library(blavaan) # standard package
library(rstan)
library(mvtnorm) # mv normal data
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
res1 <- bsem(model, data=as.data.frame(x),#target="jags",     ##################CORRECTLY SCPECIFIED????
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
                  phi=phi) #diagonal matrix

params <- c("ly","sigma.eps","sigma.eta")
model_real_jags <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=2000, n.chains=4,n.thin=1,n.burnin = 1000,
                        model.file="model1_cfa.txt")
model_real_jags
est1 <- model_real_jags$BUGSoutput$summary

# 3rd model: stan
stan_data <- list(
  N=N,
  x=x
)

stan_model_code <- '
data {
  int<lower=0> N;
  matrix[N, 6] x;  
}

parameters {
  vector[3] eta1;
  vector[3] eta2;
  real<lower=0> sigma;
}

model {
  eta1 ~ normal(0, 10);
  eta2 ~ normal(0, 10);
  for (i in 1:N) {
    x[i, 1:3] ~ normal(eta1, sigma);
    x[i, 4:6] ~ normal(eta2, sigma);
  }
}

generated quantities {
  
  /* define predictions and residuals */
  vector[N] y_hat;
  vector[N] resid;
  
  /* calculate predictions and residuals */
  y_hat = X * beta;
  resid = y - y_hat;
  
}

'

find.package('rstan')
model_real_stan <- stan(model_code = stan_model_code, iter = 1000, verbose = FALSE, data=stan_data ) 
posterior_samples <- extract(model_real_stan)
posterior_samples
fit_indices <- blavFitIndices(posterior_samples$eta1, posterior_samples$eta2)
posterior_samples$eta1
print(fit_indices)



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
# 3rd model: stan



'
Questions to ask:
1) How to make BRMSEA  non NaN?
3) Am I doing it correctly so far?
4) issue with jags instalation
5) cloud tempo Solution?
6) Group Partner
7) # fitmeasures(res2) is not working anymore? 
8) Stan Solution problematic? not sure if correct, several resources were used.
9) Feeling the Form

What to do for Step 1 Project:
4) Stan write the code

''



