# install.packages("rstan")
# # install.packages("R2jags")
# install.packages("blavaan")
# install.packages("mvtnorm")

# library(rjags) # issue with jags instalation
# library(R2jags)
library(blavaan) # standard package
library(rstan)
library(mvtnorm) # mv normal data
packageVersion("blavaan")
packageVersion("rjags")
packageVersion("rstan")

# BNFI
# BTLI
# BCFI
# BMc       (Found)
# BΓadj     (Found)
# BΓˆ       (Found)
# BRMSEA    (Found)
# D(θ ̄)  already in BRMSEA 
# pDLOO
# PPPχ2
###########################################
# generate real data
###########################################
# lv
N <- 250
phi <- diag(2)
eta <- rmvnorm(N,sigma=phi)
epsilon <- rmvnorm(N,sigma=diag(6))
lambda <- matrix(c(1,1,1,0,0,0,
                   0,0,0,1,1,1),6,2,byrow=F)
# misspecify lambda -> 1 where a zero is (one of them)
#observed data
x <- eta%*%t(lambda)+epsilon
###########################################
# generate fake data
###########################################
# lv
lambda_fake <- matrix(c(1,1,1,0,0,0,
                   1,0,0,1,1,1),6,2,byrow=F)
# misspecify lambda -> 1 where a zero is (one of them)
#observed data
x_fake <- eta%*%t(lambda_fake)+epsilon

###########################################
# data analysis with real data
###########################################
# 1st model 
model <- '
eta1 =~ V1+V2+V3
eta2 =~ V4+V5+V6
'
# run the bayesian model
res1 <- bsem(model, data=as.data.frame(x),#target="jags",
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1", std.lv = TRUE) # this is wrong
# extract model fit
res1
'
  Statistic                                 MargLogLik         PPP
  Value                                      -2501.918       0.558
'
summary(res1)
blavFitIndices(res1)
'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.016        0.998        0.994        0.997 
'

# 2nd model: jags
# 3rd model: stan
scode <- "
parameters {
  real y[2];
}
model {
  y[1] ~ normal(0, 1);
  y[2] ~ double_exponential(0, 2);
}
"


# scode <-"
# data {
#   int<lower=0> N; // Number of observations
#   int<lower=0> K; // Number of items (variables)
#   matrix[N, K] Y; // Observed data
# }
# 
# parameters {
#   matrix[K, K] L; // Factor loadings
#   vector[K] mu;   // Factor means
#   vector<lower=0>[K] sigma; // Residual standard deviations
# }
# 
# model {
#   // Prior distributions for the model parameters
#   L ~ normal(0, 2);  // Priors for factor loadings
#   mu ~ normal(0, 5); // Priors for factor means
#   sigma ~ cauchy(0, 5); // Priors for residual standard deviations
# 
#   // Likelihood: CFA model
#   for (n in 1:N) {
#     Y[n] ~ multi_normal(mu + L * (eta1[n] - eta2[n]), diag_matrix(sigma));
#   }
# }
# 
# generated quantities {
#   vector[K] eta1[N]; // Latent variable scores for eta1
#   vector[K] eta2[N]; // Latent variable scores for eta2
#   for (n in 1:N) {
#     eta1[n] = mu + L * eta1[n];
#     eta2[n] = mu + L * eta2[n];
#   }
# }
# "
fit1 <- stan(model_code = scode, iter = 10, verbose = FALSE)
print(fit1)
fit2 <- stan(fit = fit1, iter = 10000, verbose = FALSE)
###########################################
# data analysis with fake data
###########################################

# run the bayesian model
res2 <- bsem(model, data=as.data.frame(x_fake),# target="jags",
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1",std.lv = TRUE) # this is wrong
# extract model fit
res2
'
  Statistic                                 MargLogLik         PPP
  Value                                      -2552.209       0.000
'
summary(res2)


blavInspect(res2, 'rhat')
blavInspect(res2, 'neff')
fitMeasures(res2)

blavFitIndices(res2)

'
      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.182        0.917        0.787        0.873 
'

fitmeasures(res2)


# 2nd model: jags


# 3rd model: stan



'
Questions to ask:
1) How to make BRMSEA  non NaN?
2) Are four indicies enough for the project?
3) Am I doing it correctly so far?
4) issue with jags instalation
5) cloud tempo Solution?


What to do for Step 1 Project:
1) Garnier Paper read again
2) Second Paper
3) 01_cfa go through
4) Stan write the code

''



