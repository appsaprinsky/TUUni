install.packages("R2jags")
install.packages("mvtnorm")

library(R2jags)
library(mvtnorm) # mv normal data

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
'
Inference for Bugs model at "model1_cfa.txt", fit using jags,
 4 chains, each with 2000 iterations (first 1000 discarded)
 n.sims = 4000 iterations saved
                mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat
ly[1]             1.139   0.130    0.903    1.046    1.136    1.228    1.395 1.014
ly[2]             1.100   0.132    0.866    1.008    1.091    1.185    1.377 1.007
ly[3]             0.982   0.112    0.771    0.906    0.978    1.060    1.203 1.005
ly[4]             0.909   0.111    0.706    0.832    0.904    0.982    1.141 1.009
sigma.eps[1]      1.280   0.152    1.000    1.178    1.273    1.377    1.600 1.004
sigma.eps[2]      0.852   0.152    0.571    0.751    0.846    0.945    1.164 1.009
sigma.eps[3]      0.994   0.154    0.701    0.891    0.991    1.093    1.305 1.006
sigma.eps[4]      1.025   0.157    0.732    0.919    1.016    1.125    1.353 1.005
sigma.eps[5]      1.114   0.153    0.832    1.010    1.108    1.212    1.425 1.001
sigma.eps[6]      0.863   0.127    0.621    0.778    0.860    0.944    1.113 1.002
sigma.eta[1,1]    1.000   0.192    0.669    0.862    0.989    1.121    1.406 1.013
sigma.eta[2,1]   -0.064   0.091   -0.246   -0.123   -0.062   -0.004    0.112 1.006
sigma.eta[1,2]   -0.064   0.091   -0.246   -0.123   -0.062   -0.004    0.112 1.006
sigma.eta[2,2]    1.186   0.212    0.816    1.038    1.166    1.318    1.652 1.010
deviance       4258.208  52.112 4179.329 4227.740 4255.619 4283.727 4345.906 1.001
               n.eff
ly[1]            200
ly[2]            490
ly[3]            520
ly[4]            310
sigma.eps[1]     670
sigma.eps[2]     330
sigma.eps[3]     760
sigma.eps[4]     580
sigma.eps[5]    4000
sigma.eps[6]    1800
sigma.eta[1,1]   240
sigma.eta[2,1]   480
sigma.eta[1,2]   480
sigma.eta[2,2]   310
deviance        4000

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1358.2 and DIC = 5616.4
DIC is an estimate of expected predictive error (lower deviance is better).
'

###########################################
# data analysis with fake data
###########################################

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
'
For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1358.2 and DIC = 5616.4
DIC is an estimate of expected predictive error (lower deviance is better).
'



###########################################
# STEP 2
###########################################
library(R2jags)
install.packages("loo")
library(loo)

###############################################################
# DATA
# simulated data set based on Flückiger et al. (2022)
###############################################################
data.jags <- list(N=N,           #sample size
                  y=x,         #data
                  mu.0=c(0,0),   #zero-vector
                  phi.0=phi) #diagonal matrix

params <- c("ly","sigma.eps","sigma.eta")

###############################################################
# RUN MODELS
###############################################################
model1 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="model1_cfa.txt")

data.jags <- list(N=N,           #sample size
                  y=x_fake,         #data
                  mu.0=c(0,0),   #zero-vector
                  phi.0=phi) #diagonal matrix


model2 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="model1_cfa.txt")

###############################################################
samps <- as.mcmc(model1)
#xyplot, densityplot
pdf(file="xyplot.pdf",pointsize=6)
plot(samps)
dev.off()

pdf("gelmanplot.pdf")
par(new=F)
for (v in 1:nvar(samps)){
  gelman.plot(samps[,v],ylim=c(0,7),col=c("black",NA))
  par(new=T)
}
abline(h=1.1)
dev.off()
###############################################################

est1 <- model1$BUGSoutput$summary
est2 <- model2$BUGSoutput$summary
round(est1,2)
'
                  mean    sd    2.5%     25%     50%     75%   97.5% Rhat n.eff
deviance       4408.32 49.81 4318.44 4376.90 4407.37 4437.00 4500.67 1.00  6000
ly[1]             1.14  0.18    0.83    1.01    1.12    1.25    1.53 1.01   400
ly[2]             1.18  0.19    0.86    1.05    1.16    1.29    1.59 1.01   270
ly[3]             1.00  0.12    0.78    0.91    0.99    1.08    1.25 1.00  1300
ly[4]             0.92  0.11    0.71    0.84    0.91    0.99    1.15 1.00   540
sigma.eps[1]      2.23  0.25    1.76    2.07    2.22    2.39    2.73 1.00  1200
sigma.eps[2]      0.97  0.18    0.62    0.85    0.96    1.08    1.31 1.00  1900
sigma.eps[3]      0.93  0.19    0.55    0.80    0.93    1.05    1.29 1.00  1200
sigma.eps[4]      1.04  0.15    0.76    0.94    1.03    1.15    1.35 1.00   920
sigma.eps[5]      1.10  0.16    0.80    0.99    1.09    1.20    1.41 1.00  4000
sigma.eps[6]      0.86  0.13    0.62    0.78    0.86    0.95    1.12 1.00  1400
sigma.eta[1,1]    0.94  0.24    0.53    0.76    0.92    1.08    1.48 1.01   360
sigma.eta[2,1]    0.10  0.10   -0.08    0.03    0.09    0.16    0.32 1.01   240
sigma.eta[1,2]    0.10  0.10   -0.08    0.03    0.09    0.16    0.32 1.01   240
sigma.eta[2,2]    1.17  0.22    0.79    1.02    1.16    1.31    1.64 1.00   660

'
round(est2,2)

'
                  mean    sd    2.5%     25%     50%     75%   97.5% Rhat n.eff
deviance       4408.32 49.81 4318.44 4376.90 4407.37 4437.00 4500.67 1.00  6000
ly[1]             1.14  0.18    0.83    1.01    1.12    1.25    1.53 1.01   400
ly[2]             1.18  0.19    0.86    1.05    1.16    1.29    1.59 1.01   270
ly[3]             1.00  0.12    0.78    0.91    0.99    1.08    1.25 1.00  1300
ly[4]             0.92  0.11    0.71    0.84    0.91    0.99    1.15 1.00   540
sigma.eps[1]      2.23  0.25    1.76    2.07    2.22    2.39    2.73 1.00  1200
sigma.eps[2]      0.97  0.18    0.62    0.85    0.96    1.08    1.31 1.00  1900
sigma.eps[3]      0.93  0.19    0.55    0.80    0.93    1.05    1.29 1.00  1200
sigma.eps[4]      1.04  0.15    0.76    0.94    1.03    1.15    1.35 1.00   920
sigma.eps[5]      1.10  0.16    0.80    0.99    1.09    1.20    1.41 1.00  4000
sigma.eps[6]      0.86  0.13    0.62    0.78    0.86    0.95    1.12 1.00  1400
sigma.eta[1,1]    0.94  0.24    0.53    0.76    0.92    1.08    1.48 1.01   360
sigma.eta[2,1]    0.10  0.10   -0.08    0.03    0.09    0.16    0.32 1.01   240
sigma.eta[1,2]    0.10  0.10   -0.08    0.03    0.09    0.16    0.32 1.01   240
sigma.eta[2,2]    1.17  0.22    0.79    1.02    1.16    1.31    1.64 1.00   660
'

# extract likelihoods
niter <- 2000
nchain <- 3
mc_ll1 <- mc_ll2 <- array(NA,c(niter,nchain,4))

samps1 <-  as.mcmc(model1)
samps2 <-  as.mcmc(model2)


for(r in 1:nchain){
  mc_ll1[,r,] <-  samps1[[r]][,paste0("ly[",1:4,"]")] # samps1[[r]][,paste0("ly[",1:4,"]")]
  mc_ll2[,r,] <-  samps2[[r]][,paste0("ly[",1:4,"]")]
}

# model fit indices are similar to information criteria (AIC, BIC)
waic1 <- waic(mc_ll1)
waic2 <- waic(mc_ll2)
loo1 <- loo(mc_ll1)
loo2 <- loo(mc_ll2)

# WAIC
print(loo_compare(waic1, waic2), digits = 2)
# model 3 is fitting best
'
       elpd_diff se_diff
model2  0.00      0.00  
model1 -0.08      0.07  

'

# loo
print(loo_compare(loo1, loo2), digits = 2)
# same here

# DIC
model1$BUGSoutput$DIC
model2$BUGSoutput$DIC









