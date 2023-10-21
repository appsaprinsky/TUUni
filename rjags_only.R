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
                  phi=phi) #diagonal matrix

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
# data analysis with fake data
###########################################

# 2nd model: jags
data.jags <- list(N=N,           #sample size
                  y=x_fake,         #data
                  mu.0=c(0,0),   #zero-vector
                  phi=phi) #diagonal matrix

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
                  phi=phi) #diagonal matrix

params <- c("ly","sigma.eps","sigma.eta")

###############################################################
# RUN MODELS
###############################################################
model1 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="model1_dsem.txt")

data.jags <- list(N=N,           #sample size
                  y=x_fake,         #data
                  mu.0=c(0,0),   #zero-vector
                  phi=phi) #diagonal matrix


model2 <- jags.parallel(data=data.jags, 
                        parameters.to.save=params,
                        n.iter=4000, n.chains=3,n.thin=1,n.burnin = 2000,
                        model.file="model1_dsem.txt")

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
round(est2,2)


# extract likelihoods
niter <- 2000
nchain <- 3
mc_ll1 <- mc_ll2 <- array(NA,c(niter,nchain,N))

samps1 <-  as.mcmc(model1)
samps2 <-  as.mcmc(model2)


for(r in 1:nchain){
  mc_ll1[,r,] <-  samps1[[r]][,paste0("loglik1[",1:N,"]")]
  mc_ll2[,r,] <-  samps2[[r]][,paste0("loglik1[",1:N,"]")]
}

# model fit indices are similar to information criteria (AIC, BIC)
waic1 <- waic(mc_ll1)
waic2 <- waic(mc_ll2)
loo1 <- loo(mc_ll1)
loo2 <- loo(mc_ll2)

# WAIC
print(loo_compare(waic1, waic2), digits = 2)
# model 3 is fitting best

# loo
print(loo_compare(loo1, loo2), digits = 2)
# same here

# DIC
model1$BUGSoutput$DIC
model2$BUGSoutput$DIC









