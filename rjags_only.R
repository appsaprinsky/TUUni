# install.packages("rstan")
install.packages("R2jags")
install.packages("rjags")
# install.packages("blavaan")
# install.packages("mvtnorm")

library(rjags) # issue with jags instalation
library(R2jags)
# library(blavaan) # standard package
# library(rstan)
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

###########################################
# data analysis with fake data
###########################################


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



