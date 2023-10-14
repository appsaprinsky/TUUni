library(mvtnorm) # mv normal data
library(blavaan) # standard package
# install.packages("blavaan")

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
res1 <- bsem(model, data=as.data.frame(x),target="jags",
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1")
             # ,std.lv = TRUE) # this is wrong
# extract model fit
res1
'
  Statistic                                 MargLogLik         PPP
  Value                                      -2509.620       0.639
'
summary(res1)
blavFitIndices(res1)
'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.011        0.999        0.996        0.998 
'

# 2nd model: jags
# 3rd model: stan
###########################################
# data analysis with fake data
###########################################

# run the bayesian model
res2 <- bsem(model, data=as.data.frame(x_fake),target="jags",
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1")
# ,std.lv = TRUE) # this is wrong
# extract model fit
res2
'
  Statistic                                 MargLogLik         PPP
  Value                                      -3044.599       0.000
'
summary(res2)

blavFitIndices(res2)

'BRMSEA    BGammaHat adjBGammaHat          BMc 
NaN        0.896        1.091        0.840 '

# 2nd model: jags


# 3rd model: stan





'
Questions to ask:
1) How to make BRMSEA  non NaN?
2) Are four indicies enough for the project?
3) Am I doing it correctly so far?


What to do for Step 1 Project:
1) Garnier Paper read again
2) Second Paper
3) 01_cfa go through
4) Creat account for ovrleaf

''



