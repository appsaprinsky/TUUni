library(mvtnorm) # mv normal data
library(blavaan) # standard package
#install.packages("blavaan")

###########################################
# generate data
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
# data analysis
###########################################
# 1st model 
model <- '
eta1 =~ x1+x2+x3
eta2 =~ x4+x5+x6
'

# run the bayesian model
res1 <- bcfa(x,model) # this is wrong

# extract model fit



# 2nd model: jags



# 3rd model: stan
