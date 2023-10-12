library(mvtnorm) # mv normal data
library(blavaan) # standard package
# install.packages("blavaan")

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
eta1 =~ V1+V2+V3
eta2 =~ V4+V5+V6
'


# run the bayesian model
res1 <- bsem(model, data=as.data.frame(x),target="jags",
             burnin = 1000, n.chains = 4, sample = 1000,
             mcmcfile = "model1",
             std.lv = TRUE) # this is wrong

# extract model fit

res1
# blavaan 0.5.1 ended normally after 1000 iterations
# 
# Estimator                                      BAYES
# Optimization method                             MCMC
# Number of model parameters                        13
# 
# Number of observations                           250
# 
# Statistic                                 MargLogLik         PPP
# Value                                      -2495.627       0.613
summary(res1)
# blavaan 0.5.1 ended normally after 1000 iterations
# 
# Estimator                                      BAYES
# Optimization method                             MCMC
# Number of model parameters                        13
# 
# Number of observations                           250
# 
# Statistic                                 MargLogLik         PPP
# Value                                      -2495.627       0.613
# 
# Parameter Estimates:
#   
#   
#   Latent Variables:
#   Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
# eta1 =~                                                                      
#   V1                1.088    0.090    0.912    1.270    0.999    normal(0,10)
# V2                1.003    0.091    0.835    1.178    1.000    normal(0,10)
# V3                1.102    0.093    0.925    1.289    1.000    normal(0,10)
# eta2 =~                                                                      
#   V4                0.882    0.101    0.684    1.078    0.999    normal(0,10)
# V5                1.095    0.111    0.877    1.318    0.999    normal(0,10)
# V6                0.880    0.099    0.683    1.070    0.999    normal(0,10)
# 
# Covariances:
#   Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
# eta1 ~~                                                                      
#   eta2             -0.055    0.085   -0.219    0.107    0.999     lkj_corr(1)
# 
# Variances:
#   Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
# .V1                0.785    0.131    0.537    1.052    1.000 gamma(1,.5)[sd]
# .V2                1.006    0.131    0.770    1.284    1.000 gamma(1,.5)[sd]
# .V3                0.921    0.138    0.661    1.198    0.999 gamma(1,.5)[sd]
# .V4                1.191    0.151    0.912    1.495    1.000 gamma(1,.5)[sd]
# .V5                0.992    0.191    0.601    1.367    0.999 gamma(1,.5)[sd]
# .V6                0.976    0.142    0.718    1.263    1.000 gamma(1,.5)[sd]
# eta1              1.000                                                    
# eta2              1.000





# 2nd model: jags



# 3rd model: stan
