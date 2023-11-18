# short sim

# https://discourse.mc-stan.org/t/specification-of-bayesian-sem-models-with-a-data-augmentation-approach/19208
library(loo)
library(dplyr)
library(blavaan) # standard package
library(lavaan) # standard package
library(rstan)
library(mvtnorm) # mv normal data
packageVersion("blavaan")
packageVersion("rjags")
packageVersion("rstan")

source("helperfuns.R")

###########################################
# data analysis with real data
###########################################
# 1st model 
model <- '
eta1 =~ V1+V2+V3
eta2 =~ V4+V5+V6
'

setwd("C:\\holger\\Uni Tuebingen\\lehre\\gutachten\\Saprinsky")
rt1 <- stanc("current_model_short_HB.stan") #proper & high      (b1pop,.1)
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)

pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"

R <- 50
results <- matrix(NA,R,10)

# start loop
set.seed(1024)
for(r in 1:R){
  ###########################################
  # generate real data
  ###########################################
  N <- 250
  phi <- diag(2)
  eta <- rmvnorm(N,sigma=phi)
  epsilon <- rmvnorm(N,sigma=diag(6))
  lambda <- matrix(c(1,1,1,0,0,0,
                     0,0,0,1,1,1),6,2,byrow=F) # misspecify lambda -> 1 where a zero is (one of them)
  x <- eta%*%t(lambda)+epsilon #observed data
  x <- scale(x)
  
  
  ###########################################################
  # BLAVAAN
  ###########################################################
  res1 <- bsem(model, data=as.data.frame(x),#target="jags",    
               burnin = 500, n.chains = 1, sample = 1000, 
               mcmcfile = "model1", std.lv = TRUE,
               dp=dpriors(lambda="normal(0,1)"))#,itheta="dt(0,1/5,1)[sd]"))
  
  ###########################################################
  # STAN
  ###########################################################
  data2 <- list(
    N=N,
    x=as.data.frame(x),
    mu_0=c(0,0), 
    x_mean=colMeans(x),
    cov_x_DATA=cov(x),
    identity_matrix = diag(6),
    #LY_TRIAL = c(1.01855482, 1.14085758, 1.10058595, 1.10067097, 1.06861697, 0.98243182),
    #VARIANCE_TRIAL = c(1.291, 0.859, 1.005, 1.038, 1.106, 0.875),
    zero_vector = c(0,0,0,0,0,0)
  )
  
  
  fit1 <- sampling(sm1, data=data2,chains=1,iter=1500,warmup=500)
  posterior_samples <- extract(fit1)
  
  ddd <- data.frame(posterior_samples$log_lik)
  
  ddd1 <- data.frame(posterior_samples$log_lik_sat)
  sum(ddd1[1,])
  rowSums (ddd1)
  
  
  ################### FIT Measures blavaan #####################
  # likelihood
  ll1 <- mean(res1@external$samplls[,,1])
  
  # chisq.
  chisqs <- as.numeric(apply(res1@external$samplls, 2,
                             function(x) 2*(x[,2] - x[,1])))   
  fit_pd <- fitMeasures(res1, paste0('p_', pD))              
  
  a0 <- BayesChiFit(obs = chisqs, 
                    nvar = 6, pD = fit_pd[1],
                    N = 250,
                    fit.measures = fit.measures, ms = FALSE, null_model = FALSE)#@details
  
  out1 <- c(mean(a0@indices$BRMSEA),
            mean(a0@indices$BGammaHat),
            mean(a0@indices$adjBGammaHat),
            mean(a0@indices$BMc))
  
  ################### FIT Measures for STAN #####################
  # likelihood
  ll2 <- mean(rowSums (ddd))
  
  loo(fit1)$p_loo[1]
  
  my_array1 <- Create2DimMatrix(posterior_samples$log_lik, posterior_samples$log_lik_sat, 1000)
  chisqs1 <- as.numeric(apply(my_array1, 2,
                              function(x) 2*(x[,2] - x[,1]))) 
  a1 <- BayesChiFit(obs = chisqs1, 
                    nvar = 6, pD = loo(fit1)$p_loo[1],
                    N = 250,
                    fit.measures = fit.measures, ms = FALSE, null_model = FALSE)
  
  out2 <- c(mean(a1@indices$BRMSEA),
            mean(a1@indices$BGammaHat),
            mean(a1@indices$adjBGammaHat),
            mean(a1@indices$BMc))
  
  
  results[r,] <- c(ll1,ll2,out1,out2)
}

colnames(results) <- c("logB","logS",paste0("fitsB",1:4),paste0("fitsS",1:4))
results <- data.frame(results)
matrix(apply(results[,-(1:2)],2,mean),nrow=2,byrow=T)

plot(results$logB,results$logS)
plot(density(results$logB))
par(new=T)
plot(density(results$logS))





