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

###########################################
# generate real data
###########################################
set.seed(1024)
N <- 250
phi <- diag(2)
eta <- rmvnorm(N,sigma=phi)
epsilon <- rmvnorm(N,sigma=diag(6))
lambda <- matrix(c(1,1,1,0,0,0,
                   1,0,0,1,1,1),6,2,byrow=F) # misspecify lambda -> 1 where a zero is (one of them)
x <- eta%*%t(lambda)+epsilon #observed data
x <- scale(x)
target_correlation_1 <- 0.6
target_correlation_2 <- 0.6

###########################################
# generate real data new DATA
###########################################

# original_series <- rnorm(1000)
# target_correlation <- 0.6
# independent_series <- rnorm(1000)
# new_series <- target_correlation * original_series + sqrt(1 - target_correlation^2) * independent_series
# cor(original_series, new_series)

Tt <- 5 # num of Series
rows <- dim(epsilon)[1]
cols <- dim(epsilon)[2]
depth <- Tt
matrix_X <- array(0, dim = c(rows, cols, depth))
matrix_ETA <- array(0, dim = c(rows, dim(phi)[2], depth))
for (tt in 1:Tt){
  eta <- rmvnorm(N,sigma=phi)
  matrix_ETA[,,tt] <- eta[1:rows, 1:dim(phi)[2]]
  if (tt != 1){
    eta_1a <- target_correlation_1 * matrix_ETA[,,tt-1][,1] + sqrt(1 - target_correlation_1^2) * eta[,1]
    eta_2a <- target_correlation_2 * matrix_ETA[,,tt-1][,2] + sqrt(1 - target_correlation_2^2) * eta[,2]
    # print(cor(matrix_ETA[,,tt-1], matrix_ETA[,,tt]))
    eta[,1] <- eta_1a
    eta[,2] <- eta_2a
    print(cor(matrix_ETA[,,tt-1],eta))
    matrix_ETA[,,tt] <- eta[1:rows, 1:dim(phi)[2]]
  }
  epsilon <- rmvnorm(N,sigma=diag(6))
  x <- eta%*%t(lambda)+epsilon
  x <- scale(x)  
  matrix_X[,,tt] <- x[1:rows, 1:cols]
}

cor(matrix_ETA[,,1], matrix_ETA[,,3])

summary(matrix_ETA[,,5])


mu_0 <- array(0, dim = c(2))
x_mean <-array(0, dim = c(6, Tt))
cov_x_DATA <- array(0, dim = c(Tt, 6, 6)) 
dim(mu_0)
  
for (tt in 1:Tt){
  x_mean[,tt] <- colMeans(matrix_X[,,tt])
  cov_x_DATA[tt,,] <- cov(matrix_X[,,tt])
  
}

# dim(matrix_X) <- c(1000, 6, 15)
data2 <- list(
  Tt=Tt,
  N=N,
  x=matrix_X,
  mu_0=mu_0, 
  x_mean=x_mean,
  cov_x_DATA=cov_x_DATA
)

dim(matrix_X)


getwd()
setwd("/Users/admin/Desktop/Research Project/CODE/TUUni")
rt1 <- stanc("current_model_short_HB.stan") #proper & high      (b1pop,.1)
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)


fit1 <- sampling(sm1, data=data2,chains=1,iter=1500,warmup=500)

stan_rhat(fit1,bins=60)
posterior_samples <- extract(fit1)
ddd <- data.frame(posterior_samples$log_lik)

################### Measures for STAN #####################
dim(posterior_samples$log_lik)

big_my_array_1 <- array(0, dim = c(1000*Tt, 1, 2))
for (i in 1:Tt){
  print(i)
  print(1000*(i-1)+1)
  print(1000*i)
  my_array1 <- Create2DimMatrix(posterior_samples$log_lik[,,i], posterior_samples$log_lik_sat[,,i], 1000)
  big_my_array_1[(1000*(i-1)+1):(1000*i),,] <- my_array1
}

big_my_array_1

fit.measures = "all"

chisqs1 <- as.numeric(apply(big_my_array_1, 2,
                            function(x) 2*(x[,2] - x[,1]))) 
a1 <- BayesChiFit(obs = chisqs1, 
                  nvar = 6, pD = loo(fit1)$p_loo[1],
                  N = 250,
                  fit.measures = fit.measures, ms = FALSE, null_model = FALSE)

'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
         NaN        0.907        1.039        0.857 
'

out2 <- c(mean(a1@indices$BRMSEA),
          mean(a1@indices$BGammaHat),
          mean(a1@indices$adjBGammaHat),
          mean(a1@indices$BMc))






big_my_array_1[4500:5000,,]



























########################################### ###################################
########################################### ###################################
########################################### ###################################
########################################### ###################################
########################################### ###################################
# BLAVAN COMPARISON
########################################### ###################################
########################################### ###################################
########################################### ###################################
########################################### ###################################
########################################### ###################################

 ###########################################
# data analysis with real data
###########################################
# 1st model 
model <- '
eta1 =~ V1+V2+V3
eta2 =~ V4+V5+V6
'

res1 <- bsem(model, data=as.data.frame(x),#target="jags",    
             burnin = 500, n.chains = 1, sample = 1000, 
             mcmcfile = "model1", std.lv = TRUE,
             dp=dpriors(lambda="normal(0,1)"))#,itheta="dt(0,1/5,1)[sd]"))
res1 
summary(res1)
blavFitIndices(res1)


#?dpriors

'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.017        0.998        0.994        0.996 
'

#########################################################
#########################################################

###########################################################
# STAN

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
getwd()
setwd("/Users/admin/Desktop/Research Project/CODE/TUUni")
rt1 <- stanc("current_model_short_HB.stan") #proper & high      (b1pop,.1)
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)


fit1 <- sampling(sm1, data=data2,chains=1,iter=1500,warmup=500)

stan_rhat(fit1,bins=60)



params <- c("epsilon","sigma","ly")
parlog <- c("log_lik")


print(fit1,pars=params)
print(fit1,pars=parlog)
print(fit1,pars="curr_cov_ll")
print(fit1,pars="mu_together_ml")


#stan_trace(fit1,"b1")
#stan_dens(fit1,"b1",separate_chains =T)

#extract parameters for model-implied covariance matrix
par0 <- c(paste0("curr_cov_ll[",1:6,",",1,"]"),
          paste0("curr_cov_ll[",1:6,",",2,"]"),
          paste0("curr_cov_ll[",1:6,",",3,"]"),
          paste0("curr_cov_ll[",1:6,",",4,"]"),
          paste0("curr_cov_ll[",1:6,",",5,"]"),
          paste0("curr_cov_ll[",1:6,",",6,"]"))

a2 <- as.data.frame(fit1)
cov.stan <- matrix(apply(a2[,par0],2,mean),6,6)
cov.emp <- cov(x)
cov.blav <- res1@Fit@Sigma.hat[[1]]

round(cov.emp-cov.stan,3)
round(cov.blav-cov.stan,3)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#model_real_stan <- stan(model_code = stan_model_code, iter = 1500, verbose = FALSE, data=stan_data, chains=1, warmup = 500 )
posterior_samples <- extract(fit1)

ddd <- data.frame(posterior_samples$log_lik)
sum(ddd[1,])
ll2 <- mean(rowSums (ddd))
ll1 <- mean(res1@external$samplls[,,1])

#qqplot(rowSums(dddalt),rowSums(ddd))


qqplot(rowSums(ddd), res1@external$samplls[,,1])
abline(0,1)

plot(density(rowSums(ddd)),xlim=c(min(rowSums(ddd)),max(rowSums(ddd))),ylim=c(0,.2))
par(new=T)
plot(density(res1@external$samplls[,,1]),xlim=c(min(rowSums(ddd)),max(rowSums(ddd))),ylim=c(0,.2),lty=2)
abline(v=c(mean(rowSums(ddd)),mean(res1@external$samplls[,,1])))



ddd1 <- data.frame(posterior_samples$log_lik_sat)
sum(ddd1[1,])
rowSums (ddd1)


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#"
#Please run first # IMPORTANT PRE USED FUNCTIONS
#
#Create2DimMatrix
#BayesChiFit
#"



################### Measures for BLAVAAN #####################
pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"


chisqs <- as.numeric(apply(res1@external$samplls, 2,
                           function(x) 2*(x[,2] - x[,1])))   
fit_pd <- fitMeasures(res1, paste0('p_', pD))              
#fitMeasures(res1)
#'
#      npar       logl        ppp        bic        dic      p_dic       waic     p_waic    se_waic      looic 
#    13.000  -2495.968      0.541   5063.664   5017.958     13.011   5018.188     13.005     56.695   5018.271 
#     p_loo     se_loo margloglik 
#    13.047     56.699         NA 
#'

a0 <- BayesChiFit(obs = chisqs, 
            nvar = 6, pD = fit_pd[1],
            N = 250,
            fit.measures = fit.measures, ms = FALSE, null_model = FALSE)#@details


out1 <- c(mean(a0@indices$BRMSEA),
          mean(a0@indices$BGammaHat),
          mean(a0@indices$adjBGammaHat),
          mean(a0@indices$BMc))


#'
#Posterior mean (EAP) of devm-based fit indices:
#
#      BRMSEA    BGammaHat adjBGammaHat          BMc 
#       0.017        0.998        0.994        0.996 
#'

########################################

################### Measures for STAN #####################
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


#'
#Posterior mean (EAP) of devm-based fit indices:
#
#      BRMSEA    BGammaHat adjBGammaHat          BMc 
#       0.021        0.996        0.994        0.995 
#'

################### COMPARE#####################
qqplot(chisqs, chisqs1)
abline(0,1)
########################################

#######################################################################

# IMPORTANT PRE USED FUNCTIONS

#######################################################################


BayesChiFit <- function(obs, reps = NULL, nvar, pD, N, Ngr = 1,
                        ms = TRUE, Min1 = FALSE,
                        rescale = c("devM","ppmc"), fit.measures = "all",
                        null_model = TRUE, obs_null = NULL,
                        reps_null = NULL, pD_null = NULL) {
  if (!is.character(fit.measures)) {
    stop('blavaan ERROR: fit.measures must be a character vector')
  }
  fit.measures <- tolower(fit.measures)
  if (any(fit.measures == "all")) {
    fit.measures <- c("brmsea","bgammahat","adjbgammahat","bmc")
    if (null_model) fit.measures <- c(fit.measures, "bcfi","btli","bnfi")
  }
  
  if (Min1) N <- N - Ngr
  
  rescale <- tolower(as.character(rescale[1]))
  if (rescale == "devm") {
    reps <- pD
    if (!is.null(null_model)) reps_null <- pD_null
  }
  if (rescale == "ppmc" && (is.null(reps) || (null_model && is.null(reps_null)))) {
    stop('blavaan ERROR: rescale="ppmc" requires non-NULL reps argument (and reps_null, if applicable).')
  }
  
  ## ensure number of variables is a vector with length == Ngr
  #FIXME: This shouldn't be necessary.  
  #       Is object@Model@nvar always a vector, even when equal across groups?
  if (Ngr > 1L) {
    if (length(nvar) == 1L) nvar <- rep(nvar, Ngr)
  }
  ## Compute number of modeled moments
  p <- sum(sapply(nvar, function(nv) {
    nMoments <- nv * (nv + 1) / 2      # sample (co) variances
    if (ms) nMoments <- nMoments + nv  # plus means
    nMoments
  } ))
  ## Difference between number of moments and effective number of parameters
  dif.ppD <- p - pD
  
  if (dif.ppD[1] < 0) warning("blavaan WARNING: The effective number of parameters exceeds the number of sample statistics (covariances, etc.), so fit index calculations may lead to uninterpretable results.", call. = FALSE)
  
  nonc <- obs - reps - dif.ppD # == obs - p when rescale == "devm" because reps = pD
  ## Correct if numerator is smaller than zero
  nonc[nonc < 0] <- 0
  
  ## assemble results in a vector
  result <- list()
  
  ## Compute BRMSEA
  if ("brmsea" %in% fit.measures) {
    result[["BRMSEA"]] <- sqrt(nonc / (dif.ppD * N)) * sqrt(Ngr)
  }
  
  ## compute GammaHat and adjusted GammaHat
  if ("bgammahat" %in% fit.measures) {
    result[["BGammaHat"]] <- sum(nvar) / (sum(nvar) + 2*nonc/N)
    if ("adjbgammahat" %in% fit.measures) {
      result[["adjBGammaHat"]] <- 1 - (p / dif.ppD) * (1 - result[["BGammaHat"]])
    }
  } else if ("adjbgammahat" %in% fit.measures) {
    gammahat <- sum(nvar) / (sum(nvar) + 2*nonc/N)
    result[["adjBGammaHat"]] <- 1 - (p / dif.ppD) * (1 - gammahat)
  }
  
  ## compute McDonald's centrality index
  if ("bmc" %in% fit.measures) {
    result[["BMc"]] <- exp(-.5 * nonc/N)
  }
  
  ## calculate incremental fit when null model is provided
  if (null_model) {
    dif.ppD_null <- p - pD_null
    nonc_null <- (obs_null - reps_null) - dif.ppD_null
    
    if ("bcfi" %in% fit.measures) {
      result[["BCFI"]] <- 1 - (nonc / nonc_null)
    }
    if ("btli" %in% fit.measures) {
      tli_null_part <- (obs_null - reps_null) / dif.ppD_null
      result[["BTLI"]] <- (tli_null_part - (obs - reps) / dif.ppD) / (tli_null_part - 1)
    }
    if ("bnfi" %in% fit.measures) {
      result[["BNFI"]] <- ((obs_null - reps_null) - (obs - reps)) / (obs_null - reps_null)
    }
  }
  
  out <- new("blavFitIndices",
             details = list(chisq = obs - reps, df = dif.ppD,
                            pD = pD, rescale = rescale),
             indices = result)
  ## for print methods
  class(out@details$chisq) <- c("lavaan.vector","numeric")
  class(out@details$df) <- c("lavaan.vector","numeric")
  for (i in seq_along(out@indices)) {
    class(out@indices[[i]]) <- c("lavaan.vector","numeric")
  }
  
  out
}






Create2DimMatrix <- function(object_ll, object1_ll_sat, dim=1000){
  
  mat1 <- data.frame(object_ll)
  print(sum(mat1[1,]))

  mat2 <- data.frame(object1_ll_sat)
  print(sum(mat2[1,]))

  my_matrix1 <- matrix(rowSums (mat1), nrow = dim, ncol = 1)
  my_array1 <- array(data = my_matrix1, dim = c(dim, 1, 2))
  my_array1[,,2] <- rowSums (mat2)
  my_array1
}




