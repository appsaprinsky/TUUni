###########################################
# HB
###########################################
# https://discourse.mc-stan.org/t/specification-of-bayesian-sem-models-with-a-data-augmentation-approach/19208
setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/MultitimeVSoneHB")
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
N <- 250 # 250
phi <- diag(2)
eta <- rmvnorm(N,sigma=phi)
epsilon <- rmvnorm(N,sigma=diag(6))
lambda <- matrix(c(1,1,1,0,0,0,
                   1,0,0,1,1,1),6,2,byrow=F) # misspecify lambda -> 1 where a zero is (one of them)
x <- eta%*%t(lambda)+epsilon #observed data
# x <- scale(x)
target_correlation_1 <- 0.6
target_correlation_2 <- 0.6



Tt <- 1 # num of Series
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
  # x <- scale(x)  
  matrix_X[,,tt] <- x[1:rows, 1:cols]
}




mu_0 <- array(0, dim = c(2))
x_mean <-array(0, dim = c(6, Tt))
cov_x_DATA <- array(0, dim = c(Tt, 6, 6)) 
dim(mu_0)

for (tt in 1:Tt){
  x_mean[,tt] <- colMeans(matrix_X[,,tt])
  cov_x_DATA[tt,,] <- cov(matrix_X[,,tt])
  
}

''
cov(x)
cov_x_DATA[,,]
''

# dim(matrix_X) <- c(1000, 6, 15)
data2_N <- list(
  Tt=Tt,
  N=N,
  x=matrix_X,
  mu_0=mu_0, 
  x_mean=x_mean,
  cov_x_DATA=cov_x_DATA
)

dim(matrix_X)
getwd()
rt1_N <- stanc("Multivariate.stan") #proper & high      (b1pop,.1)
sm1_N <- stan_model(stanc_ret = rt1_N, verbose=FALSE)


fit1_N <- sampling(sm1_N, data=data2_N,chains=1,iter=1500,warmup=500)

posterior_samples_N <- extract(fit1_N)

################### Measures for STAN #####################


big_my_array_1_N <- array(0, dim = c(1000*Tt, 1, 2))
for (i in 1:Tt){
  print(i)
  print(1000*(i-1)+1)
  print(1000*i)
  my_array1 <- Create2DimMatrix(posterior_samples_N$log_lik[,,i], posterior_samples_N$log_lik_sat[,,i], 1000)
  big_my_array_1_N[(1000*(i-1)+1):(1000*i),,] <- my_array1
}

'
[1] -2557.899
[1] -2512.573
'
fit.measures = "all"

chisqs1_N <- as.numeric(apply(big_my_array_1_N, 2,
                            function(x) 2*(x[,2] - x[,1]))) 
a1_N <- BayesChiFit(obs = chisqs1_N, 
                  nvar = 6, pD = loo(fit1_N)$p_loo[1],
                  N = 250,
                  fit.measures = fit.measures, ms = FALSE, null_model = FALSE)


'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.201        0.914        0.741        0.868 

'













###########################################
# generate real data
###########################################
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
rt1 <- stanc("HB.stan") #proper & high      (b1pop,.1)
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)


fit1 <- sampling(sm1, data=data2,chains=1,iter=1500,warmup=500)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#model_real_stan <- stan(model_code = stan_model_code, iter = 1500, verbose = FALSE, data=stan_data, chains=1, warmup = 500 )
posterior_samples <- extract(fit1)


################### Measures for BLAVAAN #####################
pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"


########################################

################### Measures for STAN #####################
loo(fit1)$p_loo[1]

my_array1 <- Create2DimMatrix(posterior_samples$log_lik, posterior_samples$log_lik_sat, 1000)
'
[1] -2557.808
[1] -2512.573
'
chisqs1 <- as.numeric(apply(my_array1, 2,
                           function(x) 2*(x[,2] - x[,1]))) 
a1 <- BayesChiFit(obs = chisqs1, 
            nvar = 6, pD = loo(fit1)$p_loo[1],
            N = 250,
            fit.measures = fit.measures, ms = FALSE, null_model = FALSE)
'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.206        0.914        0.729        0.868 
'







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




