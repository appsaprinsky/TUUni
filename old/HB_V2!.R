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
                   0,0,0,1,1,1),6,2,byrow=F) # misspecify lambda -> 1 where a zero is (one of them)
x <- eta%*%t(lambda)+epsilon #observed data
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
             mcmcfile = "model1", std.lv = TRUE) 
res1 
summary(res1)
blavFitIndices(res1)




'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.017        0.998        0.994        0.996 
'
res1@external$samplls[,,2]
res1@external$samplls[,,1]

-2498.066res1@external#$stansumm




stan_model_code <- '
data {
  int<lower=0> N;
  matrix[N, 6] x; 
  vector[2] mu_0;
  vector[6] x_mean;
  matrix[6, 6] cov_x_DATA;
  matrix[6, 6] identity_matrix;
  
  vector[6] LY_TRIAL;
  vector[6] VARIANCE_TRIAL;  
  // matrix[N, 2] eta;
  vector[6] zero_vector;
}

parameters {
  matrix[N, 2]  eta;
  vector<lower=0>[6] epsilon; 
  real sigma; 
  vector[6] ly;
}

model {

  matrix[2,2] COV_MATRIX;  
  vector[6] TEMP_MEAN;
  matrix[6, 6] TEPM_VAR;    
  
  COV_MATRIX[1,1] = 1;
  COV_MATRIX[2,2] = 1;
  COV_MATRIX[1,2] = sigma;
  COV_MATRIX[2,1] = sigma;
  

  for(k in 1:6){
    ly[k] ~ normal(0, 1);  // DELETE LATERRRRRRR
  }
//  sigma ~ lkj_corr(1);
  for(k in 1:6){
    epsilon[k] ~ cauchy(0, 1);   // DELETE LATERRRRRRR gamma  0, 1
  }

  for (i in 1:N) {
  
//    x[i,1] ~ normal(ly[1]*eta[i, 1], epsilon[1]);
//    x[i,2] ~ normal(ly[2]*eta[i, 1], epsilon[2]);
//    x[i,3] ~ normal(ly[3]*eta[i, 1], epsilon[3]);
//    x[i,4] ~ normal(ly[4]*eta[i, 2], epsilon[4]);
//    x[i,5] ~ normal(ly[5]*eta[i, 2], epsilon[5]);
//    x[i,6] ~ normal(ly[6]*eta[i, 2], epsilon[6]);


    TEMP_MEAN[1] = ly[1]*eta[i, 1];
    TEMP_MEAN[2] = ly[2]*eta[i, 1];    
    TEMP_MEAN[3] = ly[3]*eta[i, 1];    
    TEMP_MEAN[4] = ly[4]*eta[i, 2];    
    TEMP_MEAN[5] = ly[5]*eta[i, 2];    
    TEMP_MEAN[6] = ly[6]*eta[i, 2];
    
    for (www1 in 1:6){
      for (www2 in 1:6){
        if (www1==www2){
          TEPM_VAR[www1, www2] = epsilon[www1];
        } else{
          TEPM_VAR[www1, www2] = 0;        
        }
      }
    }
    
    x[i,1:6] ~ multi_normal(TEMP_MEAN,TEPM_VAR);    
    eta[i, 1:2] ~ multi_normal(mu_0,COV_MATRIX);    

    
  }
  // epsilon ~ cauchy(1, 1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] log_lik_sat;
  
  matrix[6, 6] curr_cov_ll;
  
  vector[6] x_together_ml;
  vector[6] mu_together_ml;  
  
  vector[6] mu_zeros_ml;
  
  mu_together_ml[1] = 0;   
  mu_together_ml[2] = 0;  
  mu_together_ml[3] = 0;  
  mu_together_ml[4] = 0; 
  mu_together_ml[5] = 0; 
  mu_together_ml[6] = 0; 
  
  for (i in 1:N){
    mu_together_ml[1] += ly[1]*eta[i, 1];   
    mu_together_ml[2] += ly[2]*eta[i, 1];  
    mu_together_ml[3] += ly[3]*eta[i, 1];  
    mu_together_ml[4] += ly[4]*eta[i, 2]; 
    mu_together_ml[5] += ly[5]*eta[i, 2]; 
    mu_together_ml[6] += ly[6]*eta[i, 2]; 
  }
  

  
  mu_together_ml = mu_together_ml/N;
  
  for (i in 1:N) {
    log_lik[i] = 0;  
    log_lik[i] += normal_lpdf(x[i,1] | ly[1]*eta[i, 1], epsilon[1]);
    log_lik[i] += normal_lpdf(x[i,2] | ly[2]*eta[i, 1], epsilon[2]);
    log_lik[i] += normal_lpdf(x[i,3] | ly[3]*eta[i, 1], epsilon[3]);
    log_lik[i] += normal_lpdf(x[i,4] | ly[4]*eta[i, 2], epsilon[4]);
    log_lik[i] += normal_lpdf(x[i,5] | ly[5]*eta[i, 2], epsilon[5]);
    log_lik[i] += normal_lpdf(x[i,6] | ly[6]*eta[i, 2], epsilon[6]);
    
  
    for (aaa in 1:6) {
      for (sss in 1:6) {
        if (aaa == sss) {
          curr_cov_ll[aaa, sss] = ly[aaa]*ly[aaa] + epsilon[aaa]; //epsilon[aaa];  VARIANCE_TRIAL[aaa]; TRIAL
        } else {
          curr_cov_ll[aaa, sss] = ly[aaa]*ly[sss]*sigma; // curr_cov_ll[aaa, sss] = cov_x_DATA[aaa, sss]; 
        }
      }
    }
    
    for (aaa in 1:6) {
      for (sss in 1:6) {
        if (aaa == sss) {
          curr_cov_ll[aaa, sss] = ly[aaa]*ly[aaa] + VARIANCE_TRIAL[aaa]; //epsilon[aaa];  VARIANCE_TRIAL[aaa]; TRIAL
        } else {
          if ((aaa <= 3 && sss <= 3) || (aaa >= 4 && sss >= 4)){
            curr_cov_ll[aaa, sss] = ly[aaa]*ly[sss]; // curr_cov_ll[aaa, sss] = cov_x_DATA[aaa, sss]; 
          } else {
            curr_cov_ll[aaa, sss] = ly[aaa]*ly[sss]*sigma; // curr_cov_ll[aaa, sss] = cov_x_DATA[aaa, sss]; 
          }
        }
      }
    }    

    

    ////////////
    
    for (jjj in 1:3){
      x_together_ml[jjj] = x[i,jjj]; 
//      mu_together_ml[jjj] = ly[jjj]*eta[i, 1];
      mu_zeros_ml[jjj] = LY_TRIAL[jjj]*eta[i, 1];   // TRIAL
    }
    for (jjj in 4:6){
      x_together_ml[jjj] = x[i,jjj]; 
//      mu_together_ml[jjj] = ly[jjj]*eta[i, 2];
      mu_zeros_ml[jjj] = LY_TRIAL[jjj]*eta[i, 2];   // TRIAL
    }  
    

    mu_together_ml[1] = ly[1]*eta[i, 1];    
    mu_together_ml[2] = ly[2]*eta[i, 1];  
    mu_together_ml[3] = ly[3]*eta[i, 1];  
    mu_together_ml[4] = ly[4]*eta[i, 2]; 
    mu_together_ml[5] = ly[5]*eta[i, 2]; 
    mu_together_ml[6] = ly[6]*eta[i, 2]; 
    // print("COV =", curr_cov_ll);


    
    log_lik[i] = 0;  
    log_lik_sat[i] = 0;
    log_lik[i] = multi_normal_lpdf(x_together_ml | zero_vector , curr_cov_ll); // (mu_together_ml, curr_cov_ll) mu_zeros_ml
    log_lik_sat[i] = multi_normal_lpdf(x_together_ml | x_mean, cov_x_DATA);  // 
    ////////////
    
  }
}
'


stan_data <- list(
  N=N,
  x=as.data.frame(x),
  mu_0=c(0,0), 
  x_mean=colMeans(x),
  cov_x_DATA=cor(x),
  identity_matrix = diag(6)
)

COV =[[2.31736,1.16268,1.12192,-0.0605836,-0.0588228,-0.0540355],[1.16268,2.15488,1.25624,-0.067837,-0.0658654,-0.0605049],[1.12192,1.25624,2.2132,-0.0654589,-0.0635563,-0.0583838],[-0.0605836,-0.067837,-0.0654589,2.2322,1.17697,1.08118],[-0.0588228,-0.0658654,-0.0635563,1.17697,2.24476,1.04976],[-0.0540355,-0.0605049,-0.0583838,1.08118,1.04976,1.84632]]

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

stan_data <- list(
  N=N,
  x=as.data.frame(x),
  mu_0=c(0,0), 
  x_mean=colMeans(x),
  cov_x_DATA=cov(x),
  identity_matrix = diag(6),
  LY_TRIAL = c(1.026, 1.147, 1.110, 1.109, 1.078,  0.993),
  VARIANCE_TRIAL = c(1.27935035, 0.85319543, 1.00105014, 1.02030526, 1.10168394, 0.88158536),
  zero_vector = c(0,0,0,0,0,0)
)

stan_data <- list(
  N=N,
  x=as.data.frame(x),
  mu_0=c(0,0), 
  x_mean=colMeans(x),
  cov_x_DATA=cov(x),
  identity_matrix = diag(6),
  LY_TRIAL = c(1.01855482, 1.14085758, 1.10058595, 1.10067097, 1.06861697, 0.98243182),
  VARIANCE_TRIAL = c(1.291, 0.859, 1.005, 1.038, 1.106, 0.875),
  zero_vector = c(0,0,0,0,0,0)
)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


model_real_stan <- stan(model_code = stan_model_code, iter = 1500, verbose = FALSE, data=stan_data, chains=1, warmup = 500 )
posterior_samples <- extract(model_real_stan)
posterior_samples$log_lik
posterior_samples$log_lik_sat
posterior_samples$sigma
posterior_samples$cov_x_generated[1,,]
posterior_samples$epsilon[1000,]
posterior_samples$curr_cov_ll

# posterior_samples$target

ddd <- data.frame(posterior_samples$log_lik)
sum(ddd[1,])
rowSums (ddd)


ddd1 <- data.frame(posterior_samples$log_lik_sat)
sum(ddd1[1,])
rowSums (ddd1)


loo(model_real_stan)

plot(rowSums(ddd), res1@external$samplls[,,1])

ddd <- ddd %>%
  mutate(ll = rowSums(., na.rm=TRUE))

ddd$ll
res1@external$samplls[1:500,,1]
res1@external$samplls[,,1]

sum(data.frame(posterior_samples$log_lik_sat)[1,])
sum(data.frame(posterior_samples$log_lik)[1,])




res <- rstan::read_stan_csv(model_real_stan $output_files())



loo::extract_log_lik(model_real_stan, parameter_name = "log_lik_sat")


llsat <- loo::extract_log_lik(lavjags, parameter_name = "log_lik_sat")


loo::extract_log_lik(model_real_stan, parameter_name = "log_lik")
loo::extract_log_lik(res1, parameter_name = "log_lik")




posterior_samples
res1@external

res1






posterior_samples
posterior_samples$log_lik
posterior_samples$lp

posterior_samples$eta
posterior_samples$sigma
posterior_samples$epsilon
posterior_samples$mu_0
posterior_samples$mu_x
posterior_samples$ly


res1@external$samplls


summary(posterior_samples)


rstan::expose_stan_functions("model1/sem.stan")





# for(j in 1:4){
#   csdist <- rep(NA, 1000)
#   for(i in 1:1000){
#     chisq.obs <- -2*(res1@external$samplls[i, j, 1] - res1@external$samplls[i, j, 2])
#     csdist[i] <- chisq.obs
#   }
#   list(csdist = csdist)
#   csdist <- unlist(lapply(res, function(x) x$csdist))
# }
# 
# csdist
# 
# chisq.obs <- -2*(samplls[i, j, 1] -
#                    samplls[i, j, 2])

# res1 <- bcfa(model, data=as.data.frame(x),#target="jags",     ##################COnvergence rate , rhat stat 1, ess above 400
#              burnin = 1000, n.chains = 4, sample = 1000,
#              mcmcfile = "model1") # this is wrong
# postpred_mgv(4, 1000, res1@external$samplls)$chisqs[,1][1:1000]


# sampler_params <- get_sampler_params(model_real_stan, inc_warmup = FALSE)
# sampler_params







aaa <- postpred_mgv(4, 1000, res1@external$samplls)$chisqs[,1]#[1:1000]

res1@external$samplls











plot(res1@external$samplls[,,1], )


#######################################################################

res1 <- bsem(model, data=as.data.frame(x),#target="jags",    
             burnin = 500, n.chains = 1, sample = 1000, 
             mcmcfile = "model1", estimator = "ML") 
# extract model fit
res1 # summary(res1)
summary(res1)
blavFitIndices(res1)

'
Posterior mean (EAP) of devm-based fit indices:

      BRMSEA    BGammaHat adjBGammaHat          BMc 
       0.017        0.998        0.994        0.996 
'


res1@external$samplls[,,2]



BayesRelFit(obs=res1@external$samplls[1,,1], rep=res1@external$samplls[1,,2], nvar=6, 
            pD=1358.2, N=250)


res1@external$p_waic

res1@ParTable
# semTools::fitmeasures(res1, c("dic", "pd"))
# semTools::fitmeasures(res1)


BayesFIT (fit1=res1)
loo::loo(res1)
library(coda)
mcmc_samples <- extract(res1)
as.mcmc.list(mcmc_samples)
loo(res1)


chisqs <- as.numeric(apply(res1@external$samplls, 2,
                           function(x) 2*(x[,2] - x[,1])))

chisqs






blavFitIndices(res1)





blavFitIndices1(res1)





blavInspect(res1, 'ntotal') 





################### BY HANDS #####################
pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"


chisqs <- as.numeric(apply(res1@external$samplls, 2,
                           function(x) 2*(x[,2] - x[,1])))    # $$$$$$$$$$$$$$$$$$$$$$ NEED
fit_pd <- fitMeasures(res1, paste0('p_', pD))                 # $$$$$$$$$$$$$$$$$$$$$$ NEED

fitMeasures(res1)


BayesChiFit(obs = chisqs, 
            nvar = 6, pD = fit_pd[1],
            N = 250,
            fit.measures = fit.measures, ms = FALSE, null_model = FALSE)#@details
loo(model_real_stan)$p_loo[1]

my_array1 <- Create2DimMatrix(posterior_samples$log_lik, posterior_samples$log_lik_sat, 1000)
my_array1[,,2]
chisqs1 <- as.numeric(apply(my_array1, 2,
                           function(x) 2*(x[,2] - x[,1]))) 
BayesChiFit(obs = chisqs1, 
            nvar = 6, pD = loo(model_real_stan)$p_loo[1],
            N = 250,
            fit.measures = fit.measures, ms = FALSE, null_model = FALSE)


plot(chisqs, chisqs1)

#######################################################################

## Public function
blavFitIndices1 <- function(object, thin = 1, pD = c("loo","waic","dic"),
                           rescale = c("devM","ppmc"),
                           fit.measures = "all", baseline.model = NULL) {
  
    out <- BayesChiFit(obs = chisqs, reps = reps,
                       nvar = object@Model@nvar, pD = fit_pd,
                       N = blavInspect(object, 'ntotal'),
                       Ngr = blavInspect(object, 'ngroups'),
                       Min1 = blavInspect(object, 'options')$mimic == "EQS",
                       ms = blavInspect(object, 'meanstructure'),
                       rescale = rescale, fit.measures = fit.measures,
                       null_model = null_model, obs_null = chisq_null,
                       reps_null = reps_null, pD_null = pD_null)
  
  nChains <- blavInspect(object, 'n.chains')
  out@details <- c(out@details, list(n.chains = nChains))
  out
}





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




