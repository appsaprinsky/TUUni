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
res1 <- bsem(model, data=as.data.frame(x),#target="jags",    
             burnin = 500, n.chains = 1, sample = 1000, 
             mcmcfile = "model1", std.lv = TRUE) 
# extract model fit
res1 # summary(res1)
summary(res1)
blavFitIndices(res1)
res1@external$samplls[,,2]
res1@external$samplls[,,1]

stan_model_code <- '
data {
  int<lower=0> N;
  matrix[N, 6] x; 
  // vector[4] ly;
}

parameters {
  matrix[N, 2]  eta;
  vector<lower=0>[6] epsilon; 
  cov_matrix[2] sigma; 
  vector[2] mu_0; 
  vector[4] ly;
}

transformed parameters {
  matrix[N, 6] mu_x;

  for (i in 1:N) {
    mu_x[i, 1] = eta[i, 1];
    mu_x[i, 2] = ly[1]*eta[i, 1];
    mu_x[i, 3] = ly[2]*eta[i, 1];

    mu_x[i, 4] = eta[i, 2];
    mu_x[i, 5] = ly[3]*eta[i, 2];
    mu_x[i, 6] = ly[4]*eta[i, 2];
  }
}

model {

  for(k in 1:4){
    ly[k] ~ normal(0, 1);  // DELETE LATERRRRRRR
  }

  for(k in 1:6){
    epsilon[k] ~ cauchy(1, 1);   // DELETE LATERRRRRRR gamma
  }
  
  for (i in 1:N) {
    for (j in 1:6) {
      x[i,j] ~ normal(mu_x[i,j], epsilon[j]);
    }
    eta[i, 1:2] ~ multi_normal(mu_0,sigma);
  }
  // epsilon ~ cauchy(1, 1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] log_lik_sat;
  for (i in 1:N) {
    log_lik[i] = 0;  // normal_lpdf(x[i,1] | mu_x[i,1], epsilon[1]);
    for (j in 1:6) {
      log_lik[i] += normal_lpdf(x[i,j] | mu_x[i,j], epsilon[j]);
    }
    
    log_lik_sat[i] = 0;
    for (j in 1:6) {
      log_lik_sat[i] += normal_lpdf(x[i,j] | mu_x[i,j], 1);
    }
    
  }
}
'
#######################################STAN CODE###########################

stan_model_code <- '
data {
  int<lower=0> N;
  matrix[N, 6] x; 
  // vector[4] ly;
}

parameters {
  matrix[N, 2]  eta;
  vector<lower=0>[6] epsilon; 
  cov_matrix[2] sigma; 
  vector[2] mu_0; 
  vector[4] ly;
}

transformed parameters {
  matrix[N, 6] mu_x;

  for (i in 1:N) {
    mu_x[i, 1] = eta[i, 1];
    mu_x[i, 2] = ly[1]*eta[i, 1];
    mu_x[i, 3] = ly[2]*eta[i, 1];

    mu_x[i, 4] = eta[i, 2];
    mu_x[i, 5] = ly[3]*eta[i, 2];
    mu_x[i, 6] = ly[4]*eta[i, 2];
  }
}

model {

  for(k in 1:4){
    ly[k] ~ normal(0, 1);  // DELETE LATERRRRRRR
  }

  for(k in 1:6){
    epsilon[k] ~ cauchy(1, 1);   // DELETE LATERRRRRRR gamma  0, 1
  }
  
  for (i in 1:N) {
    for (j in 1:6) {
      x[i,j] ~ normal(mu_x[i,j], epsilon[j]);
    }
    eta[i, 1:2] ~ multi_normal(mu_0,sigma);
  }
  // epsilon ~ cauchy(1, 1);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = 0;  // normal_lpdf(x[i,1] | mu_x[i,1], epsilon[1]);
    for (j in 1:6) {
      log_lik[i] += normal_lpdf(x[i,j] | mu_x[i,j], epsilon[j]);
    }
  }
}
'




stan_model_code <- '
data {
  int<lower=0> N;
  matrix[N, 6] x; 
  vector[2] mu_0;
  // vector[4] ly;
}

parameters {
  matrix[N, 2]  eta;
  vector<lower=0>[6] epsilon; 
  corr_matrix[2] sigma; 
  // vector[2] mu_0;
  vector[6] ly;
}

model {

  for(k in 1:6){
    ly[k] ~ normal(0, 1);  // DELETE LATERRRRRRR
  }
  
  sigma ~ lkj_corr(1);

  for(k in 1:6){
    epsilon[k] ~ cauchy(0, 1);   // DELETE LATERRRRRRR gamma  0, 1
  }
  
  for (i in 1:N) {
  
    x[i,1] ~ normal(ly[5]*eta[i, 1], epsilon[1]);
    x[i,2] ~ normal(ly[1]*eta[i, 1], epsilon[2]);
    x[i,3] ~ normal(ly[2]*eta[i, 1], epsilon[3]);
    x[i,4] ~ normal(ly[6]*eta[i, 2], epsilon[4]);
    x[i,5] ~ normal(ly[3]*eta[i, 2], epsilon[5]);
    x[i,6] ~ normal(ly[4]*eta[i, 2], epsilon[6]);
    
    eta[i, 1:2] ~ multi_normal(mu_0,sigma);
  }
  // epsilon ~ cauchy(1, 1);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = 0;  
    log_lik[i] += normal_lpdf(x[i,1] | ly[5]*eta[i, 1], epsilon[1]);
    log_lik[i] += normal_lpdf(x[i,2] | ly[1]*eta[i, 1], epsilon[2]);
    log_lik[i] += normal_lpdf(x[i,3] | ly[2]*eta[i, 1], epsilon[3]);
    log_lik[i] += normal_lpdf(x[i,4] | ly[6]*eta[i, 2], epsilon[4]);
    log_lik[i] += normal_lpdf(x[i,5] | ly[3]*eta[i, 2], epsilon[5]);
    log_lik[i] += normal_lpdf(x[i,6] | ly[4]*eta[i, 2], epsilon[6]);
  }
}
'



stan_data <- list(
  N=N,
  x=x
)


stan_data <- list(
  N=N,
  x=x,
  mu_0=c(0,0)
)


model_real_stan <- stan(model_code = stan_model_code, iter = 1000, verbose = FALSE, data=stan_data, chains=1, warmup = 500 )
posterior_samples <- extract(model_real_stan)
posterior_samples$log_lik
posterior_samples$log_lik_sat

# posterior_samples$target

ddd <- data.frame(posterior_samples$log_lik)
sum(ddd[1,])
rowSums (ddd)
plot(rowSums(ddd), res1@external$samplls[501:1000,,1])

ddd <- ddd %>%
  mutate(ll = rowSums(., na.rm=TRUE))

ddd$ll
res1@external$samplls[1:500,,1]
res1@external$samplls[,,2]



sum(data.frame(posterior_samples$log_lik_sat)[1,])
sum(data.frame(posterior_samples$log_lik)[1,])




res <- rstan::read_stan_csv(model_real_stan $output_files())



loo::extract_log_lik(model_real_stan, parameter_name = "log_lik_sat")


llsat <- loo::extract_log_lik(lavjags, parameter_name = "log_lik_sat")








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


