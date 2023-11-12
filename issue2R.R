###########################################
# Set Up for Issue
###########################################
library(loo)
library(dplyr)
library(blavaan) # standard package
library(lavaan) # standard package
library(rstan)
library(mvtnorm) # mv normal data
packageVersion("blavaan")
packageVersion("rjags")
packageVersion("rstan")
set.seed(1024)
N <- 250
phi <- diag(2)
eta <- rmvnorm(N,sigma=phi)
epsilon <- rmvnorm(N,sigma=diag(6))
lambda <- matrix(c(1,1,1,0,0,0,
                   0,0,0,1,1,1),6,2,byrow=F) 
x <- eta%*%t(lambda)+epsilon 

model <- '
eta1 =~ V1+V2+V3
eta2 =~ V4+V5+V6
'

res1 <- bsem(model, data=as.data.frame(x),#target="jags",    
             burnin = 500, n.chains = 1, sample = 1000, 
             mcmcfile = "model1", std.lv = TRUE) 
res1 
summary(res1)

LY_USED <- as.vector(res1@external$stansumm[1:6,1])
VARIANCE_USED <- as.vector(res1@external$stansumm[7:12,1])
Psi_cov_USED <- as.vector(res1@external$stansumm[13,1])
LOG_LIKE_SAT <- as.vector(res1@external$samplls[,,2])[1000]
LOG_LIKE <- as.vector(res1@external$samplls[,,1])[1000]


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
  real sigma_USE;
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
  
    x[i,1] ~ normal(ly[1]*eta[i, 1], epsilon[1]);
    x[i,2] ~ normal(ly[2]*eta[i, 1], epsilon[2]);
    x[i,3] ~ normal(ly[3]*eta[i, 1], epsilon[3]);
    x[i,4] ~ normal(ly[4]*eta[i, 2], epsilon[4]);
    x[i,5] ~ normal(ly[5]*eta[i, 2], epsilon[5]);
    x[i,6] ~ normal(ly[6]*eta[i, 2], epsilon[6]);


//    TEMP_MEAN[1] = ly[1]*eta[i, 1];
//    TEMP_MEAN[2] = ly[2]*eta[i, 1];    
//    TEMP_MEAN[3] = ly[3]*eta[i, 1];    
//    TEMP_MEAN[4] = ly[4]*eta[i, 2];    
//    TEMP_MEAN[5] = ly[5]*eta[i, 2];    
//    TEMP_MEAN[6] = ly[6]*eta[i, 2];
    

//    for (aaa in 1:6) {
//      for (sss in 1:6) {
//        if (aaa == sss) {
//          TEPM_VAR[aaa, sss] = ly[aaa]*ly[aaa] + epsilon[aaa]; 
//        } else {
//          if ((aaa <= 3 && sss <= 3) || (aaa >= 4 && sss >= 4)){
//            TEPM_VAR[aaa, sss] = ly[aaa]*ly[sss];  
//          } else {
//            TEPM_VAR[aaa, sss] = ly[aaa]*ly[sss]*sigma; 
//          }
//        }
//      }
//    }       
    
    
    
//    x[i,1:6] ~ multi_normal(TEMP_MEAN,TEPM_VAR);    
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
          curr_cov_ll[aaa, sss] = LY_TRIAL[aaa]*LY_TRIAL[aaa] + VARIANCE_TRIAL[aaa]; //epsilon[aaa];  VARIANCE_TRIAL[aaa]; TRIAL
        } else {
          if ((aaa <= 3 && sss <= 3) || (aaa >= 4 && sss >= 4)){
            curr_cov_ll[aaa, sss] = LY_TRIAL[aaa]*LY_TRIAL[sss]; // curr_cov_ll[aaa, sss] = cov_x_DATA[aaa, sss]; 
          } else {
            curr_cov_ll[aaa, sss] = LY_TRIAL[aaa]*LY_TRIAL[sss]*sigma_USE; // curr_cov_ll[aaa, sss] = cov_x_DATA[aaa, sss]; 
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
  cov_x_DATA=cov(x),
  identity_matrix = diag(6),
  LY_TRIAL = LY_USED,
  VARIANCE_TRIAL = VARIANCE_USED,
  zero_vector = c(0,0,0,0,0,0),
  sigma_USE = Psi_cov_USED
)
model_real_stan <- stan(model_code = stan_model_code, iter = 1500, verbose = FALSE, data=stan_data, chains=1, warmup = 500 )
posterior_samples <- extract(model_real_stan)


ddd <- data.frame(posterior_samples$log_lik)
sum(ddd[1,])
rowSums (ddd)
ddd1 <- data.frame(posterior_samples$log_lik_sat)
sum(ddd1[1,])
rowSums (ddd1)

print('Compare Log Saturated Models: ')
sum(ddd1[1,]) #=== MINE
LOG_LIKE_SAT  # BLAVAAN
# -2492.285 vs -2492.279

print('Compare Log LIK Models: ')
sum(ddd[1,]) #=== MINE
LOG_LIKE  # BLAVAAN
# -2497.986 vs -2501.028

## There is still a difference!!!

'
Is this still appropriate or am I missing Something?
Is My Log Like Is still wrong?

For a full procedure without pre computed values look at 

CurrentVersion.R
'




