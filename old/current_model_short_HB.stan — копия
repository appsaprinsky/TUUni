data {
  int<lower=0> N;
  matrix[N, 6] x; 
  vector[2] mu_0;
  vector[6] x_mean;
  matrix[6, 6] cov_x_DATA;
  matrix[6, 6] identity_matrix;
  vector[6] zero_vector;
}

parameters {
  matrix[N, 2]  eta;
  vector<lower=0>[6] epsilon; 
  corr_matrix[2] COV_MATRIX;
  vector<lower=0>[6] ly;
}

model {
  vector[6] TEMP_MEAN;
  matrix[6, 6] TEPM_VAR;    

  for (i in 1:N) {
    x[i,1] ~ normal(ly[1]*eta[i, 1], epsilon[1]);
    x[i,2] ~ normal(ly[2]*eta[i, 1], epsilon[2]);
    x[i,3] ~ normal(ly[3]*eta[i, 1], epsilon[3]);
    x[i,4] ~ normal(ly[4]*eta[i, 2], epsilon[4]);
    x[i,5] ~ normal(ly[5]*eta[i, 2], epsilon[5]);
    x[i,6] ~ normal(ly[6]*eta[i, 2], epsilon[6]);

    eta[i,1:2] ~ multi_normal(mu_0,COV_MATRIX);    
    //note: use cholesky for computational reasons later
  }
  
  //priors
  ly ~ normal(0,10);
  epsilon ~ cauchy(0,5);
  COV_MATRIX ~ lkj_corr(1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] log_lik_sat;
  matrix[6, 6] curr_cov_ll;
  vector[6] x_together_ml;
  vector[6] mu_together_ml;  
  real sigma;
  matrix[6,6] epsm;
  matrix[6,2] lym;

  sigma = COV_MATRIX[2,1];  
  
  // this is the mean across persons
  // below is the person-specific mean (with same name)
  for(j in 1:6){mu_together_ml[j] = 0;}  

  for (i in 1:N){
    for(j in 1:3){mu_together_ml[j] += ly[j]*eta[i, 1];}
    for(j in 4:6){mu_together_ml[j] += ly[j]*eta[i, 2];}
  }
  
  mu_together_ml = mu_together_ml/N;
  
  // this is the model-implied covariance matrix
    // general form is cov=ly%*%Sigma%*%t(ly)+eps_M
    // with ly vector
    
    for(j in 1:6){
      epsm[j,j] = epsilon[j]*epsilon[j];
      for(k in 1:6){
        if(k!=j){epsm[k,j] = 0;}
      }
    }
    
    for(j in 1:3){
      lym[j,1] =ly[j];
      lym[j,2] =0;
    }  
    for(j in 4:6){
      lym[j,2] =ly[j];
      lym[j,1] =0;
    }  

    curr_cov_ll = lym * COV_MATRIX * lym'+ epsm;
    
    ////////////
    for (i in 1:N) {  
    for (j in 1:6){
      x_together_ml[j] = x[i,j]; 
    }
    
    // log_like
    log_lik[i] = 0;  
    log_lik[i] = multi_normal_lpdf(x_together_ml | mu_together_ml , curr_cov_ll); // (mu_together_ml, curr_cov_ll) mu_zeros_ml
    
    log_lik_sat[i] = 0;
    log_lik_sat[i] = multi_normal_lpdf(x_together_ml | x_mean, cov_x_DATA);  // 
    ////////////
    
    
    
  }
}
