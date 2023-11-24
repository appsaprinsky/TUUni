data {
  int<lower=0> Tt;
  int<lower=0> N;
  real x[N, 6, Tt]; 
  //  matrix[N,6] x[Tt];
  // matrix[2, Tt] mu_0;
  vector[2] mu_0;
  matrix[6, Tt] x_mean;
  // real cov_x_DATA[6, 6, Tt];
  matrix[6,6] cov_x_DATA[Tt];  
}

parameters {
  // real eta[N, 2, Tt];
  matrix[N, 2] eta[Tt];
  matrix<lower=0>[6, Tt] epsilon; 
  // real COV_MATRIX[2, 2, Tt];
  corr_matrix[2] COV_MATRIX[Tt];
  // cholesky_factor_corr[2] COV_MATRIX[Tt];
  matrix<lower=0>[6, Tt] ly;
}


model {
  for (t in 1:Tt) {
    for (i in 1:N) {
      x[i,1,t] ~ normal(ly[1, t]*eta[t][i, 1], epsilon[1, t]);
      x[i,2,t] ~ normal(ly[2, t]*eta[t][i, 1], epsilon[2, t]);
      x[i,3,t] ~ normal(ly[3, t]*eta[t][i, 1], epsilon[3, t]);
      x[i,4,t] ~ normal(ly[4, t]*eta[t][i, 2], epsilon[4, t]);
      x[i,5,t] ~ normal(ly[5, t]*eta[t][i, 2], epsilon[5, t]);
      x[i,6,t] ~ normal(ly[6, t]*eta[t][i, 2], epsilon[6, t]);
      //eta[t][i,1:2] ~ multi_normal(mu_0, multiply_lower_tri_self_transpose(COV_MATRIX[t]));   
      eta[t][i,1:2] ~ multi_normal(mu_0, COV_MATRIX[t]); 
      //note: use cholesky for computational reasons later
    }
  }
  
  //priors
  //ly ~ normal(0,10);
  //epsilon ~ cauchy(0,5);
  //COV_MATRIX ~ lkj_corr(1);
}

generated quantities {
  matrix[N,Tt] log_lik;
  matrix[N,Tt] log_lik_sat;
  //real curr_cov_ll[6, 6, Tt];
  matrix[6, 6] curr_cov_ll[Tt];
  matrix[6,Tt] x_together_ml;
  matrix[6,Tt] mu_together_ml;  
  // real epsm[6,6,Tt];
  matrix[6, 6] epsm[Tt];
  // real lym[6,2,Tt];
  matrix[6, 2] lym[Tt];
  
  // this is the mean across persons
  // below is the person-specific mean (with same name)
  for (t in 1:Tt) {   
  for(j in 1:6){mu_together_ml[j, t] = 0;}  

  for (i in 1:N){
    for(j in 1:3){mu_together_ml[j, t] += ly[j, t]*eta[t][i, 1];}
    for(j in 4:6){mu_together_ml[j, t] += ly[j, t]*eta[t][i, 2];}
  }
  
  mu_together_ml[, t] = mu_together_ml[, t]/N;
  
  // this is the model-implied covariance matrix
    // general form is cov=ly%*%Sigma%*%t(ly)+eps_M
    // with ly vector
    
    for(j in 1:6){
      epsm[t][j,j] = epsilon[j, t]*epsilon[j, t];
      for(k in 1:6){
        if(k!=j){epsm[t][k,j] = 0;}
      }
    }
    
    for(j in 1:3){
      lym[t][j,1] =ly[j,t];
      lym[t][j,2] =0;
    }  
    for(j in 4:6){
      lym[t][j,2] =ly[j,t];
      lym[t][j,1] =0;
    }  

    ////////////  print("dd", ddd);
    //curr_cov_ll[t] = lym[t] * multiply_lower_tri_self_transpose(COV_MATRIX[t]) * lym[t]'+ epsm[t];
    curr_cov_ll[t] = lym[t] * COV_MATRIX[t] * lym[t]'+ epsm[t];
    for (i in 1:N) {  
    for (j in 1:6){
      x_together_ml[j,t] = x[i,j,t]; 
    }
    
    // log_like
    log_lik[i,t] = 0;  
    log_lik[i,t] = multi_normal_lpdf(x_together_ml[, t] | mu_together_ml[, t] , curr_cov_ll[t]); 

    log_lik_sat[i,t] = 0;
    log_lik_sat[i,t] = multi_normal_lpdf(x_together_ml[, t] | x_mean[, t], cov_x_DATA[t]);
    ////////////
    
    
    
  }
}
}
