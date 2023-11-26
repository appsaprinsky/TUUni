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
  //matrix[N, 2] eta[Tt];
  matrix[N, 2] eta;
  matrix<lower=0>[6, Tt] epsilon; 
  // real COV_MATRIX[2, 2, Tt];
  corr_matrix[2] COV_MATRIX[Tt];
  corr_matrix[2] COV_MATRIX_1[Tt];
  corr_matrix[2] COV_MATRIX_2[Tt];
  // cholesky_factor_corr[2] COV_MATRIX[Tt];
  matrix<lower=0>[6, Tt] ly;
  real<lower=0.5, upper=0.7> target_correlation;
  //matrix<lower=0>[2, Tt] error_1;
  matrix<lower=0>[N, 2] error_1[Tt];
  //corr_matrix[2] COV_MATRIX_TEMP[Tt];
}


model {
  vector[2] temp_v1;
  vector[2] temp_v2;
  vector[2] temp_v;
  matrix[N, 2] eta_TEMP[Tt];
  //corr_matrix[2] COV_MATRIX_TEMP[Tt];

  for (t in 1:Tt) {

    if (t > 1){
      //eta_TEMP[t][,1] = eta_TEMP[t-1][,1] * target_correlation + error_1[t][,1];
      //eta_TEMP[t][,2] = eta_TEMP[t-1][,2] * target_correlation + error_1[t][,2];

      for (i in 1:N) {
        vector[2] temp_vector = multi_normal_rng(mu_0, COV_MATRIX[t-1] );
        eta_TEMP[t][i,] = temp_vector;
      }


      eta_TEMP[t][,1] = target_correlation * eta_TEMP[t-1][,1] + sqrt(1 - target_correlation^2) * eta_TEMP[t][,1];
      eta_TEMP[t][,2] = target_correlation * eta_TEMP[t-1][,2] + sqrt(1 - target_correlation^2) * eta_TEMP[t][,2];

      //print(eta_TEMP[t][,1]);
      for (i in 1:N) {
        x[i,1,t] ~ normal(ly[1, t]*eta_TEMP[t][i, 1], epsilon[1, t]);
        x[i,2,t] ~ normal(ly[2, t]*eta_TEMP[t][i, 1], epsilon[2, t]);
        x[i,3,t] ~ normal(ly[3, t]*eta_TEMP[t][i, 1], epsilon[3, t]);
        x[i,4,t] ~ normal(ly[4, t]*eta_TEMP[t][i, 2], epsilon[4, t]);
        x[i,5,t] ~ normal(ly[5, t]*eta_TEMP[t][i, 2], epsilon[5, t]);
        x[i,6,t] ~ normal(ly[6, t]*eta_TEMP[t][i, 2], epsilon[6, t]);
        eta_TEMP[t][i,1:2] ~ multi_normal(mu_0, COV_MATRIX[t]); 
      }       
    } else {
      eta_TEMP[t][,1] = eta[,1];
      eta_TEMP[t][,2] = eta[,2];
      for (i in 1:N) {
        x[i,1,t] ~ normal(ly[1, t]*eta[i, 1], epsilon[1, t]);
        x[i,2,t] ~ normal(ly[2, t]*eta[i, 1], epsilon[2, t]);
        x[i,3,t] ~ normal(ly[3, t]*eta[i, 1], epsilon[3, t]);
        x[i,4,t] ~ normal(ly[4, t]*eta[i, 2], epsilon[4, t]);
        x[i,5,t] ~ normal(ly[5, t]*eta[i, 2], epsilon[5, t]);
        x[i,6,t] ~ normal(ly[6, t]*eta[i, 2], epsilon[6, t]);
        eta[i,1:2] ~ multi_normal(mu_0, COV_MATRIX[t]);
      } 
    }

    //for (i in 1:N) {
    //  x[i,1,t] ~ normal(ly[1, t]*eta_TEMP[t][i, 1], epsilon[1, t]);
    //  x[i,2,t] ~ normal(ly[2, t]*eta_TEMP[t][i, 1], epsilon[2, t]);
    //  x[i,3,t] ~ normal(ly[3, t]*eta_TEMP[t][i, 1], epsilon[3, t]);
    //  x[i,4,t] ~ normal(ly[4, t]*eta_TEMP[t][i, 2], epsilon[4, t]);
    //  x[i,5,t] ~ normal(ly[5, t]*eta_TEMP[t][i, 2], epsilon[5, t]);
    //  x[i,6,t] ~ normal(ly[6, t]*eta_TEMP[t][i, 2], epsilon[6, t]);
    //  eta_TEMP[t][i,1:2] ~ multi_normal(mu_0, COV_MATRIX[t]); 

      //eta[t][i,1:2] ~ multi_normal(mu_0, multiply_lower_tri_self_transpose(COV_MATRIX[t]));   
      //temp_v[1] = eta[t][i,1];
      //temp_v[2] = eta[t][i,2];
      //if (t > 1){
        //eta[t][i,1] = eta[t-1][i,1] * target_correlation + error_1[1,t];
        //eta[t][i,2] = eta[t-1][i,2] * target_correlation + error_1[2,t];


        //temp_v[1] = target_correlation * eta[t-1][i,1] + sqrt(1-target_correlation^2)*eta[t][i,1];
        //temp_v[2] = target_correlation * eta[t-1][i,2] + sqrt(1-target_correlation^2)*eta[t][i,2];
        //temp_v1[1] = target_correlation * eta[t-1][i,1] + sqrt(1-target_correlation^2)*eta[t][i,1]; //eta[t-1][i,1];
        //temp_v1[2] = eta[t][i,1];
        //temp_v1 ~ multi_normal(mu_0, COV_MATRIX_1[Tt]);
        //temp_v2[1] = target_correlation * eta[t-1][i,2] + sqrt(1-target_correlation^2)*eta[t][i,2]; //eta[t-1][i,2];
        //temp_v2[2] = eta[t][i,2];
        //temp_v2 ~ multi_normal(mu_0, COV_MATRIX_2[Tt]);
        //eta[t-1:t][i,1] ~ multi_normal(mu_0, COV_MATRIX_1[t]);
        //eta[t-1:t][i,2] ~ multi_normal(mu_0, COV_MATRIX_2[t]); 

      //}
      //note: use cholesky for computational reasons later
      //temp_v ~ multi_normal(mu_0, COV_MATRIX[t]);
    //}
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
  matrix[N, 2] eta_TEMP[Tt];
  
  // this is the mean across persons
  // below is the person-specific mean (with same name)
  for (t in 1:Tt) {   
  for(j in 1:6){mu_together_ml[j, t] = 0;} 

  if (t > 1){
    eta_TEMP[t][,1] = eta_TEMP[t-1][,1] * target_correlation + error_1[t][,1];
    eta_TEMP[t][,2] = eta_TEMP[t-1][,2] * target_correlation + error_1[t][,2];
  } else {
    eta_TEMP[t][,1] = eta[,1];
    eta_TEMP[t][,2] = eta[,2];
  }

  for (i in 1:N){
    for(j in 1:3){mu_together_ml[j, t] += ly[j, t]*eta_TEMP[t][i, 1];}
    for(j in 4:6){mu_together_ml[j, t] += ly[j, t]*eta_TEMP[t][i, 2];}
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
