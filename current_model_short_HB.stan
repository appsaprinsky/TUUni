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
  //matrix[N, 2] eta;
  matrix<lower=0>[6, Tt] epsilon; 
  // real COV_MATRIX[2, 2, Tt];
  corr_matrix[2] COV_MATRIX[Tt];
  //corr_matrix[2] COV_MATRIX_1[Tt];
  //corr_matrix[2] COV_MATRIX_2[Tt];
  // cholesky_factor_corr[2] COV_MATRIX[Tt];
  matrix<lower=0>[6, Tt] ly;
  real<lower=0.1, upper=0.9> target_correlation;
  //matrix<lower=0>[2, Tt] error_1;
  matrix<lower=0>[N, 2] error_1[Tt];
  //corr_matrix[2] COV_MATRIX_TEMP[Tt];



  vector[2]      beta; //ar(1) coefficients
  corr_matrix[2] phi;
  matrix[N,2] eta[Tt];
}


model {
  vector[2] temp_v1;
  vector[2] temp_v2;
  vector[2] temp_v;
  matrix[N, 2] eta_TEMP[Tt];
  //corr_matrix[2] COV_MATRIX_TEMP[Tt];
  matrix[2,2] COV_MATRIX_1[Tt];
  matrix[2,2] COV_MATRIX_2[Tt];


  matrix[N,2] mueta[Tt];
  for(i in 1:N){
    //t==1
    mueta[1,i,1] = 0;
    mueta[1,i,2] = 0;
  
    //t==rest
    for(t in 2:Tt){
      mueta[t,i,1] = beta[1]*eta[t-1,i,1];
      mueta[t,i,2] = beta[2]*eta[t-1,i,2];
    }
  
  }


  for (i in 1:N) {
    for(t in 1:Tt){
      x[i,1,t] ~ normal(ly[1, 1]*eta[t,i,1], epsilon[1, 1]);
      x[i,2,t] ~ normal(ly[2, 1]*eta[t,i,1], epsilon[2, 1]);
      x[i,3,t] ~ normal(ly[3, 1]*eta[t,i,1], epsilon[3, 1]);
      x[i,4,t] ~ normal(ly[4, 1]*eta[t,i,2], epsilon[4, 1]);
      x[i,5,t] ~ normal(ly[5, 1]*eta[t,i,2], epsilon[5, 1]);
      x[i,6,t] ~ normal(ly[6, 1]*eta[t,i,2], epsilon[6, 1]);
        
      eta[t,i,1:2] ~ multi_normal(mueta[t,i,1:2],phi);    

    }
    
  }

  beta ~ normal(0,1);
  phi ~ lkj_corr(1);

}

