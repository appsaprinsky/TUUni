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
  //corr_matrix[2] COV_MATRIX[Tt];
  //corr_matrix[2] COV_MATRIX_1[Tt];
  //corr_matrix[2] COV_MATRIX_2[Tt];
  // cholesky_factor_corr[2] COV_MATRIX[Tt];
  matrix<lower=0>[6, Tt] ly;
  real<lower=0.1, upper=0.9> target_correlation;
  //matrix<lower=0>[2, Tt] error_1;
  //matrix<lower=0>[N, 2] error_1[Tt];
  //corr_matrix[2] COV_MATRIX_TEMP[Tt];

  matrix<lower=0>[6, Tt] MULTI_epsilon[Tt]; 
  corr_matrix[2] COV_MATRIX;

}


model {
  vector[2] temp_v1;
  vector[2] temp_v2;
  vector[2] temp_v;
  matrix[N, 2] eta_TEMP[Tt];
  //corr_matrix[2] COV_MATRIX_TEMP[Tt];
  matrix[2,2] COV_MATRIX_1[Tt];
  matrix[2,2] COV_MATRIX_2[Tt];
  matrix[6, Tt] MULTI_mean;
  matrix[6, Tt] MULTI_X;




  for (i in 1:N) {

    for (ttt in 1:Tt){
      MULTI_mean[1,ttt] = ly[1, ttt] * eta[i, 1]*(target_correlation^(ttt-1)) ;
      MULTI_mean[2,ttt] = ly[2, ttt] * eta[i, 1]*(target_correlation^(ttt-1)) ;
      MULTI_mean[3,ttt] = ly[3, ttt] * eta[i, 1]*(target_correlation^(ttt-1)) ;
      MULTI_mean[4,ttt] = ly[4, ttt] * eta[i, 2]*(target_correlation^(ttt-1)) ;
      MULTI_mean[5,ttt] = ly[5, ttt] * eta[i, 2]*(target_correlation^(ttt-1)) ;
      MULTI_mean[6,ttt] = ly[6, ttt] * eta[i, 2]*(target_correlation^(ttt-1)) ;

      MULTI_X[1,ttt] = x[i,1,Tt];
      MULTI_X[2,ttt] = x[i,2,Tt];
      MULTI_X[3,ttt] = x[i,3,Tt];
      MULTI_X[4,ttt] = x[i,4,Tt];
      MULTI_X[5,ttt] = x[i,5,Tt];
      MULTI_X[6,ttt] = x[i,6,Tt];

      //MULTI_X[1:6, ttt] ~ multi_normal(MULTI_mean[1:6, ttt], MULTI_epsilon[ttt]);


      
    }

    //x[i,1:6,1:Tt] ~ multi_normal(MULTI_mean, MULTI_epsilon[1:6, 1:6]);
    //print(MULTI_X );
    //print(MULTI_mean);
    //print(MULTI_epsilon[1:6, 1:6]);
    //MULTI_X ~ multi_normal(MULTI_mean, MULTI_epsilon[1:6, 1:6]);
    eta[i,1:2] ~ multi_normal(mu_0, COV_MATRIX); 


    //x[i,1,1] ~ normal(ly[1, 1]*eta[i, 1], epsilon[1, 1]);
    //x[i,2,1] ~ normal(ly[2, 1]*eta[i, 1], epsilon[2, 1]);
    //x[i,3,1] ~ normal(ly[3, 1]*eta[i, 1], epsilon[3, 1]);
    //x[i,4,1] ~ normal(ly[4, 1]*eta[i, 2], epsilon[4, 1]);
    //x[i,5,1] ~ normal(ly[5, 1]*eta[i, 2], epsilon[5, 1]);
    //x[i,6,1] ~ normal(ly[6, 1]*eta[i, 2], epsilon[6, 1]);
    //eta_TEMP[1][i, 1] = eta[i, 1];
    //eta_TEMP[1][i, 1] = eta[i, 1];
    //eta_TEMP[1][i, 1] = eta[i, 1];
    //eta_TEMP[1][i, 2] = eta[i, 2];
    //eta_TEMP[1][i, 2] = eta[i, 2];
    //eta_TEMP[1][i, 2] = eta[i, 2];
    //eta[i,1:2] ~ multi_normal(mu_0, COV_MATRIX[1]); 

    //if (Tt>1){
      //for (t in 2:Tt){ 
        //eta_TEMP[t][i,1] = eta_TEMP[t-1][i,1] * 0.6;//target_correlation; //+ error_1[t][i,1];
        //eta_TEMP[t][i,2] = eta_TEMP[t-1][i,2] * 0.6;//target_correlation; //+ error_1[t][i,2];
        //x[i,1,t] ~ normal(ly[1, t]*eta_TEMP[t][i, 1], epsilon[1, t]);
        //x[i,2,t] ~ normal(ly[2, t]*eta_TEMP[t][i, 1], epsilon[2, t]);
        //x[i,3,t] ~ normal(ly[3, t]*eta_TEMP[t][i, 1], epsilon[3, t]);
        //x[i,4,t] ~ normal(ly[4, t]*eta_TEMP[t][i, 2], epsilon[4, t]);
        //x[i,5,t] ~ normal(ly[5, t]*eta_TEMP[t][i, 2], epsilon[5, t]);
        //x[i,6,t] ~ normal(ly[6, t]*eta_TEMP[t][i, 2], epsilon[6, t]);
        //eta[i,1:2] ~ multi_normal(mu_0, COV_MATRIX[t]); 
        //temp_v1[1] = eta_TEMP[t][i,1];
        //temp_v1[2] = eta_TEMP[t-1][i,1];
        //temp_v2[1] = eta_TEMP[t][i,2];
        //temp_v2[2] = eta_TEMP[t-1][i,2];  
        //COV_MATRIX_1[Tt][2,1] = target_correlation;
        //COV_MATRIX_1[Tt][1,2] = target_correlation;
        //COV_MATRIX_1[Tt][2,2] = COV_MATRIX[1][2,2];
        //COV_MATRIX_1[Tt][1,1] = COV_MATRIX[1][1,1];
        //COV_MATRIX_2[Tt][2,1] = target_correlation;
        //COV_MATRIX_2[Tt][1,2] = target_correlation; 
        //COV_MATRIX_2[Tt][2,2] = COV_MATRIX[1][2,2];
        //COV_MATRIX_2[Tt][1,1] = COV_MATRIX[1][1,1]; 
        //temp_v1 ~ multi_normal(mu_0, COV_MATRIX_1[Tt]);  
        //temp_v2 ~ multi_normal(mu_0, COV_MATRIX_2[Tt]);  
        //eta_TEMP[t][i,1:2] ~ multi_normal(mu_0, COV_MATRIX[1]);   

      //}
    //}
  }



  //for (t in 1:Tt) {
    //for (i in 1:N) {
      //x[i,1,t] ~ normal(ly[1, t]*eta_TEMP[t][i, 1], epsilon[1, t]);
      //x[i,2,t] ~ normal(ly[2, t]*eta_TEMP[t][i, 1], epsilon[2, t]);
      //x[i,3,t] ~ normal(ly[3, t]*eta_TEMP[t][i, 1], epsilon[3, t]);
      //x[i,4,t] ~ normal(ly[4, t]*eta_TEMP[t][i, 2], epsilon[4, t]);
      //x[i,5,t] ~ normal(ly[5, t]*eta_TEMP[t][i, 2], epsilon[5, t]);
      //x[i,6,t] ~ normal(ly[6, t]*eta_TEMP[t][i, 2], epsilon[6, t]);
      //eta_TEMP[t][i,1:2] ~ multi_normal(mu_0, COV_MATRIX[t]); 
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
  //}
  
  //priors
  //ly ~ normal(0,10);
  //epsilon ~ cauchy(0,5);
  //COV_MATRIX ~ lkj_corr(1);
}


