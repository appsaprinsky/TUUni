data {
  int<lower=0> N;
  int<lower=0> Nt;
  //real y[N,Nt,6]; 
  matrix[N,6] y[Nt];
  vector[Nt*6] x_mean;
  matrix[Nt*6,Nt*6] x_cov;
}

parameters {
  vector<lower=0>[6] ly;
  vector<lower=0>[6] sigmaeps; //resvar for obs var
  vector[2]      beta; //ar(1) coefficients
  corr_matrix[2] phi;
  
  //real  eta[N,Nt,2];
  matrix[N,2] eta[Nt]; // Eta over Tt data points
  
}

model {
  matrix[N,2] mueta[Nt];
  
    // mean structure for eta, simple AR(1), no crossloadings
  for(i in 1:N){
    //t==1
    mueta[1,i,1] = 0;
    mueta[1,i,2] = 0;
  
    //t==rest
    for(j in 2:Nt){
      mueta[j,i,1] = beta[1]*eta[j-1,i,1];
      mueta[j,i,2] = beta[2]*eta[j-1,i,2];
    }
  
  }


  for (i in 1:N) {
    for(j in 1:Nt){
      for(k in 1:3){y[j,i,k] ~ normal(ly[k]*eta[j,i,1], sigmaeps[k]);}
      for(k in 4:6){y[j,i,k] ~ normal(ly[k]*eta[j,i,2], sigmaeps[k]);}

      eta[j,i,1:2] ~ multi_normal(mueta[j,i,1:2],phi);    
    //note: use cholesky for computational reasons later
    //note: phi is time-invariant for now
    //note: no meanstructure for items for now
    }
    
  }

  //priors
  ly ~ normal(0,1);
  beta ~ normal(0,1);
  sigmaeps ~ cauchy(0,5);
  phi ~ lkj_corr(1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] log_lik_sat;
  
  matrix[6*Nt, 6*Nt] sigmayhat;//total cov
  
  
  matrix[N,6*Nt] x_together;
  vector[6*Nt]   mu_together= rep_vector(0, 6*Nt);  
  
  matrix[2*Nt, 2*Nt] sigmaeta; //cov over time
  
  // start with zero matrix
  matrix[2*Nt, 2*Nt] sigmazeta= rep_matrix(0, 2*Nt, 2*Nt); //residual variance
  matrix[2*Nt, 2*Nt] D = rep_matrix(0, 2*Nt, 2*Nt); //regression coefficient matrix
  matrix[6*Nt,6*Nt] epsm= rep_matrix(0, 6*Nt, 6*Nt);
  matrix[6*Nt,2*Nt] lym = rep_matrix(0, 6*Nt, 2*Nt);
  
  
  // this is the mean across persons
  // below is the person-specific mean (with same name)
  //for(j in 1:(6*Nt)){mu_together[j] = 0;}  
  // sequence is 1:6 for t=1, 1:6 for t=2 etc.
  
  for (i in 1:N){
    for(j in 1:Nt){
      for(k in 1:6){x_together[i,k+(j-1)*6] = y[j,i,k];}
      
      for(k in 1:3){mu_together[k+(j-1)*6] += ly[k]*eta[j,i,1];}
      for(k in 4:6){mu_together[k+(j-1)*6] += ly[k]*eta[j,i,2];}
    }
  }
  
  mu_together = mu_together/(N*Nt);
  /////////////////////////////////////////////////////////
  
  
  
  // this is the model-implied covariance matrix
    // general form for a cfa is cov=ly%*%Sigma%*%t(ly)+eps_M
    // with ly vector
    //      eps_M=diag(epsilon)
    //      Sigma=(I-B)^{-1}(G%*%Phi%*%G+Psi)((I-B)^{-1})^T
    //           =D^{-1}\Psi D^{-1}^T
    // with D=(I-B)
    //////////////////////////////////////////////////
    // add all repeating diagonal elements
    for(j in 1:Nt){//j<-1
      for(k in 1:6){
        epsm[k+(j-1)*6,k+(j-1)*6] += sigmaeps[k]*sigmaeps[k];
      }
    }

    //set all off-diagonals to zero
    //for(j in 1:(Nt*6)){//j<-1
    //  for(k in 1:(Nt*6)){
    //    if(k!=j){epsm[k,j] = 0}
    //  }
    //}
    //////////////////////////////////////////////////
    
    //////////////////////////////////////////////////
    //set factor loadings in right places
    //rest should be zeros from definition above
    for(j in 1:Nt){
      for(k in 1:3){lym[k+(j-1)*6,(j-1)*2+1] += ly[k];}
      for(k in 4:6){lym[k+(j-1)*6,(j-1)*2+2] += ly[k];}
    } 
    //////////////////////////////////////////////////
    
    
    //////////////////////////////////////////////////
    // build covariance matrix for eta
    //////////////////////////////////////////////////
    for(j in 1:(Nt-1)){//j<-2
      D[(j-1)*2+2+1,(j-1)*2+1] = -beta[1];
      D[(j-1)*2+2+2,(j-1)*2+2] = -beta[2];
    }
    
    for(j in 1:(2*Nt)){
      D[j,j] = 1;
      sigmazeta[j,j] = 1;
    }

    for(j in 1:Nt){//j<-1
      sigmazeta[(j-1)*2+2,(j-1)*2+1] = phi[2,1];
      sigmazeta[(j-1)*2+1,(j-1)*2+2] = phi[2,1];
    }

    
    sigmaeta = inverse(D) * sigmazeta * inverse(D)';
    
    //////////////////////////////////////////////////
    // build covariance matrix for y
    //////////////////////////////////////////////////
    sigmayhat = lym * sigmaeta * lym'+ epsm;
    
    
    
    //////////////////////////////////////////////////
    // loglikelihood
    //////////////////////////////////////////////////
    for (i in 1:N) {
    // defined above, or use from here as vector defined per person
    // within the 1:N loop; change x_together[i,] below to x_togeether
    //for(j in 1:Nt){for (k in 1:6){
    //  x_together_ml[j+(k-1)*6] = x[j,i,k]; 
    //}}
    

    log_lik[i] = 0;  
    //log_lik[i] = multi_normal_lpdf(x_together_ml | zero_vector , curr_cov_ll); // (mu_together_ml, curr_cov_ll) mu_zeros_ml
    log_lik[i] = multi_normal_lpdf(x_together[i,] | mu_together , sigmayhat); // (mu_together_ml, curr_cov_ll) mu_zeros_ml
    
    log_lik_sat[i] = 0;
    log_lik_sat[i] = multi_normal_lpdf(x_together[i,] | x_mean, x_cov);  // 
    ////////////
    }
    
    
    
  
}
