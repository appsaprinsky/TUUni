# setwd("C:\\holger\\SEM\\modelfit\\stanversion")      HB!!!!!!!!!!!!!
setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/FINAL")
getwd()
# Sample size N = 25, 50, 75, 150, 300, 500
# Number of time points: Nt = 1, 2, 3, 4, 5, 10, 15
library(mvtnorm)
library(matrixcalc)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(blavaan)



######## FUNCTIONS IMPORTANT
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

###################################

rt1 <- stanc("dsem02.stan") 
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)

36*2

person_size_SIMULATE <- c(25, 50, 75, 150, 300, 500) # N
time_point_SIMULATE <- c(2, 3, 4, 5, 10, 15) # Nt
model_TRUE_MISS_SIMULATE  <- c(0, 0.3) 

##### SHORTER #####
person_size_SIMULATE <- c(25, 50, 75) # N
time_point_SIMULATE <- c(2, 3, 4) # Nt


# population parameters
# constant across conditins
phi0 <- diag(2)*.5+.5 # cov(eta)
mu0  <- c(0,0)        #mean(eta)
ar0  <- c(.5,.5)      # ar(1) structure
ly0  <- matrix(c(1,1,1,0,0,0,
                 0,0,0,1,1,1),6,2,byrow=F) # factor loadings
td   <- diag(6)*.25 # cond. var(y) -> res var

global_SIMULATE_Info <- data.frame(
  PersonSize = numeric(),
  TimePoint = numeric(), 
  ModelTorF = numeric(), 
  BRMSEA = numeric(),
  BGammaHat = numeric(),
  adjBGammaHat = numeric(),
  BMc = numeric()
)
is_positive_def <- FALSE
rowCheck <- 1
pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"

#HB: Set.seed for replicability
set.seed(131212023)

for (model_TRUE_MISS in model_TRUE_MISS_SIMULATE){
  for (person_size in person_size_SIMULATE) {
    for (time_point in time_point_SIMULATE) {
      
      is_positive_def <- FALSE
      

      N <- person_size  # persons
      Nt <- time_point #time points
      ly1  <- matrix(c(1,1,1,0,0,model_TRUE_MISS,
                       0,0,0,1,1,1),6,2,byrow=F) # factor loadings
      
      #######################################

      eta <- array(NA,c(N,Nt,2))
      y   <- array(NA,c(Nt,N,6))
      eta[,1,] <- rmvnorm(N,mu0,phi0)
      for(j in 2:Nt){#j<-2
        zeta <- rmvnorm(N,c(0,0),phi0*(1-ar0^2)) 
        eta[,j,1] <- mu0[1] + ar0[1]*eta[,j-1,1] + zeta[,1]
        eta[,j,2] <- mu0[2] + ar0[2]*eta[,j-1,2] + zeta[,2]
      }
      
      for(j in 1:Nt){
        y[j,,] <- eta[,j,]%*%t(ly1)+rmvnorm(N,sigma = td)
      }
      
      y0 <- matrix(NA,N,6*Nt)
      for(i in 1:N){#i<-1
        for(j in 1:Nt){#j<-1
          y0[i,(j-1)*6+1:6] <- y[j,i,]  
        }
      }
      y0 <- data.frame(y0)
      cnom <- paste0("y",1:6,"t",1)
      for(j in 2:Nt){
        cnom <- c(cnom,paste0("y",1:6,"t",j))
      }
      colnames(y0)<-cnom
      x_cov  <- cov(y0)
      x_mean <- apply(y0,2,mean)
      is_positive_def <- is.positive.definite(x_cov)
      is_positive_def
      
      print(is_positive_def)
      #}
      #######################################
      #######################################
      
      if(is_positive_def==TRUE){
        
        global_SIMULATE_Info[rowCheck, 1] <- person_size
        global_SIMULATE_Info[rowCheck, 2] <- time_point
        global_SIMULATE_Info[rowCheck, 3] <- model_TRUE_MISS      
        data1 <- list(y=y,Nt=Nt,N=N,x_mean=x_mean,x_cov=x_cov)
        
        
        fit1 <- sampling(sm1, data=data1)
        posterior_samples <- extract(fit1)
        ###########################
        my_array1 <- Create2DimMatrix(posterior_samples$log_lik, posterior_samples$log_lik_sat, nrow(posterior_samples$log_lik))#typo
        chains_number <- 4
        my_array1 <- array(my_array1, dim = c(dim(my_array1)[1]/chains_number, chains_number, 2)) ###### TRNASFORM BASED ON CHAINS
        chisqs1 <- as.numeric(apply(my_array1, 2,
                                    function(x) 2*(x[,2] - x[,1]))) 
        a1 <- BayesChiFit(obs = chisqs1, 
                          nvar = 6*Nt, pD = loo(fit1)$p_loo[1],
                          N = N,
                          fit.measures = fit.measures, ms = FALSE, null_model = FALSE)
        global_SIMULATE_Info[rowCheck, 4] <- mean(a1@indices$BRMSEA)
        global_SIMULATE_Info[rowCheck, 5] <- mean(a1@indices$BGammaHat)
        global_SIMULATE_Info[rowCheck, 6] <- mean(a1@indices$adjBGammaHat)
        global_SIMULATE_Info[rowCheck, 7] <- mean(a1@indices$BMc)
        write.csv(global_SIMULATE_Info, file = "global_SIMULATE_Info.csv", row.names = FALSE)  # Write constantly to monitor the output
        rowCheck = rowCheck+1
      }#end of if-loop
      
      # name_local_SIMULATE_Info <- paste("local", as.character(person_size), as.character(time_point), as.character(model_TRUE_MISS), ".csv", sep = "_")
      # write.csv(local_SIMULATE_Info, file = name_local_SIMULATE_Info, row.names = FALSE) # Save samples run

    }
  }
}


