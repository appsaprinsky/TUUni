setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/SIMULATION")
getwd()
# Sample size N = 25, 50, 75
# Number of time points: Nt = 5, 10, 25
# For each condition run R = 100 samples



library(mvtnorm)
library(matrixcalc)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



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


# data_generate <- function(x, y) {
#   result <- x * y
#   return(result)
# }

###################################





















person_size_SIMULATE <- c(25, 50, 75) # N
time_point_SIMULATE <- c(5, 10, 25) # Nt
run_Samples_SIMULATE <- 100
model_TRUE_MISS <- c(0, 0.3) 



##### For shorter check
person_size_SIMULATE <- c(75) # N
time_point_SIMULATE <- c(10) # Nt
run_Samples_SIMULATE <- 2
model_TRUE_MISS_SIMULATE <- c(0, 0.3) 

global_SIMULATE_Info <- data.frame(
  PersonSize = numeric(),
  TimePoint = numeric(), 
  ModelTorF = numeric(), 
  RunSamples = numeric(),
  BRMSEA = numeric(),
  BGammaHat = numeric(),
  adjBGammaHat = numeric(),
  BMc = numeric()
)


# (75 and 10 ) is ok
# (75 and 10 ) is ok

is_positive_def <- FALSE
rowCheck <- 1

for (model_TRUE_MISS in model_TRUE_MISS_SIMULATE){
  for (person_size in person_size_SIMULATE) {
    for (time_point in time_point_SIMULATE) {
      is_positive_def <- FALSE
      global_SIMULATE_Info[rowCheck, 1] <- person_size
      global_SIMULATE_Info[rowCheck, 2] <- time_point
      global_SIMULATE_Info[rowCheck, 3] <- model_TRUE_MISS      
      global_SIMULATE_Info[rowCheck, 4] <- run_Samples_SIMULATE
      
      while (is_positive_def != TRUE) {   # Rerun data creation part until cov matrix is positive definite
        N <- person_size  # persons
        Nt <- time_point #time points
        phi0 <- diag(2)*.5+.5 # cov(eta)
        mu0  <- c(0,0)        #mean(eta)
        ar0  <- c(.5,.5)      # ar(1) structure
        ly0  <- matrix(c(1,1,1,0,0,0,
                         0,0,0,1,1,1),6,2,byrow=F) # factor loadings
        ly1  <- matrix(c(1,1,1,0,0,model_TRUE_MISS,
                         0,0,0,1,1,1),6,2,byrow=F) # factor loadings
        td   <- diag(6)*.25 # cond. var(y) -> res var
        # empty matrices for lvs and ovs
        eta <- array(NA,c(N,Nt,2))
        #matrix[N,6] y[Nt]; declared
        y   <- array(NA,c(Nt,N,6))
        # latent variables
        # time point 1
        eta[,1,] <- rmvnorm(N,mu0,phi0)
        cov(eta[,1,]) # check cov
        # NOTE: the total variance of the latent factors is currently not 1
        for(j in 2:Nt){#j<-2
          zeta <- rmvnorm(N,c(0,0),phi0*(1-ar0^2)) # this is a residual with var (1-phi^2) sp tjat var(eta)==1
          eta[,j,1] <- mu0[1] + ar0[1]*eta[,j-1,1] + zeta[,1]
          eta[,j,2] <- mu0[2] + ar0[2]*eta[,j-1,2] + zeta[,2]
        }
        # observed data
        
        for(j in 1:5){
          y[j,,] <- eta[,j,]%*%t(ly1)+rmvnorm(N,sigma = td)
        }
        if (time_point >= 6){          # If time points is 5, we do not run this part
          for(j in 6:Nt){
            y[j,,] <- eta[,j,]%*%t(ly0)+rmvnorm(N,sigma = td)
          }
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
        dim(x_cov)
        #round(eigen(x_cov)$values,3)
        x_mean <- apply(y0,2,mean)
        is_positive_def <- is.positive.definite(x_cov)
        print(is_positive_def)
      }
      
      pD = c("loo","waic","dic")
      rescale = c("devm","ppmc")
      fit.measures = "all"
      data1 <- list(y=y,Nt=Nt,N=N,x_mean=x_mean,x_cov=x_cov)
      rt1 <- stanc("dsem02.stan") 
      sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)
      
      
      
      local_SIMULATE_Info <- data.frame(
        BRMSEA = numeric(),
        BGammaHat = numeric(),
        adjBGammaHat = numeric(),
        BMc = numeric()
      )
      for (SAMPLING in 1:run_Samples_SIMULATE){ ### All Samples run
        fit1 <- sampling(sm1, data=data1)
        posterior_samples <- extract(fit1)
        my_array1 <- Create2DimMatrix(posterior_samples$log_lik, posterior_samples$log_lik_sat, nrow(logliki))
        chains_number <- 4
        my_array1 <- array(my_array1, dim = c(dim(my_array1)[1]/chains_number, chains_number, 2)) ###### TRNASFORM BASED ON CHAINS
        chisqs1 <- as.numeric(apply(my_array1, 2,
                                    function(x) 2*(x[,2] - x[,1]))) 
        a1 <- BayesChiFit(obs = chisqs1, 
                          nvar = 6*Nt, pD = loo(fit1)$p_loo[1],
                          N = N,
                          fit.measures = fit.measures, ms = FALSE, null_model = FALSE)
        # out2 <- c(mean(a1@indices$BRMSEA),
        #           mean(a1@indices$BGammaHat),
        #           mean(a1@indices$adjBGammaHat),
        #           mean(a1@indices$BMc))
        local_SIMULATE_Info[SAMPLING,1] <- mean(a1@indices$BRMSEA)
        local_SIMULATE_Info[SAMPLING,2] <- mean(a1@indices$BGammaHat)
        local_SIMULATE_Info[SAMPLING,3] <- mean(a1@indices$adjBGammaHat)
        local_SIMULATE_Info[SAMPLING,4] <- mean(a1@indices$BMc)
      }
      name_local_SIMULATE_Info <- paste("local", as.character(person_size), as.character(time_point), as.character(model_TRUE_MISS), ".csv", sep = "_")
      write.csv(local_SIMULATE_Info, file = name_local_SIMULATE_Info, row.names = FALSE) # Save samples run
      
      global_SIMULATE_Info[rowCheck, 5] <- mean(local_SIMULATE_Info$BRMSEA)
      global_SIMULATE_Info[rowCheck, 6] <- mean(local_SIMULATE_Info$BGammaHat)
      global_SIMULATE_Info[rowCheck, 7] <- mean(local_SIMULATE_Info$adjBGammaHat)     
      global_SIMULATE_Info[rowCheck, 8] <- mean(local_SIMULATE_Info$BMc)
      
      write.csv(global_SIMULATE_Info, file = "global_SIMULATE_Info.csv", row.names = FALSE)  # Write constantly to monitor the output
      rowCheck = rowCheck+1
      
    }
  }
}









# # test if the model converged
# # 1. Rhat (first general plot: should be smaller 1.1 [or 1.2])
# stan_rhat(fit1)
# # 2. trace plots (check if the lines overlap and when)
# # You should investigate all parameters
# stan_trace(fit1,"sigmaeps")
# stan_trace(fit1,"beta")
# stan_trace(fit1,"ly")
# stan_trace(fit1,"phi")
# # 3. density plots (check if the distributions are as expected (e.g. normal) and overlapping)
# # Again, you should investigate all parameters
# stan_dens(fit1,"Sigma",separate_chains =T)
# stan_dens(fit1,"sigmahd",separate_chains =T)
# # extract the results in a table
# params <- c("beta","ly","phi","sigmaeps")
# print(fit1,pars=params)
# print(fit1,pars="D")
# print(fit1,pars="lym")
# print(fit1,pars="epsm")
# print(fit1,pars="sigmazeta")
# print(fit1,pars="sigmaeta")
# blub <- as.matrix(fit1, pars = c("sigmazeta"))
# blub <- as.matrix(fit1, pars = c("sigmaeta"))
# blub <- as.matrix(fit1, pars = c("D"))
# blub2 <- apply(blub,2,mean)
# Dmat <- matrix(blub2,Nt*2,Nt*2)
# round(Dmat,2)
# blub <- as.matrix(fit1, pars = c("sigmayhat"))
# blub2 <- apply(blub,2,mean)
# ymat <- matrix(blub2,Nt*6,Nt*6)
# round(ymat-x_cov,2)
# round(ymat,2)
# round(x_cov,2)
# blub <- as.matrix(fit1, pars = c("lym"))
# blub2 <- apply(blub,2,mean)
# lm <- matrix(blub2,Nt*6,Nt*2,byrow=F)
# round(lm,2)
# blub <- as.matrix(fit1, pars = c("mu_together"))
# blub2 <- apply(blub,2,mean)
# x_mean-blub2
# x_mean
# apply(as.matrix(fit1, pars = c("log_lik")),2,mean)

########################################################
########################################################
# logliki <- data.frame(posterior_samples$log_lik)
# ll <- mean(rowSums (logliki))
# 
# loglikisat <- data.frame(posterior_samples$log_lik_sat)
# llsat <- mean(rowSums (loglikisat))
# 
# ll.blavaan <- mean(res1@external$samplls[,,1])
# 
# qqplot(rowSums(logliki), c(res1@external$samplls[,,1]))
# abline(0,1)
# 
# plot(density(rowSums(logliki)),xlim=c(min(rowSums(logliki)),max(rowSums(logliki))),ylim=c(0,.2))
# par(new=T)
# plot(density(res1@external$samplls[,,1]),xlim=c(min(rowSums(logliki)),max(rowSums(logliki))),ylim=c(0,.2),lty=2)
# abline(v=c(mean(rowSums(logliki)),mean(res1@external$samplls[,,1])))
# 


################### Measures for STAN #####################
#loo(fit1)$p_loo[1]

