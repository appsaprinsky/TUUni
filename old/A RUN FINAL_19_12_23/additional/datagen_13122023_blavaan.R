setwd("C:\\holger\\SEM\\modelfit\\stanversion")
#setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/SIMULATION")
getwd()
# Sample size N = 25, 50, 75
# Number of time points: Nt = 5, 10, 25
# For each condition run R = 100 samples

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


# load model once here
rt1 <- stanc("dsem02_.stan") 
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)

rt2 <- stanc("dsem02_chol.stan") 
sm2 <- stan_model(stanc_ret = rt2, verbose=FALSE)



###################################
###    SIMULATION STARTS HERE
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

# these are global variables for the condition
pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"
# HB: Do not reload the model each time in the loop
#rt1 <- stanc("dsem02.stan") 
#sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)



is_positive_def <- FALSE
rowCheck <- 1



dsem <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
eta1t4 =~ ly1*y1t4 + ly2*y2t4 + ly3*y3t4
eta2t4 =~ ly4*y4t4 + ly5*y5t4 + ly6*y6t4
#
eta1t5 =~ ly1*y1t5 + ly2*y2t5 + ly3*y3t5
eta2t5 =~ ly4*y4t5 + ly5*y5t5 + ly6*y6t5
#
eta1t6 =~ ly1*y1t6 + ly2*y2t6 + ly3*y3t6
eta2t6 =~ ly4*y4t6 + ly5*y5t6 + ly6*y6t6
#
eta1t7 =~ ly1*y1t7 + ly2*y2t7 + ly3*y3t7
eta2t7 =~ ly4*y4t7 + ly5*y5t7 + ly6*y6t7
#
eta1t8 =~ ly1*y1t8 + ly2*y2t8 + ly3*y3t8
eta2t8 =~ ly4*y4t8 + ly5*y5t8 + ly6*y6t8
#
eta1t9 =~ ly1*y1t9 + ly2*y2t9 + ly3*y3t9
eta2t9 =~ ly4*y4t9 + ly5*y5t9 + ly6*y6t9
#
eta1t10 =~ ly1*y1t10 + ly2*y2t10 + ly3*y3t10
eta2t10 =~ ly4*y4t10 + ly5*y5t10 + ly6*y6t10
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
y1t5 ~~ td1*y1t5
y1t6 ~~ td1*y1t6
y1t7 ~~ td1*y1t7
y1t8 ~~ td1*y1t8
y1t9 ~~ td1*y1t9
y1t10 ~~ td1*y1t10
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
y2t5 ~~ td2*y2t5
y2t6 ~~ td2*y2t6
y2t7 ~~ td2*y2t7
y2t8 ~~ td2*y2t8
y2t9 ~~ td2*y2t9
y2t10 ~~ td2*y2t10
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
y3t5 ~~ td3*y3t5
y3t6 ~~ td3*y3t6
y3t7 ~~ td3*y3t7
y3t8 ~~ td3*y3t8
y3t9 ~~ td3*y3t9
y3t10 ~~ td3*y3t10
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
y4t5 ~~ td4*y4t5
y4t6 ~~ td4*y4t6
y4t7 ~~ td4*y4t7
y4t8 ~~ td4*y4t8
y4t9 ~~ td4*y4t9
y4t10 ~~ td4*y4t10
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
y5t5 ~~ td5*y5t5
y5t6 ~~ td5*y5t6
y5t7 ~~ td5*y5t7
y5t8 ~~ td5*y5t8
y5t9 ~~ td5*y5t9
y5t10 ~~ td5*y5t10
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
y6t5 ~~ td6*y6t5
y6t6 ~~ td6*y6t6
y6t7 ~~ td6*y6t7
y6t8 ~~ td6*y6t8
y6t9 ~~ td6*y6t9
y6t10 ~~ td6*y6t10
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
eta1t5 ~ beta1*eta1t4
eta1t6 ~ beta1*eta1t5
eta1t7 ~ beta1*eta1t6
eta1t8 ~ beta1*eta1t7
eta1t9 ~ beta1*eta1t8
eta1t10 ~ beta1*eta1t9
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
eta2t5 ~ beta2*eta2t4
eta2t6 ~ beta2*eta2t5
eta2t7 ~ beta2*eta2t6
eta2t8 ~ beta2*eta2t7
eta2t9 ~ beta2*eta2t8
eta2t10 ~ beta2*eta2t9
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
eta1t5 ~~ phi21*eta2t5
eta1t6 ~~ phi21*eta2t6
eta1t7 ~~ phi21*eta2t7
eta1t8 ~~ phi21*eta2t8
eta1t9 ~~ phi21*eta2t9
eta1t10 ~~ phi21*eta2t10
'





#HB: Set.seed for replicability
set.seed(131212023)


for (model_TRUE_MISS in model_TRUE_MISS_SIMULATE){
  for (person_size in person_size_SIMULATE) {
    for (time_point in time_point_SIMULATE) {
      
      ######################################
      # check for now
      model_TRUE_MISS <- .3#model_TRUE_MISS_SIMULATE[1]
      person_size <- 75#person_size_SIMULATE[1]
      time_point <- 10#time_point_SIMULATE[1]
      ######################################
      
      is_positive_def <- FALSE
      global_SIMULATE_Info[rowCheck, 1] <- person_size
      global_SIMULATE_Info[rowCheck, 2] <- time_point
      global_SIMULATE_Info[rowCheck, 3] <- model_TRUE_MISS      
      global_SIMULATE_Info[rowCheck, 4] <- run_Samples_SIMULATE
      
      
      #HB: this goes before the sampling loop
      local_SIMULATE_Info <- data.frame(
        BRMSEA = numeric(),
        BGammaHat = numeric(),
        adjBGammaHat = numeric(),
        BMc = numeric()
      )
      
      #while (is_positive_def != TRUE) {   # Rerun data creation part until cov matrix is positive definite
      # HB: I think the pos. def. is a matter of identification depending on no. of time points, 
      # I dont think looping through that makes much sense -> also this while loop will run forever, you need a break of the loop
      
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
      
      #######################################
      # start generating data here (the rest is constant within the condition)
      
      for (SAMPLING in 1:run_Samples_SIMULATE){ ### All Samples run !!!!!!!!!!!!!!!!!
        #SAMPLING <- 1
        
        # empty matrices for lvs and ovs
        eta <- array(NA,c(N,Nt,2))
        #matrix[N,6] y[Nt]; declared
        y   <- array(NA,c(Nt,N,6))
        # latent variables
        # time point 1
        eta[,1,] <- rmvnorm(N,mu0,phi0)
        #cov(eta[,1,]) # check cov
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
        #dim(x_cov)
        #round(eigen(x_cov)$values,3)
        x_mean <- apply(y0,2,mean)
        is_positive_def <- is.positive.definite(x_cov)
        is_positive_def
        
        #  print(is_positive_def)
        #}
        #######################################
        #######################################
        
        if(is_positive_def==TRUE){
          
          a0<- Sys.time()
          res1 <- bsem(dsem, data=y0,
                       burnin = 1000, n.chains = 4, sample = 1000, 
                       #mcmcfile = "model1", 
                       std.lv = TRUE)#,
                       #dp=dpriors(lambda="normal(0,1)"))
          #summary(res1)
          a1 <- blavFitIndices(res1)
          
          #chisqs <- as.numeric(apply(res1@external$samplls, 2,
          #                           function(x) 2*(x[,2] - x[,1])))   
          #fit_pd <- fitMeasures(res1, paste0('p_', pD))              
          #  
          #a0 <- BayesChiFit(obs = chisqs, 
          #                  nvar = 6*Nt, pD = fit_pd[1],
          #                  N = N,
          #                  fit.measures = fit.measures, ms = FALSE, null_model = FALSE)#@details
          a0[2]<- Sys.time()
          a0[2]-a0[1]
          
          local_SIMULATE_Info[SAMPLING,1] <- mean(a1@indices$BRMSEA)
          local_SIMULATE_Info[SAMPLING,2] <- mean(a1@indices$BGammaHat)
          local_SIMULATE_Info[SAMPLING,3] <- mean(a1@indices$adjBGammaHat)
          local_SIMULATE_Info[SAMPLING,4] <- mean(a1@indices$BMc)
        }#end of if-loop
        
        name_local_SIMULATE_Info <- paste("local", as.character(person_size), as.character(time_point), as.character(model_TRUE_MISS), ".csv", sep = "_")
        write.csv(local_SIMULATE_Info, file = name_local_SIMULATE_Info, row.names = FALSE) # Save samples run
      }# end of sampling
        
        
        
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

