setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/MAIN")
getwd()


'
NEWEST!
'
# generate data for time series 
library(mvtnorm)
library(matrixcalc)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# install.packages("matrixcalc")
# setwd("C:\\holger\\SEM\\modelfit\\stanversion")



N <- 100  # persons
Nt <- 5 #time points


phi0 <- diag(2)*.5+.5 # cov(eta)
mu0  <- c(0,0)        #mean(eta)
ar0  <- c(.5,.5)      # ar(1) structure
ly0  <- matrix(c(1,1,1,0,0,0,
                 0,0,0,1,1,1),6,2,byrow=F) # factor loadings
td   <- diag(6)*.25 # cond. var(y) -> res var


# empty matrices for lvs and ovs
eta <- array(NA,c(N,Nt,2))

#matrix[N,6] y[Nt]; declared
#10,20,6 -> Nt, N, 6
y   <- array(NA,c(Nt,N,6))

# latent variables
# time point 1
eta[,1,] <- rmvnorm(N,mu0,phi0)
cov(eta[,1,]) # check cov
# rest
# NOTE: the total variance of the latent factors is currently not 1
for(j in 2:Nt){#j<-2
  zeta <- rmvnorm(N,c(0,0),phi0*(1-ar0^2)) # this is a residual with var (1-phi^2) sp tjat var(eta)==1
  eta[,j,1] <- mu0[1] + ar0[1]*eta[,j-1,1] + zeta[,1]
  eta[,j,2] <- mu0[2] + ar0[2]*eta[,j-1,2] + zeta[,2]
  #cov(zeta)
  #cov(eta[,2,]) # check cov
}

apply(eta[,,2],2,mean)

# observed data
for(j in 1:Nt){
  y[j,,] <- eta[,j,]%*%t(ly0)+rmvnorm(N,sigma = td)
}

###########################################################
# data analysis
# this is only copypaste from exercise 07 LDA

# data file
#x_mean <- c()
#dim(y)
#for(j in 1:Nt){
#  x_mean <- c(x_mean,apply(y[j,,],2,mean))
#}

# restructure data for overall covariance matrix
y0 <- matrix(NA,N,6*Nt)
for(i in 1:N){#i<-1
  for(j in 1:Nt){#j<-1
    y0[i,(j-1)*6+1:6] <- y[j,i,]  
  }
}
head(y0)
dim(y0)

y0 <- data.frame(y0)

cnom <- paste0("y",1:6,"t",1)
for(j in 2:Nt){
  cnom <- c(cnom,paste0("y",1:6,"t",j))
}
colnames(y0)<-cnom

x_cov  <- cov(y0)
dim(x_cov)
#round(eigen(x_cov)$values,3)
is.positive.definite(x_cov)

x_mean <- apply(y0,2,mean)

data1 <- list(y=y,Nt=Nt,N=N,x_mean=x_mean,x_cov=x_cov)

# read in input file (checks for mistakes)
rt1 <- stanc("dsem02.stan") 
# c++ compiler
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)
# fit the model to the data
# NOTE: currently 4 chains
fit1 <- sampling(sm1, data=data1)


# test if the model converged
# 1. Rhat (first general plot: should be smaller 1.1 [or 1.2])
stan_rhat(fit1)
# 2. trace plots (check if the lines overlap and when)
# You should investigate all parameters
stan_trace(fit1,"sigmaeps")
stan_trace(fit1,"beta")
stan_trace(fit1,"ly")
stan_trace(fit1,"phi")
# 3. density plots (check if the distributions are as expected (e.g. normal) and overlapping)
# Again, you should investigate all parameters
stan_dens(fit1,"Sigma",separate_chains =T)
stan_dens(fit1,"sigmahd",separate_chains =T)

# some other plots
#stan_plot(fit1,pars="phi")
#stan_plot(fit1, point_est = "mean", show_density = TRUE, fill_color = "maroon",pars="phi")
#stan_plot(fit1, point_est = "mean", show_density = TRUE, fill_color = "maroon",pars="lx")
#stan_hist(fit1,pars="lx")
#stan_diag(fit1)
#stan_par(fit1,par="lx[1]")
#stan_ac(fit1) # autocorrelation

# or use some nice plotting options here
#library(shinystan)
#launch_shinystan(fit1)


# extract the results in a table
params <- c("beta","ly","phi","sigmaeps")
print(fit1,pars=params)


print(fit1,pars="D")
print(fit1,pars="lym")
print(fit1,pars="epsm")
print(fit1,pars="sigmazeta")
print(fit1,pars="sigmaeta")

blub <- as.matrix(fit1, pars = c("sigmazeta"))
blub <- as.matrix(fit1, pars = c("sigmaeta"))
blub <- as.matrix(fit1, pars = c("D"))
blub2 <- apply(blub,2,mean)
Dmat <- matrix(blub2,Nt*2,Nt*2)

round(Dmat,2)

blub <- as.matrix(fit1, pars = c("sigmayhat"))
blub2 <- apply(blub,2,mean)
ymat <- matrix(blub2,Nt*6,Nt*6)

round(ymat-x_cov,2)

round(ymat,2)
round(x_cov,2)

blub <- as.matrix(fit1, pars = c("lym"))
blub2 <- apply(blub,2,mean)
lm <- matrix(blub2,Nt*6,Nt*2,byrow=F)
round(lm,2)

blub <- as.matrix(fit1, pars = c("mu_together"))
blub2 <- apply(blub,2,mean)

x_mean-blub2
x_mean

apply(as.matrix(fit1, pars = c("log_lik")),2,mean)

########################################################
########################################################
source("modelfitfuns.R")
library(blavaan)
pD = c("loo","waic","dic")
rescale = c("devm","ppmc")
fit.measures = "all"

########################################################
posterior_samples <- extract(fit1)

logliki <- data.frame(posterior_samples$log_lik)
ll <- mean(rowSums (logliki))

loglikisat <- data.frame(posterior_samples$log_lik_sat)
llsat <- mean(rowSums (loglikisat))

ll.blavaan <- mean(res1@external$samplls[,,1])

qqplot(rowSums(logliki), c(res1@external$samplls[,,1]))
abline(0,1)

plot(density(rowSums(logliki)),xlim=c(min(rowSums(logliki)),max(rowSums(logliki))),ylim=c(0,.2))
par(new=T)
plot(density(res1@external$samplls[,,1]),xlim=c(min(rowSums(logliki)),max(rowSums(logliki))),ylim=c(0,.2),lty=2)
abline(v=c(mean(rowSums(logliki)),mean(res1@external$samplls[,,1])))



################### Measures for STAN #####################
#loo(fit1)$p_loo[1]
fit.measures = "all"
my_array1 <- Create2DimMatrix(posterior_samples$log_lik, posterior_samples$log_lik_sat, nrow(logliki))
chains_number <- 4
my_array1 <- array(my_array1, dim = c(dim(my_array1)[1]/chains_number, chains_number, 2)) ###### TRNASFORM BASED ON CHAINS
chisqs1 <- as.numeric(apply(my_array1, 2,
                            function(x) 2*(x[,2] - x[,1]))) 

a1 <- BayesChiFit(obs = chisqs1, 
                  nvar = 6*Nt, pD = loo(fit1)$p_loo[1],
                  N = N,
                  fit.measures = fit.measures, ms = FALSE, null_model = FALSE)

out2 <- c(mean(a1@indices$BRMSEA),
          mean(a1@indices$BGammaHat),
          mean(a1@indices$adjBGammaHat),
          mean(a1@indices$BMc))

out2

########################################################
########################################################
# lavaan code
library(blavaan)
library(lavaan)

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
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
y1t5 ~~ td1*y1t5
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
y2t5 ~~ td2*y2t5
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
y3t5 ~~ td3*y3t5
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
y4t5 ~~ td4*y4t5
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
y5t5 ~~ td5*y5t5
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
y6t5 ~~ td6*y6t5
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
eta1t5 ~ beta1*eta1t4
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
eta2t5 ~ beta2*eta2t4
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
eta1t5 ~~ phi21*eta2t5
'


res1 <- bsem(dsem, data=y0,
             burnin = 1000, n.chains = 4, sample = 1000, 
             mcmcfile = "model1", std.lv = TRUE,
             dp=dpriors(lambda="normal(0,1)"))

summary(res1)

blavFitIndices(res1)

chisqs <- as.numeric(apply(res1@external$samplls, 2,
                           function(x) 2*(x[,2] - x[,1])))   
fit_pd <- fitMeasures(res1, paste0('p_', pD))              

a0 <- BayesChiFit(obs = chisqs, 
                  nvar = 6*Nt, pD = fit_pd[1],
                  N = N,
                  fit.measures = fit.measures, ms = FALSE, null_model = FALSE)#@details


out1 <- c(mean(a0@indices$BRMSEA),
          mean(a0@indices$BGammaHat),
          mean(a0@indices$adjBGammaHat),
          mean(a0@indices$BMc))

out1

posterior_samples$log_lik
sum(posterior_samples$log_lik_sat)
########################################################
### WORKING ON FUNCTIONS!!!
########################################################


chisqs1 <- as.numeric(apply(my_array1, 2,
                            function(x) 2*(x[,2] - x[,1]))) 



dim(chisqs1)
dim(chisqs)
summary(chisqs)
chisqs


res1@external$samplls[,,1]

ddd_sat <- data.frame(posterior_samples$log_lik)
ddd_sat
ddd_sat <- ddd_sat %>%
  mutate(ll = rowSums(., na.rm=TRUE))
ddd_sat <- data.frame(res1@external$samplls[,,1] )%>%
  mutate(ll = rowSums(., na.rm=TRUE))
# -3275.492
# -3031.602 vs -3031.526, Very close!
##########################
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




