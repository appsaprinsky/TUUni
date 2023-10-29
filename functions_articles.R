
library(blavaan)
library(semTools)
library(loo)

postpred_mgv <- function(lavpartable, lavmodel, lavoptions, 
                     lavsamplestats, lavdata, lavcache, lavjags,
                     samplls, Fit, measure = "logl", thin = 5) { 
  #### MGV: added Fit to extract fx for the effective number of parameters for the DIC
  ## if we keep the pD from LOO keeping Fit is not neccesary
  
  ## run through lavjags$mcmc, generate data from various posterior
  ## samples. thin like we do in samp_lls
  lavmcmc <- blavaan:::make_mcmc(lavjags)
  samp.indices <- blavaan:::sampnums(lavjags, thin=thin)
  n.chains <- length(lavmcmc)
  psamp <- length(samp.indices)
  
  ## parallel across chains if we can
  ncores <- NA
  loop.comm <- "lapply"
  if(.Platform$OS.type != "windows" & requireNamespace("parallel", quietly = TRUE)){
    ncores <- min(n.chains, parallel::detectCores())
    loop.comm <- "mclapply"
  }
  
  origlavmodel <- lavmodel
  origlavdata <- lavdata
  
  loop.args <- list(X = 1:n.chains, FUN = function(j){
    ### MGV: added chi-square boots to later calculate relative fits
    ind <- csdist <- csboots <- rep(NA, psamp)
    for(i in 1:psamp){
      ## translate each posterior sample to a model-implied mean vector +
      ## cov matrix.
      lavmodel <- blavaan:::fill_params(lavmcmc[[j]][samp.indices[i],],
                              origlavmodel, lavpartable)
      
      ## generate data (some code from lav_bootstrap.R)
      implied <- lav_model_implied(lavmodel)
      Sigma.hat <- implied$cov
      Mu.hat <- implied$mean
      dataeXo <- lavdata@eXo
      
      ## TODO? this generates complete cases; maybe we want missing
      ## observations to stay missing in the generated data:
      dataX <- vector("list", length=lavdata@ngroups)
      for(g in 1:lavsamplestats@ngroups) {
        dataX[[g]] <- MASS::mvrnorm(n     = lavsamplestats@nobs[[g]],
                                    Sigma = Sigma.hat[[g]],
                                    mu    = Mu.hat[[g]])
        dataX[[g]][is.na(origlavdata@X[[g]])] <- NA
      }
      
      ## compute (i) X2 of generated data and model-implied
      ## moments, along with (ii) X2 of real data and model-implied
      ## moments.
      chisq.obs <- -2*(samplls[i, j, 1] -
                         samplls[i, j, 2])
      #get_ll(lavmodel = lavmodel,
      #    lavpartable = lavpartable,
      #    lavsamplestats = lavsamplestats,
      #    lavoptions = lavoptions,
      #    lavcache = lavcache,
      #    lavdata = origlavdata,
      #    measure = measure)
      
      ## check for missing, to see if we can easily get baseline ll for chisq
      mis <- FALSE
      if(any(is.na(unlist(lavdata@X)))) mis <- TRUE
      
      if(!mis){
        lavdata@X <- dataX
        
        chisq.boot <- 2*diff(blavaan:::get_ll(lavmodel = lavmodel,
                                    lavsamplestats = lavsamplestats,
                                    lavdata = lavdata,
                                    measure = measure))
      } else {
        ## we need lavaan to get the saturated log-l for missing data (EM)
        
        # YR: ugly hack to avoid lav_samplestats_from_data:
        # reconstruct data + call lavaan()
        # ed: if we need lavaan() anyway, might as well
        # get the chisq while we're here:
        DATA.X <- do.call("rbind", dataX)
        colnames(DATA.X) <- lavdata@ov.names[[1L]]
        DATA.eXo <- do.call("rbind", dataeXo)
        if(!is.null(DATA.eXo)) {
          colnames(DATA.eXo) <- lavdata@ov.names.x[[1L]]
          DATA <- cbind(DATA.X, DATA.eXo)
        } else {
          DATA <- DATA.X
        }
        DATA <- as.data.frame(DATA)
        
        lavoptions2 <- lavoptions
        lavoptions2$verbose <- FALSE
        lavoptions2$estimator <- "ML"
        lavoptions2$se <- "none"
        lavoptions2$test <- "standard"
        lavoptions2$optim.method <- "none"
        lavmodel2 <- lavmodel
        if("control" %in% slotNames(lavmodel2)){
          lavmodel2@control <- list(optim.method="none")
        }
        if(lavsamplestats@ngroups > 1L) {
          DATA$.g. <- rep(1:lavdata@ngroups, 
                          times = unlist(lavdata@nobs))
          out <- lavaan(slotOptions = lavoptions2, 
                        slotParTable = lavpartable,
                        slotSampleStats = NULL, slotData = NULL, 
                        slotModel = lavmodel2, slotCache = lavcache, 
                        data = DATA, group = ".g.")
        } else {
          out <- lavaan(slotOptions = lavoptions2, 
                        slotParTable = lavpartable,
                        slotSampleStats = NULL, slotData = NULL, 
                        slotModel = lavmodel2, slotCache = lavcache, 
                        data = DATA)
        }
        # bootSampleStats <- out@SampleStats
        # end of ugly hack
        
        if(measure %in% c("logl", "chisq")){
          chisq.boot <- fitMeasures(out, "chisq")
        } else {
          chisq.boot <- fitMeasures(out, measure)
        }
        
        ## see lines 286-298 of lav_bootstrap to avoid fixed.x errors?
        ## chisq.boot <- 2*diff(get_ll(lavmodel = lavmodel,
        ##                             lavpartable = lavpartable,
        ##                             lavsamplestats = bootSampleStats,
        ##                             lavoptions = lavoptions,
        ##                             lavcache = lavcache,
        ##                             lavdata = lavdata,
        ##                             measure = measure))
      }
      ## record whether observed value is larger
      ind[i] <- chisq.obs < chisq.boot
      csdist[i] <- chisq.obs
      csboots[i] <- chisq.boot
      ## MGV: saved boots to calculate relative fits
    } # i
    
    list(ind = ind, csdist = csdist, csboots = csboots)
  })
  
  if(loop.comm == "mclapply"){
    loop.args <- c(loop.args, list(mc.cores = ncores))
    res <- do.call(parallel::mclapply, loop.args)
  } else {
    res <- do.call(lapply, loop.args)
  }
  
  ind <- unlist(lapply(res, function(x) x$ind))
  csdist <- unlist(lapply(res, function(x) x$csdist))
  csboots <- unlist(lapply(res, function(x) x$csboots))
  
  ## pD from LOO
  ## MGV: I am recalculating the casell to get the pD from loo, is there a way to save this time?
  ## also have the calculation from the pD from DIC, 
  ## is it worth to keep this or is the DIC pD good enough?
  pd <- loo(blavaan:::case_lls(lavjags, origlavmodel, lavpartable, 
                               lavsamplestats, lavoptions, lavcache, 
                               origlavdata))$p_loo
  ## pD from DIC
  ## got this from blav_fit_measures
  #pd <- 2*(Fit@fx - mean(as.numeric(samplls[,,1])))
  
  ## MGV: this function calculates rmsea, gammahat, adjusted gammahat
  ## set to export all the estimates. To later saved mean and quantiles of interest
  
  ppval <- mean(as.numeric(ind))
  cspi <- quantile(as.numeric(csdist), c(.025,.975))
  ## MGV: estimate the mean and 90% CI for RMSEA, gammahat, adjusted gammahat
  
  list(ppval=ppval, cspi=cspi, pd=pd, chisqs=cbind(obs=csdist, reps=csboots))
}



### Bayesian RMSEA from Rens
### MGV: modified to also calculate gammahat, adjusted gammahat
### TDJ: added McDonald's centrality index

### Bayesian RMSEA from Rens
### MGV: modified to also calculate gammahat, adjusted gammahat
### TDJ: added McDonald's centrality index
### if a null model information provided, it calculates CFI, TLI, NFI
BayesRelFit <-function(obs, rep, nvar, pD, N, ms = TRUE, Min1 = TRUE, Ngr = 1,
                       null_model=TRUE, obs_null=NULL, rep_null=NULL, pD_null=NULL){
  
  # # Compute number of parameters
  if(ms) p <- (((nvar * (nvar + 1)) / 2) + nvar)
  if(!ms) p <- (((nvar * (nvar + 1))/ 2) + 0)
  p <- p * Ngr
  # # Substract parameters and estimated parameters
  dif.ppD <- p - pD
  nonc <- ( ( obs-pD ) - dif.ppD )
  # # Correct if numerator is smaller than zero
  nonc[nonc < 0] <- 0
  # # Compute BRMSEA (with or without the -1 correction)
  if(Min1)
    BRMSEA <- sqrt(nonc / (dif.ppD * (N -1)))*sqrt(Ngr)
  if(!Min1) BRMSEA <- sqrt(nonc / (dif.ppD * N ))*sqrt(Ngr)
  
  ## compute GammaHat and adjusted GammaHat
  gammahat <- nvar / ( nvar+2* ((nonc)/(N-1))  )
  adjgammahat <- 1 - (((Ngr * nvar * (nvar + 1))/2)/dif.ppD) * (1 - gammahat)
  
  ## compute McDonald's centrality index
  Mc <- exp(-.5 * nonc/(N-1) )
  
  ## calculate fit when null model is provided
  if (null_model) {
    dif.ppD_null <- p - pD_null
    nonc_null <- ( ( obs_null-pD_null ) - dif.ppD_null )
    
    cfi <- (nonc_null - nonc)/nonc_null 
    tli <- ((( obs_null-pD_null )/dif.ppD_null) - (( obs-pD )/dif.ppD)) / (((( obs_null-pD_null )/dif.ppD_null))-1) 
    nfi <- (( obs_null-pD_null ) - ( obs-pD )) / ( obs_null-pD_null )
    
    out <- cbind(BRMSEA=BRMSEA, BGammaHat=gammahat, adjBGammaHat=adjgammahat, BMc = Mc, BCFI=cfi, BTLI=tli, BNFI=nfi)
  } else {
    out <- cbind(BRMSEA=BRMSEA, BGammaHat=gammahat, adjBGammaHat=adjgammahat, BMc = Mc)
  }
  
  return(out)
}




## fit@ParTable
## fit@Model
## fit@Options
## fit@SampleStats
## fit@Data
## fit@Cache
## fit@external$mcmcout
## fit@external$samplls
## fit@Fit


## single group model
#temp1g <- postpred_mgv(lavpartable=fit1g@ParTable, lavmodel=fit1g@Model, 
#                     lavoptions=fit1g@Options, lavsamplestats=fit1g@SampleStats, 
#                     lavdata=fit1g@Data, lavcache=fit1g@Cache, 
#                     lavjags=fit1g@external$mcmcout,
#                     samplls=fit1g@external$samplls, Fit=fit1g@Fit, 
#                     measure = "logl", thin = 5)
#temp1g


#blavCompare(fit, fit1g)
#blavCompare(fit,fit_null)

### for null model
#temp_null <- postpred_mgv(lavpartable=fit_null@ParTable, lavmodel=fit_null@Model, 
#                      lavoptions=fit_null@Options, lavsamplestats=fit_null@SampleStats, 
#                      lavdata=fit_null@Data, lavcache=fit_null@Cache, 
#                      lavjags=fit_null@external$mcmcout,
#                      samplls=fit_null@external$samplls, Fit=fit_null@Fit, 
#                      measure = "logl", thin = 5)
#temp_null



#ff <- BayesRelFit(obs=temp1g$chisqs[,1], rep=temp1g$chisqs[,2], nvar=ncol(fit1g@Data@X[[1]]), 
#            pD=temp1g$pd, N=fit1g@Data@nobs[[1]], ms = TRUE, Min1 = TRUE, Ngr = fit1g@Data@ngroups,
#            null_model=T, obs_null=temp_null$chisqs[,1], rep_null=temp_null$chisqs[,2], pD_null=temp_null$pd)
#ff


#t(apply(ff, 2, function(x){c(mean=mean(x),quantile(x,probs = c(.05,.5,.95)))}))



