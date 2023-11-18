
#######################################################################

# IMPORTANT PRE USED FUNCTIONS

#######################################################################


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