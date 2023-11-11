'
In Research Article, in order to calculate Bayesian fit indicies, we need to 
use the following function:
BayesRelFit , which use the output from the another function postpred_mgv

On of the important output is the calculated chisq.obs with the formula

      chisq.obs <- -2*(samplls[i, j, 1] -
                         samplls[i, j, 2])


samplls[i, j, 1] are the matrix of log likelihood values with the fllowing 
structure:

               1         2           ...           K of chains
1
2
3
4
...
N of datapoints


samplls[i, j, 2] is more problematic, It is calculated from:

https://rdrr.io/cran/blavaan/src/R/blavaan.R

      } else {
        samplls <- samp_lls(res, lavmcmc, lavobject = LAV, standata = rjarg$data)
      }

https://rdrr.io/cran/blavaan/src/R/blav_model_loglik.R

  } else {
    ## the model log-likelihoods have already been computed in stan
    llmat <- array(NA, c(nsamps, nchain, 2))
    lls <- loo::extract_log_lik(lavjags)
    llsat <- loo::extract_log_lik(lavjags, parameter_name = "log_lik_sat")
    for(j in 1:nchain){
      idx <- (j-1)*nsamps + itnums
      llmat[itnums,j,1] <- rowSums(lls[idx,])
      llmat[itnums,j,2] <- rowSums(llsat[idx,]) + llmat[itnums,j,1]
    }
    
    


https://rdrr.io/cran/blavaan/src/R/blavaan.R

if(jag.do.fit){
...
res <- try(do.call(rjcall, rjarg))


 rjcall <- blavmod$sample

'




'

 log_lik[mm] = wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sigma[mm]);
 if (do_test) {
   log_lik_sat[mm] = -log_lik[mm] + wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sstar[mm]);


 log_lik[mm] = wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sigma[mm]);
 if (do_test) {
   log_lik_sat[mm] = -log_lik[mm] + wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sstar[mm]);
'

log_lik[i] += normal_lpdf(x[i,j] | mu_x[i,j], epsilon[j]);
log_lik_sat[i] = -log_lik[i] + normal_lpdf(x[i,j] | mu_x[i,j], epsilon[j]);












