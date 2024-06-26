library(Matrix);  library(LaplacesDemon); library(dplyr)

source("~/coding/spike_slab_fh/sparse_spatial_fh/mcmc_helper.R")

#' Description: runs one MCMC chain to fit a basic Fay-Herriot model (IID effects only)
#' @param X the covariates (from model.matrix, include intercept)
#' @param y the response 
#' @param D.i known sample *variances*
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param var.prior (Optional) Prior ("improper" or "inverse-gamma") 
#' @param hyp (Optional) list of the hyperparameters - depends on the prior choice 
#' @param ini (Optional) list of the initial values for any of the parameters
#' @param verbose (Optional) Default T. Print statements with prior information
#' @returns list of MCMC containers (matrices) for each set of parameters

fh_fit <- function(X, y, D.i, ndesired, nburn, nthin, var.prior="improper", 
                   hyp = list(), ini=list(), verbose=T){
  
  nsim <- nthin*ndesired
  # parameters
  m <- nrow(X); j <- ncol(X)
  # creates the containers to hold the results 
  Res_beta = mcmc_containers(nsim/nthin, j, colnames(X))
  Res_v = mcmc_containers(nsim/nthin, m, names=paste0("v_", 1:m))
  Res_sigma.sq = mcmc_containers(nsim/nthin, 1, names="sigma sq.")
  # one for theta (the target)
  Res_theta = mcmc_containers(nsim/nthin, m, names=paste0("theta_", 1:m) )
  
  # ----- hyperparameters -----
  if(var.prior=="inverse-gamma"){ 
    c <- ifelse(is.null(hyp$c), 0.5, hyp$c)
    d <- ifelse(is.null(hyp$d), 0.5, hyp$d)
    if(verbose){
      print(paste0("Hyperparams sigma.sq ~ IG(c, d): ", paste(c, d, sep=", "))) 
    }
  }else if (var.prior=="improper"){
    c <- -1; d <- 0
    if(verbose){
      print("Improper prior for beta & sigma.sq are assumed.")
    }
  }else{
    stop("Invalid variance prior. Implemented priors: inverse-gamma or improper")
  }  

  # ---- set initial values ----
  ls_fit = lm(y~X-1)
  beta = ifelse(rep(is.null(ini$beta), j), ls_fit$coefficients, ini$beta)
  sigma.sq = ifelse(is.null(ini$sigma.sq), 
                    mean(ls_fit$residuals^2), ini$sigma.sq)
  ls_fit = NULL
  v = ifelse(rep(is.null(ini$v), m), rnorm(m, sd=sqrt(sigma.sq)), ini$v)

  # MCMC chain 
  if(verbose) ptm <- start_chain(nsim, nthin, nburn)
  for (index in 1:(nsim+nburn)) {
    # update Beta
    Xstd <- X/sqrt(D.i)
    cov <- solve(t(Xstd)%*%Xstd) 
    mean <- apply(X*(y-v)/D.i,2,sum)
    mean <- cov%*%mean
    Beta <- MASS::mvrnorm(1, mu=mean,Sigma=cov)
    
    # update v
    mean.v <- sigma.sq*(y-X%*%Beta)/(sigma.sq+D.i)
    var.v <- sigma.sq*D.i/(sigma.sq+D.i)
    v <- rnorm(m, mean=mean.v, sd=sqrt(var.v))
    
    # update sigma.sq
    sigma.sq <- 1/rgamma(1,shape=0.5*m +c, rate=0.5*sum(v^2) + d)
    
    if (index > nburn && (index-nburn)%%nthin==0){
      Res_beta[(index-nburn)/nthin, ] = Beta
      Res_v[(index-nburn)/nthin,] = v
      Res_sigma.sq[(index-nburn)/nthin,] = sigma.sq
      # also save theta
      Res_theta[(index-nburn)/nthin, ] = X%*%Beta + v
    }	
  }
  if(verbose) print(proc.time() - ptm)
  return(list(beta=Res_beta, v=Res_v, sigma.sq=Res_sigma.sq, theta=Res_theta))
}
  

fh_ver <- "v4: edited October 19th.  Added initial values features & adjustments for modeling non-rate responses."

