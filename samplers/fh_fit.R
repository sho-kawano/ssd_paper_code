library(Matrix);  library(LaplacesDemon); library(dplyr)

source("~/coding/sparse_sae_code/utility_functions/mcmc_helper.R")

#' Description: runs one MCMC chain to fit a basic Fay-Herriot model (IID effects only)
#' @param X the covariates (from model.matrix, include intercept)
#' @param y the response 
#' @param d.var known sample *variances*
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param var.prior (Optional) Prior ("improper" or "inverse-gamma") 
#' @param hyp (Optional) list of the hyperparameters - depends on the prior choice 
#' @param ini (Optional) list of the initial values for any of the parameters
#' @param verbose (Optional) Default T. Print statements with progress bar & computation time 
#' @returns list of MCMC containers (matrices) for each set of parameters

# NOTE: this sampler can be coded more efficiently (but it's already fast)

fh_fit <- function(X, y, d.var, ndesired, nburn, nthin, var.prior="improper", 
                   hyp = list(), ini=list(), verbose=T){
  nsim <- nthin*ndesired
  m <- nrow(X)
  j <- ncol(X)
  
  # ----- matrices to hold results -----
  Res_beta <- mcmc_mat(nsim/nthin, j, colnames(X))
  Res_v <- mcmc_mat(nsim/nthin, m, names=paste0("v_", 1:m))
  Res_sigma.sq <- mcmc_mat(nsim/nthin, 1, names="sigma sq.")
  Res_theta <- mcmc_mat(nsim/nthin, m, names=paste0("theta_", 1:m))
  
  # ----- hyperparameters -----
  if(var.prior=="inverse-gamma"){ 
    c <- ifelse(is.null(hyp$c), 0.5, hyp$c)
    d <- ifelse(is.null(hyp$d), 0.5, hyp$d)
    details <- paste0("prior: sigma.sq ~ IG(c, d): ", paste(c, d, sep=", ")) 
  }else if (var.prior=="improper"){
    c <- -1; d <- 0
    details <- "Improper prior for beta & sigma.sq are assumed."
  }else{
    stop("Invalid variance prior. Implemented priors: inverse-gamma or improper")
  }  
  
  # ---- set initial values ----
  ls_fit <- lm(y~X-1) # use lm fit to set initial vals
  # if initial values are not specified, use the default values
  beta <- ifelse(rep(is.null(ini$beta), j), ls_fit$coefficients, ini$beta)
  sigma.sq <- ifelse(is.null(ini$sigma.sq), mean(ls_fit$residuals^2), ini$sigma.sq)
  v <- ifelse(rep(is.null(ini$v), m), rnorm(m, sd=sqrt(sigma.sq)), ini$v)
  
  if(verbose){
    ptm <- start_chain(nsim, nthin, nburn)
    pb <- txtProgressBar(min=0, max=nsim+nburn, style=3)
  }
  
  # Start MCMC chain 
  for (index in 1:(nsim+nburn)) {
    # update Beta
    Xstd <- X/sqrt(d.var)
    cov <- solve(t(Xstd)%*%Xstd) 
    mean <- apply(X*(y-v)/d.var,2,sum)
    mean <- cov%*%mean
    beta <- MASS::mvrnorm(1, mu=mean, Sigma=cov)
    
    # update v
    mean.v <- sigma.sq*(y-X%*%beta)/(sigma.sq+d.var)
    var.v <- sigma.sq*d.var/(sigma.sq+d.var)
    v <- rnorm(m, mean=mean.v, sd=sqrt(var.v))
    
    # update sigma.sq
    sigma.sq <- 1/rgamma(1,shape=0.5*m +c, rate=0.5*sum(v^2) + d)
    
    if (index > nburn && (index-nburn)%%nthin==0){
      Res_beta[(index-nburn)/nthin, ] = beta
      Res_v[(index-nburn)/nthin,] = v
      Res_sigma.sq[(index-nburn)/nthin,] = sigma.sq
      Res_theta[(index-nburn)/nthin, ] = X%*%beta + v
    }	
    if(verbose) setTxtProgressBar(pb, index)
  }
  if(verbose) {writeLines(""); print(proc.time() - ptm)}
  return(list(beta=Res_beta, v=Res_v, sigma.sq=Res_sigma.sq, theta=Res_theta, 
              prior_details=details))
}
  

fh_ver <- "v5: cleaned up version with working test code"

# ---- TESTING ---- 
# library(tidyverse); library(bayesplot)
# load('samplers/test_data.RDA')
# 
# y = all_data$povPerc
# d = all_data$povPercSE^2
# X = model.matrix(~., all_data[, c("degree",  "assistance")])
# post.samples = fh_fit(X, y, d, ndesired=2000, nburn=10000, nthin=2, verbose = F)
# 
# post.samples$beta %>% mcmc_dens()
# 
# eblup_fit <- sae::eblupFH(y~X-1, vardir=d)
# 
# #check against EBLUP
# plot(eblup_fit$eblup, colMeans(post.samples$theta),
#      xlab="EBLUP", ylab="Bayes", main="theta")
# 
# post.samples$prior_details
# 
