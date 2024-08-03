library(LaplacesDemon); library(dplyr)

source("~/coding/sparse_sae_code/utility_functions/mcmc_helper.R")

#' Description: runs one MCMC chain to fit a Datta-Mandal Model. 
#' @param X the covariates (from model.matrix, include intercept)
#' @param Y the response 
#' @param d.var known sample *variances*
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param hyp (Optional) dataframe of the hyperparameters. 
#' @param ini (Optional) list of the initial values for any of the parameters
#' @param verbose (Optional) Default T. Print statements with progress bar & computation time 
#' @returns list of MCMC containers (matrices) for each set of parameters and the details on fit

# NOTE: this sampler can be coded more efficiently 

dm_fit <- function(X, y, d.var, ndesired, nburn, nthin, 
                   hyp=list(), ini=list(), verbose=T) {
  nsim <- nthin*ndesired
  # default parameters/hyperparameters
  #if(length(hyp)<4) hyp = data.frame(a=0.5, b=0.5, c=1, d=2)
  if(length(hyp)<4) hyp = data.frame(a=3, b=2*mean(d.var), c=1, d=4)
  
  m <- nrow(X); j <- ncol(X)
  a <- hyp$a; b <- hyp$b; c <- hyp$c; d <- hyp$d
  details <- paste0("Priors: sigma.sq ~ IG(", a, ", ", b, ")")
  details <- paste(details, paste0("p ~ beta(", c, ", ", d, ")"), sep=" | ")
  
  # initial values
  ls_fit = lm(y~X-1)
  beta = ifelse(rep(is.null(ini$beta), j), ls_fit$coefficients, ini$beta)
  sigma.sq = ifelse(is.null(ini$sigma.sq), 
                    mean(ls_fit$residuals^2), ini$sigma.sq)
  ls_fit = NULL
  v = ifelse(rep(is.null(ini$v), m), rnorm(m, sd=sqrt(sigma.sq)), ini$v)
  delta <- ifelse(rep(is.null(ini$delta), m), rep(1, m), ini$delta)
  p <- ifelse(is.null(ini$p), 0.5, ini$p)
  
  # creates the containers to hold the results 
  # parameters 
  Res_beta = mcmc_mat(nsim/nthin, j, colnames(X))
  Res_v = mcmc_mat(nsim/nthin, m, names=paste0("v_", 1:m))
  Res_delta = mcmc_mat(nsim/nthin, m, names=paste0("delta_", 1:m))
  Res_p = mcmc_mat(nsim/nthin, 1, names="p")
  Res_p.delta = mcmc_mat(nsim/nthin, m, names=paste0("p.delta_", 1:m))
  Res_sigma.sq = mcmc_mat(nsim/nthin, 1, names="sigma sq.")
  # one for theta 
  Res_theta = mcmc_mat(nsim/nthin, m, names=paste0("theta_", 1:m) )
  
  # calculate to save computation time
  Xstd <- X/sqrt(d.var)
  cov <- solve(t(Xstd)%*%Xstd) 
  
  # MCMC chain 
  if(verbose){
    ptm <- start_chain(nsim, nthin, nburn)
    pb <- txtProgressBar(min=0, max=nsim+nburn, style=3)
  }
  
  for (index in 1:(nsim+nburn)) {
    # update sigma.sq
    sigma.sq <- 1/rgamma(1,shape=a+0.5*sum(delta),rate=b+0.5*sum(delta*v^2))
    
    # update p
    p <- rbeta(1, shape1=c+sum(delta), shape2=d+m-sum(delta))
    
    # update delta
    p.delta <- p/(p+(1-p)*sqrt((sigma.sq+d.var)/d.var)*exp(-0.5*(y-X%*%beta)^2*sigma.sq/(d.var+sigma.sq)/d.var))
    delta <- sapply(1:m, function(i){rbinom(1,1, p.delta[i])})
    
    # update beta
    # Xstd <- X/sqrt(d.var)
    # C_beta <- solve(t(Xstd)%*%Xstd) 
    mean <- apply(X*(y-delta*v)/d.var,2,sum)
    beta <- MASS::mvrnorm(1, mu=cov%*%mean,Sigma=cov)
    
    # update v
    mean.v <- sigma.sq*(y-X%*%beta)/(sigma.sq+d.var)
    var.v <- sigma.sq*d.var/(sigma.sq+d.var)
    v <- rnorm(m, mean=mean.v, sd=sqrt(var.v))*delta
    
    if (index > nburn && (index-nburn)%%nthin==0){
      # parameters 
      Res_beta[(index-nburn)/nthin, ] = beta
      Res_v[(index-nburn)/nthin,] = v
      Res_delta[(index-nburn)/nthin,  ] = delta
      Res_p[(index-nburn)/nthin, ] = p
      Res_p.delta[(index-nburn)/nthin,] = p.delta
      Res_sigma.sq[(index-nburn)/nthin,] = sigma.sq 
    }	
    if(verbose) setTxtProgressBar(pb, index)
  }
  # also save theta (a function of other params.)
  Res_theta <- Res_beta%*%t(X) + Res_delta*Res_v
  dimnames(Res_theta) <- list(NULL, paste0("theta_", 1:m))
  
  if(verbose) print(proc.time() - ptm)
  return(list(beta=Res_beta, v=Res_v, delta=Res_delta, p=Res_p,  
              sigma.sq=Res_sigma.sq, theta=Res_theta, p.delta=Res_p.delta, 
              prior_details=details))
}

dm_ver <- "v4: edited 2024-07-18. removed print code and added prior_details & working test code to chain"

# ---- TESTING ---- 
# library(tidyverse); library(bayesplot)
# load('test_data.RDA')
# 
# y = all_data$povPerc
# d = all_data$povPercSE^2
# X = model.matrix(~., all_data[, c("degree",  "assistance")])
# post.samples = dm_fit(X, y, d, ndesired=2000, nburn=10000, nthin=2)
# 
# post.samples$beta %>% mcmc_dens()
# post.samples$prior_details