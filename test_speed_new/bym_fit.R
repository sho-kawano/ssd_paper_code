library(LaplacesDemon); library(BayesLogit); library(Matrix); library(dplyr)

# this is a file with some utility functions I created 
source("~/coding/sparse_sae_code/utility_functions/mcmc_helper.R")

#' @description: runs one MCMC chain to fit a BYM random effects model
#' @param X the covariate matrix (from model.matrix, include intercept)
#' @param y the response 
#' @param d.var known design *variances*
#' @param A adjacency matrix 
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param hyp (Optional) list of the hyperparameters - depends on the prior choice 
#' @param save_latent (Optional) Save the MCMC samples of the latent variables (half-cauchy only)
#' @param ini (Optional) list of the initial values for any of the parameters 
#' @param verbose (Optional) Default T. Print statements with progress bar & computation time  
#' @returns list of MCMC containers (matrices) for each set of parameters. matrices compatible with bayesplot

bym_fit <- function(X, y, d.var, A, 
                    ndesired, nburn, nthin, 
                    hyp = list(), var.prior="inverse-gamma", save_latent=F, 
                    ini = list(), verbose=T
) {
  nsim <- nthin*ndesired
  X <- X[, -1, drop=F] # remove intercept
  n <- nrow(X); j <- ncol(X) # number of covariates (EXCLUDES intercept)
  
  # ---- set up containers ----
  Res_gamma <- mcmc_mat(nsim/nthin, (j+1)+2*n, 
                               names=c("intercept", colnames(X), paste0("v1_", 1:n), paste0("v2_", 1:n)))
  Res_beta <- mcmc_mat(nsim/nthin, j+1, c("intercept", colnames(X)))
  Res_v1 <- mcmc_mat(nsim/nthin, n, names=paste0("v1_", 1:n))
  Res_v2 <- mcmc_mat(nsim/nthin, n, names=paste0("v2_", 1:n))
  Res_sigma_1.sq <- mcmc_mat(nsim/nthin, 1, names="sigma.sq_v1")
  Res_sigma_2.sq <- mcmc_mat(nsim/nthin, 1, names="sigma.sq_v2")
  # ---- hyperparameters  ----
  # sigma.sq_vk ~ IG(c_k, d_k)
  c1 <- ifelse(is.null(hyp$c1), 5e-05, hyp$c1)
  d1 <- ifelse(is.null(hyp$d1), c1, hyp$d1)
  c2 <- ifelse(is.null(hyp$c2), 5e-05, hyp$c2)
  d2 <- ifelse(is.null(hyp$d2), c2, hyp$d2)
  # save prior details 
  details <- paste0("Priors: sigma_1.sq ~ IG(c1, d1): ", paste(c1, d1, sep=", "))
  details <- paste(details, paste0("sigma_2.sq ~ IG(c2, d2): ", 
                                   paste(c2, d2, sep=", ")), sep=" | ")
  
  # -----------------------------------------------------------------
  # ---- initial values: scale parameters ---- 
  ls_fit = lm(y~X)
  sigma_1.sq <- ifelse(is.null(ini$sigma_1.sq), 
                       mean(ls_fit$residuals^2), ini$sigma_1.sq)
  sigma_2.sq <- ifelse(is.null(ini$sigma_2.sq), 
                       mean(ls_fit$residuals^2), ini$sigma_2.sq)
  
  # ---- initial values: effects ---- 
  if(is.null(ini$beta)){
    beta <- ls_fit$coefficients[-1] # default - initial values 
    intr <- ls_fit$coefficients[1]
  }else{
    beta <- ini$beta[-1]
    intr <- ini$beta[1]
  }
  v1 <- ifelse(rep(is.null(ini$v1), n), rnorm(n, sd=sqrt(sigma_1.sq)/2), ini$v1)
  v2 <- ifelse(rep(is.null(ini$v2), n), rnorm(n, sd=sqrt(sigma_2.sq)/2), ini$v2)

  # blocked effects: theta = (1 Z)  %*% (intr gamma) = intr + X beta + v1 + v2 
  gamma <- c(beta, v1, v2)
  names(gamma) <- c(colnames(X), paste0("v1_", 1:n), paste0("v2_", 1:n))
  Z <- cbind(X,  diag(rep(1, n)), diag(rep(1, n))) 
  
  # ---- initial values: precision matrices---- 
  # precision matrix 
  ICAR <- diag(rowSums(A)) - A  # ICAR precision matrix 
  ICAR <- INLA:::inla.scale.model.bym(ICAR) # scaling
  Q_gam <- bdiag( diag(0, nrow=j), 
                  diag(1, nrow=n)/sigma_1.sq,  
                  ICAR/sigma_2.sq)
  
  # -----------------------------------------------------------------
  # ---- start MCMC chain ----  
  if(verbose){
    ptm <- start_chain(nsim, nthin, nburn)
    pb <- txtProgressBar(min=0, max=nsim+nburn, style=3)
  }
  for (index in 1:(nsim+nburn)) {
    # ---- 1A. update sigma.sq_v1 (IID)----
    v1 <- gamma[(1:n)+j]
    qform_v1 <- (t(v1)%*%v1) %>% as.numeric() #IID effects
    sigma_1.sq <- 1/rgamma(1, shape=n/2 + c1, rate=qform_v1/2  + d1)

    
    # ---- 1B. update sigma.sq_v2 (spatial) ----
    v2 <- gamma[(1:n)+n+j]
    qform_v2 <- (t(v2)%*%ICAR%*%v2) %>% as.numeric() #ICAR effects
    sigma_2.sq <- 1/rgamma(1, shape=n/2 + c2, rate=qform_v2/2  + d2)
    
    # ---- 1C. update precision matrix (depends on two sigmas) ----
    Q_gam <- bdiag( diag(0, nrow=j), 
                    diag(1, nrow=n)/sigma_1.sq,  
                    ICAR/sigma_2.sq )
    
    # ---- 2A. update gamma (block) ----
    m_gam <- t(Z/d.var)%*%y
    Zstd <- Z / sqrt(d.var)
    U_gam <- chol(forceSymmetric(t(Zstd)%*%Zstd + Q_gam))
    b <- rnorm(2*n+j)
    gamma <- backsolve(U_gam, backsolve(U_gam, m_gam, transpose=T) + b)
    
    # ---- 2B. solve for the intercept (after the sampling) ----
    v1 <- gamma[(1:n)+j]
    v2 <- gamma[(1:n)+n+j]
    intr <- mean(v1) + mean(v2)
    gamma[(1:n)+j] <- v1 - mean(v1)
    gamma[(1:n)+n+j] <- v2 - mean(v2)
    
    # ---- save results ----
    if (index > nburn && (index-nburn)%%nthin==0){
      # effects
      Res_gamma[(index-nburn)/nthin, ] <- c(intr, gamma)
      Res_beta[(index-nburn)/nthin,  ] <- c(intr, gamma[1:j])
      Res_v1[(index-nburn)/nthin,  ] <- gamma[(1:n)+j]
      Res_v2[(index-nburn)/nthin,  ] <- gamma[(1:n)+(n+j)]
      # scale parameters 
      Res_sigma_1.sq[(index-nburn)/nthin,  ] <- sigma_1.sq
      Res_sigma_2.sq[(index-nburn)/nthin,  ] <- sigma_2.sq
    }	
    if(verbose) setTxtProgressBar(pb, index)
  }
  # ---- End MCMC chain ----  
  if(verbose) print(proc.time() - ptm)
  # Save theta (small area means)
  Z <- cbind(intercept=1, X,  diag(rep(1, n)), diag(rep(1, n))) 
  Res_theta <- Res_gamma%*%t(Z)
  dimnames(Res_theta) <- list(NULL, paste0("theta_", 1:n))
  # ----- Save Results ---- 
  return(list(beta=Res_beta, v1=Res_v1, v2=Res_v2, 
    sigma_1.sq = Res_sigma_1.sq, sigma_2.sq = Res_sigma_2.sq,
    theta = Res_theta, prior_details=details)) 
}

bym_ver <- "v7. Uses backsolve. Edited Jul, '24 Wayyyyy faster."


# ---- TESTING ---- 
# library(tidyverse); library(bayesplot)
# load('test_data.RDA')
# 
# y = all_data$povPerc
# d = all_data$povPercSE^2
# X = model.matrix(~., all_data[, c("degree",  "assistance")])
# post.samples = bym_fit(X, y, d, A, ndesired=2000, nburn=1500, nthin=1)
# 
# post.samples$beta %>% mcmc_dens()
# post.samples$prior_details
