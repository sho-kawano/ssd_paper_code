library(LaplacesDemon); library(BayesLogit); library(Matrix); library(dplyr)

# this is a file with some utility functions I created 
source("samplers/mcmc_helper.R")

#' @description: runs one MCMC chain to fit a BYM random effects model
#' @param X the covariate matrix (from model.matrix, include intercept)
#' @param y the response 
#' @param d.var known design *variances*
#' @param A adjacency matrix 
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' 
#' @param hyp (Optional) list of the hyperparameters 
#' @param ini (Optional) list of the initial values for any of the parameters 
#' @param verbose (Optional) Default T. Print statements with progress bar & computation time  
#' @returns list of MCMC containers (matrices) for each set of parameters. matrices compatible with bayesplot

car_fit <- function(X, y, d.var, A, ndesired, nburn, nthin, 
                    hyp = list(), ini = list(), verbose=T
) {
  nsim <- nthin*ndesired
  n <- nrow(X); j <- ncol(X) # number of covariates (INCLUDES intercept)
  
  # ------------ Matrices to store results ------------ 
  Res_gamma <- mcmc_mat(nsim/nthin, j+n, names=c(colnames(X), paste0("v_", 1:n)))
  Res_beta <- mcmc_mat(nsim/nthin, j, c(colnames(X)))
  Res_v <- mcmc_mat(nsim/nthin, n, names=paste0("v_", 1:n))
  Res_sigma.sq <- mcmc_mat(nsim/nthin, 1, names="sigma.sq")
  Res_rho <- mcmc_mat(nsim/nthin, 1, names="rho")

  # ------------ Hyperparameters ------------ 
  # sigma.sq ~ IG(c_k, d_k)
  c <- ifelse(is.null(hyp$c), 5e-05, hyp$c)
  d <- ifelse(is.null(hyp$d), c, hyp$d)
  
  # metropolis-hastings tuning parameters 
  mh_sd <- ifelse(is.null(hyp$mh_sd), 0.25, hyp$mh_sd)
  
  # save prior details 
  details <- paste0("Priors: sigma.sq ~ IG(c, d): ", paste(c, d, sep=", "))
  details <- paste(details, paste0("MH Tuning SD's (rho): ", mh_sd), sep=" | ")
  
  # ------------ Initial Values ------------ 
  # ---- scale parameters ---- 
  ls_fit <- lm(y~X-1)
  beta <- ls_fit$coefficients  #L-S Estimate
  sigma.sq <- ifelse(is.null(ini$sigma.sq), 
                     mean(ls_fit$residuals^2), ini$sigma.sq)
  v <- rnorm(n, sd=sqrt(sigma.sq))
  #blocked effects variable
  gamma <- c(beta, v)
  Z <- X %>% cbind(diag(rep(1, n))) 
  rho <- 0.5
  
  # ---- precision matrices---- 
  Q_v <- diag(rowSums(A)) - rho * A
  Q_gam <- bdiag(diag(0, nrow=j), Q_v /sigma.sq) # improper prior for beta

  # ------------ MCMC chain ------------ 
  acs <- 0 # MH acceptance probability 
  if(verbose){
    ptm <- start_chain(nsim, nthin, nburn)
    pb <- txtProgressBar(min=0, max=nsim+nburn, style=3)
  }
  for (index in 1:(nsim+nburn)) {
    # ---- 1. update sigma.sq ----
    v <- gamma[(1:n)+j]
    qform_v <- (t(v)%*%Q_v%*%v) %>% as.numeric() 
    sigma.sq <- 1/rgamma(1, shape=n/2 + c, rate=qform_v/2  + d)
    
    # ---- 2. update rho (MH Step)----
    rho.prop <- rnorm(1, rho, mh_sd) # propose 
    rho.prop <- ifelse(rho.prop <= 1, rho.prop, 0.99) 
    rho.prop <- ifelse(rho.prop >= 0, rho.prop, 0.00) 
    
    # calculate acceptance probability
    l.prop <- dmvnp(v, 0, Omega=(diag(rowSums(A))  - rho.prop * A)/sigma.sq, log = TRUE)  + 
      dunif(rho.prop, log=TRUE)
    l.current <- dmvnp(v, 0, Omega=(diag(rowSums(A))  - rho * A)/sigma.sq, log = TRUE) + 
      dunif(rho,  log=TRUE)  
    log.q <- l.prop - l.current
    
    # accept/reject 
    if(log(runif(1)) < log.q) { 
      rho <- rho.prop
      acs <- acs + 1
    }
    
    # ---- update precision matrix ----
    Q_v <- diag(rowSums(A))  - rho * A
    Q_gam <- bdiag(diag(0, nrow=j), Q_v /sigma.sq)# improper prior for beta

    # ---- 3. update gamma (block) ----
    m_gam <- t(Z/d.var)%*%y
    Zstd <- Z / sqrt(d.var)
    U_gam <- chol(forceSymmetric(t(Zstd)%*%Zstd + Q_gam))
    b_gam <- rnorm(n+j)
    gamma <- backsolve(U_gam, backsolve(U_gam, m_gam, transpose=T) + b_gam)

    # ---- save results ----
    if (index > nburn && (index-nburn)%%nthin==0){
      Res_gamma[(index-nburn)/nthin,  ] <- gamma %>% as.numeric()
      Res_beta[(index-nburn)/nthin,  ] <- gamma[1:j]
      Res_v[(index-nburn)/nthin,  ] <- gamma[(1:n)+j]
      Res_sigma.sq[(index-nburn)/nthin,  ] <- sigma.sq
      Res_rho[(index-nburn)/nthin,] <- rho
    }	
    if(verbose) setTxtProgressBar(pb, index)
  }
  # ---- End MCMC chain ----  
  if(verbose){
    print(proc.time() - ptm)
    mh_prob_print <- paste("MH Acceptance Probabilities: ", 
                           paste0( (acs /(nsim+nburn)) %>% round(3) , 
                                   collapse =", "))
    details <- paste(details, mh_prob_print, sep=" | ")
  }
  # Save theta (small area means)
  Res_theta <- Res_gamma%*%t(Z)
  dimnames(Res_theta) <- list(NULL, paste0("theta_", 1:n))
  # ----- Save Results ---- 
  return(list(gamma=Res_gamma, beta=Res_beta, v=Res_v, 
              sigma.sq = Res_sigma.sq, rho = Res_rho, theta=Res_theta,
              prior_details=details
              )) 
}



car_ver <- "v2. Uses backsolve. Edited Aug,'24"


# ---- TESTING ---- 
# library(tidyverse); library(bayesplot)
# load('samplers/test_data.RDA')
# 
# y = all_data$povPerc
# d = all_data$povPercSE^2
# X = model.matrix(~., all_data[, c("degree",  "assistance")])
# post.samples = car_fit(X, y, d, A, ndesired=2000, nburn=1500, nthin=2)
# 
# post.samples$beta %>% mcmc_dens()
# post.samples$sigma.sq %>% mcmc_combo(nrow=2)
# post.samples$rho %>% mcmc_combo()
# post.samples$prior_details
