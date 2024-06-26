library(LaplacesDemon); library(BayesLogit); library(Matrix); library(dplyr)

source("~/coding/spike_slab_fh/sparse_spatial_fh/mcmc_helper.R")

#' @description: runs one MCMC chain to fit a BYM random effects model
#' @param X the covariate matrix (from model.matrix, include intercept)
#' @param y the response 
#' @param D.i known sample *variances*
#' @param A adjacency matrix 
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param hyp (Optional) list of the hyperparameters 
#' @param ini (Optional) list of the initial values for any of the parameters 
#' @param verbose (Optional) Default T. Print statements with prior information & iteration progress. 
#' @returns list of MCMC containers (matrices) for each set of parameters. matrices compatible with bayesplot

bym_fit <- function(X, y, D.i, A, 
                    ndesired, nburn, nthin, hyp = list(), ini = list(), verbose=T
) {
  if(verbose){
    print("------------------------------------------------------")
    print("Fitting BYM random effects model")
    print("Priors with the following Hyperparameters assumed: ")
  }
  nsim <- nthin*ndesired
  X <- X[, -1, drop=F] # remove intercept
  n <- nrow(X); j <- ncol(X) # number of covariates (EXCLUDES intercept)
  
  # ---- set up containers ----
  Res_gamma <- mcmc_containers(nsim/nthin, (j+1)+2*n, 
                               names=c("intercept", colnames(X), paste0("v1_", 1:n), paste0("v2_", 1:n)))
  Res_beta <- mcmc_containers(nsim/nthin, j+1, c("intercept", colnames(X)))
  Res_v1 <- mcmc_containers(nsim/nthin, n, names=paste0("v1_", 1:n))
  Res_v2 <- mcmc_containers(nsim/nthin, n, names=paste0("v2_", 1:n))
  Res_sigma_1.sq <- mcmc_containers(nsim/nthin, 1, names="sigma.sq_v1")
  Res_sigma_2.sq <- mcmc_containers(nsim/nthin, 1, names="sigma.sq_v2")
  # ---- hyperparameters  ----
  # sigma.sq_vk ~ IG(c_k, d_k)
  c1 <- ifelse(is.null(hyp$c1), 5e-05, hyp$c1)
  d1 <- ifelse(is.null(hyp$d1), c1, hyp$d1)
  c2 <- ifelse(is.null(hyp$c2), 5e-05, hyp$c2)
  d2 <- ifelse(is.null(hyp$d2), c2, hyp$d2)
  # print out hyperparameters
  if(verbose) {
    print(paste0("sigma_1.sq ~ IG(c1, d1): ", paste(c1, d1, sep=", ")))
    print(paste0("sigma_2.sq ~ IG(c2, d2): ", paste(c2, d2, sep=", ")))
  }

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
  if(verbose) ptm <- start_chain(nsim, nthin, nburn)
  for (index in 1:(nsim+nburn)) {
    if((index-nburn)%%1000==0 & index > nburn & verbose){print(index-nburn)}
    
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
    m_gam <- t(Z/D.i)%*%y 
    Zstd <- Z / sqrt(D.i)
    Sig_gam <- solve(t(Zstd)%*%Zstd + Q_gam, tol=1e-20) %>% 
      as.matrix()  # it's not sparse so base matrix is fine
    if(!is.symmetric.matrix(Sig_gam)){
      Sig_gam <- as.symmetric.matrix(Sig_gam)
      if(!is.positive.definite(Sig_gam)){
        if(index > nburn) {
          warning("Numerical issue: prec.matrix not pos.def after burn-in")
        }
        Sig_gam <- as.positive.definite(Sig_gam)
      }
    }
    gamma <- MASS::mvrnorm(1, mu=Sig_gam%*%m_gam, Sigma=Sig_gam, tol=1e-100)
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
  }
  # ---- End MCMC chain ----  
  if(verbose) print(proc.time() - ptm)
  # Save theta (small area means)
  Z <- cbind(intercept=1, X,  diag(rep(1, n)), diag(rep(1, n))) 
  Res_theta <- Res_gamma%*%t(Z)
  dimnames(Res_theta) <- list(NULL, paste0("theta_", 1:n))
  # ----- Save Results ---- 
  return(list(
    #gamma=Res_gamma, 
    beta=Res_beta, v1=Res_v1, v2=Res_v2, 
    sigma_1.sq = Res_sigma_1.sq, sigma_2.sq = Res_sigma_2.sq,
    theta = Res_theta)) 
}

bym_ver <- "v6.1 - removed half-cauchy prior option"
