library(LaplacesDemon)
library(BayesLogit) 
library(Matrix)
library(dplyr)

# this is a file with some utility functions I created 
source("~/coding/spike_slab_fh/sparse_spatial_fh/mcmc_helper.R")

#' @description: runs one MCMC chain to fit a Spike & Slab BYM SAE model (no logit-intercept)
#' @param X the covariate matrix (from model.matrix, include intercept)
#' @param y the response 
#' @param D.i known design *variances*
#' @param A adjacency matrix 
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param ini (Optional) list of the initial values for any of the parameters
#' @param hyp (Optional) list of the hyperparameters - depends on the prior choice 
#' @param verbose (Optional) Default T. Print statements with prior information & iteration progress. 
#' @returns list of MCMC containers (matrices) for each set of parameters. matrices compatible with bayesplot

spslBYM_fit <- function(X, y, D.i, A, ndesired, nburn, nthin, 
                        hyp = list(), ini = list(),  verbose=T) {
  if(verbose){
    print("------------------------------------------------------")
    print("Fitting SSD random effects model.")
    print("Priors with the following Hyperparameters assumed: ")
  }
  nsim <- nthin*ndesired
  n <- nrow(X); j <- ncol(X) # number of covariates INCLUDES intercept 
  
  # ---- set up containers ----
  Res_theta <- mcmc_containers(nsim/nthin, n, names=paste0("theta_", 1:n))
  # effects 
  Res_beta <- mcmc_containers(nsim/nthin, j, colnames(X))
  Res_v1 <- mcmc_containers(nsim/nthin, n, names=paste0("v1_", 1:n))
  Res_v2 <- mcmc_containers(nsim/nthin, n, names=paste0("v2_", 1:n))
  Res_delta <- mcmc_containers(nsim/nthin, n, names=paste0("delta_", 1:n))
  # variance: random effects
  Res_sigma_1.sq <- mcmc_containers(nsim/nthin, 1, names="sigma_1.sq")
  Res_sigma_2.sq <- mcmc_containers(nsim/nthin, 1, names="sigma_2.sq")
  # inclusion probs (p.delta is the actual posterior incl. prob.)
  Res_p <- mcmc_containers(nsim/nthin, n, names=paste0("p_", 1:n))
  Res_p.delta <- mcmc_containers(nsim/nthin, n, names=paste0("p.delta_", 1:n)) 
  Res_pi <- mcmc_containers(nsim/nthin, 2*n, 
                            names=c(paste0("psi1_", 1:n), paste0("psi2_", 1:n)))
  # variance: incl. prob 
  Res_s1.sq <- mcmc_containers(nsim/nthin, 1, names="s1.sq")
  Res_s2.sq <- mcmc_containers(nsim/nthin, 1, names="s2.sq")
  
  # ---- hyperparameters  ----
  # uses default values unless specified via `hyp`
  # ---- variance of the regression coefficients 
  tau.sq <- ifelse(is.null(hyp$tau.sq), 100^2, hyp$tau.sq)
  
  # ---- random effect variance 
  # sigma.sq_vk ~ IG(c_k, d_k) for k=1,2
  c1 <- ifelse(is.null(hyp$c1), 5, hyp$c1)
  d1 <- ifelse(is.null(hyp$d1), c1, hyp$d1)
  c2 <- ifelse(is.null(hyp$c2), 5, hyp$c2)
  d2 <- ifelse(is.null(hyp$d2), c2, hyp$d2)
  
  # --- s.sq_pk ~ IG(a_k, b_k) for k=1,2
  a1 <- ifelse(is.null(hyp$a1), 5, hyp$a1)
  b1 <- ifelse(is.null(hyp$b1), 10, hyp$b1)
  a2 <- ifelse(is.null(hyp$a2), 5, hyp$a2)
  b2 <- ifelse(is.null(hyp$b2), 10, hyp$b2)

  # ----  mu (mean of logit random effects - psi's)
  mu <- ifelse(is.null(hyp$mu), 0, hyp$mu) # mean
  
  # ---- print hyperparameters 
  if(verbose){
    print(paste0("beta ~ N_j(0, tau.sq*I) w/ tau.sq: ", paste(tau.sq, sep=", ")))
    
    print(paste0("sigma_1.sq ~ IG(c1, d1): ", paste(c1, d1, sep=", ")))
    print(paste0("sigma_2.sq ~ IG(c2, d2): ", paste(c2, d2, sep=", ")))
    
    print(paste0("s1.sq ~ IG(a1, b1): ", paste(a1, b1, sep=", ")))
    print(paste0("s2.sq ~ IG(a2, b2): ", paste(a2, b2, sep=", ")))
    
    print(paste0("mu (mean of logit random effects): ", mu))
    print("------------------------------------------------------")
  }
  
  # ------------ initial values ---------------------------------
  # scale parameters 
  sigma_1.sq <- ifelse(is.null(ini$sigma_1.sq), 
                       0.5, ini$sigma_1.sq)
  sigma_2.sq <- ifelse(is.null(ini$sigma_2.sq), 
                       0.5, ini$sigma_2.sq)
  s1.sq <- ifelse(is.null(ini$s1.sq), 1, ini$s1.sq)
  s2.sq <- ifelse(is.null(ini$s2.sq), 1, ini$s2.sq)
  
  # precision matrices  
  ICAR <- diag(rowSums(A)) - A  # ICAR precision matrix 
  ICAR <- INLA:::inla.scale.model.bym(ICAR) # scaling
  # will be updated per iteration. dgCMatrix class 
  Q_gam <- bdiag( diag(1/tau.sq, nrow=j), 
                  diag(1, nrow=n)/sigma_1.sq, 
                  ICAR/sigma_2.sq) 
  Q_pi <- bdiag(diag(1, nrow=n)/s1.sq,  ICAR/s2.sq) 
  
  # effects
  beta = ifelse(rep(is.null(ini$beta), j),  lm(y~X-1)$coefficients, ini$beta)
  v1 <- ifelse(rep(is.null(ini$v1), n), rnorm(n, sd=sqrt(sigma_1.sq)/2), ini$v1)
  v2 <- ifelse(rep(is.null(ini$v2), n), rnorm(n, sd=sqrt(sigma_2.sq)/2), ini$v2)
  p <- ifelse(rep(is.null(ini$p), n), runif(n), ini$p)
  delta <- ifelse(rep(is.null(ini$delta), n), rbern(n, p), ini$delta)
  
  # blocked effects (note: theta = Z %*% gamma =  X beta + delta(v1 + v2) )
  gamma <- c(beta, v1, v2) %>% as.numeric()
  Z <- cbind(X, diag(delta), diag(delta)) 
  
  # inclusion probabilities 
  q_pi <- c(rep(mu, 2*n))   # prior mean of psi1 & psi2
  psi1 <- logit(p)   # iid effects
  psi2 <- logit(p) - psi1/2 # icar effects 
  pi <- c(psi1, psi2) 
  H <- cbind(diag(rep(1, n)), diag(rep(1, n))) # H does not change.
  
  # -----  latent variables 
  W <- diag(rpg(num=n, h=1, z=0 ))  #BayesLogit:rpg (Polyá-Gamma)  
  
  # ---- start MCMC chain ----  
  if(verbose) ptm <- start_chain(nsim, nthin, nburn)
  for (index in 1:(nsim+nburn)) {
    if(verbose & (index-nburn)%%1000==0 & index > nburn) print(index-nburn)
    
    # ---- 1A. update sigma_1.sq (IID)----
    v1 <- gamma[(1:n)+j]
    qform_v1 <- (t(v1)%*%v1) %>% as.numeric() #IID effects
    sigma_1.sq <- 1/rgamma(1, shape=n/2 + c1, rate=qform_v1/2  + d1)
    
    # ---- 1B. update sigma_2.sq (Spatial) ----
    v2 <- gamma[(1:n)+n+j]
    qform_v2 <- (t(v2)%*%ICAR%*%v2) %>% as.numeric() #ICAR effects
    sigma_2.sq <- 1/rgamma(1, shape=n/2 + c2, rate=qform_v2/2  + d2)
    
    # ---- 1C. update precision matrix (depends on two sigmas) ----
    Q_gam <- bdiag(diag(1/tau.sq, nrow=j), diag(1, nrow=n)/sigma_1.sq, ICAR/sigma_2.sq )
    
    # ---- 2. update gamma (blocked effects) ----
    # if none of the effects are included... 
    if(all(delta==0)){
      Q_beta <- bdiag(diag(1/tau.sq, nrow=j))
      m_beta <- t(X/D.i)%*%y # since v1=0, v2=0
      Xstd <- X / sqrt(D.i)
      Sig_beta <- solve(t(Xstd)%*%Xstd + Q_beta)
      beta <- MASS::mvrnorm(1, mu=Sig_beta%*%m_beta, Sigma=Sig_beta, tol=1e-200)
      gamma <- c(as.numeric(beta), rep(0, n), rep(0, n))
      
    }else{
      m_gam <- t(Z/D.i)%*%y
      Zstd <- Z / sqrt(D.i)
      Sig_gam <- solve(t(Zstd)%*%Zstd + Q_gam)
      # sampling all effects (adjustments for lack of symmetry) 
      gamma <- tryCatch(
        expr={
          MASS::mvrnorm(1, mu=Sig_gam%*%m_gam, Sigma=Sig_gam, tol=1e-200)
        },
        error = function(e){ # error in case Sigma is not symmetric
          Sig_gam <- Sig_gam %>% as.matrix() %>% as.symmetric.matrix()
          MASS::mvrnorm(1, mu=Sig_gam%*%m_gam, Sigma=Sig_gam, tol=1e-200)
        }
      )
      gamma <- gamma %>% as.numeric()
      # ! sum to ZERO constraint !
      v1 <- gamma[(1:n)+j]
      v2 <- gamma[(1:n)+n+j]
      gamma[1] <- gamma[1] + mean(v1) + mean(v2)
      gamma[(1:n)+j]   <- v1 - mean(v1)
      gamma[(1:n)+n+j] <- v2 - mean(v2)
    }
    # ---- 3A. update s1.sq (IID) ----
    psi1 <- pi[(1:n)] 
    qform_p1 <- (t(psi1)%*%(psi1)) %>% as.numeric() #IID effects
    s1.sq <- 1/rgamma(1,shape=n/2 + a1, rate=qform_p1/2 + b1)
    
    # ---- 3B. update s2.sq (spatial) ----
    psi2 <- pi[(1:n)+n]
    qform_p2 <- (t(psi2)%*%ICAR%*%(psi2)) %>% as.numeric() #ICAR effects
    s2.sq <- 1/rgamma(1,shape=n/2 + a2, rate=qform_p2/2 + b2)
    
    # ---- 3C. update the precision matrix (depends on s_p1 & s_p2) ----
    Q_pi <- bdiag(diag(1, nrow=n)/s1.sq,  ICAR/s2.sq)
    
    # ---- 4A. update delta  ----
    beta <- gamma[1:j]
    v1 <- gamma[(1:n)+j]
    v2 <- gamma[(1:n)+n+j]
    p_delta <- p / (p + (1-p)*exp( 0.5*(v1+v2)^2/D.i - (y-X%*%beta)*(v1+v2)/D.i))
    delta <- rbinom(n, size=1, p_delta)
    
    # ---- 4B. update blocked design matrix  ----
    Z <- X %>% cbind(diag(delta), diag(delta)) 
    
    # ---- 5A. update latent var W (polyá-gamma)----
    W <- diag(rpg(num=n, h=rep(1, n), z=(H%*%pi) %>% as.numeric())) 
    
    # ---- 5B. update pi (block) ----
    nu <- delta - 1/2
    m_pi <- t(H)%*%nu + Q_pi%*%q_pi 
    Sig_pi <- solve(t(H)%*%W%*%H + Q_pi)
    pi <- MASS::mvrnorm(1, mu=Sig_pi%*%m_pi, Sigma=Sig_pi, tol=1e-100)
    pi <- pi %>% as.numeric()
    p <- (H%*%pi)  %>% as.numeric() %>% invlogit()
    
    # ---- save results ----
    if (index > nburn && (index-nburn)%%nthin==0){
      # save theta (small area means)...... 
      Res_theta[(index-nburn)/nthin, ] <- (Z %*% gamma) %>% as.numeric()
      # effects
      Res_beta[(index-nburn)/nthin,  ] <- gamma[1:j]
      Res_v1[(index-nburn)/nthin,  ] <- gamma[(1:n)+j]
      Res_v2[(index-nburn)/nthin,  ] <- gamma[(1:n)+(n+j)]
      Res_delta[(index-nburn)/nthin,  ] <- delta 
      
      # scale: random effects
      Res_sigma_1.sq[(index-nburn)/nthin,  ] <- sigma_1.sq
      Res_sigma_2.sq[(index-nburn)/nthin,  ] <- sigma_2.sq
      
      # inclusion probs 
      Res_p[(index-nburn)/nthin,  ] <- p
      Res_p.delta[(index-nburn)/nthin,  ] <- p_delta
      Res_pi[(index-nburn)/nthin,  ] <- pi %>% as.numeric()
      
      # scale: incl. prob 
      Res_s1.sq[(index-nburn)/nthin,  ] <- s1.sq
      Res_s2.sq[(index-nburn)/nthin,  ] <- s2.sq
    }	
  }
  # ---- End MCMC chain ----
  if(verbose) print(proc.time() - ptm)
  result_list = list(theta=Res_theta, 
                     beta=Res_beta, v1=Res_v1, v2=Res_v2, delta=Res_delta, 
                     sigma_1.sq = Res_sigma_1.sq, sigma_2.sq = Res_sigma_2.sq,
                     p=Res_p, pi=Res_pi, p.delta=Res_p.delta,
                     s1.sq = Res_s1.sq, s2.sq=Res_s2.sq
  )
 return(result_list)
}

spslBYM_ver <- "v4.1 - simplified code. removed logit intercept. changed default priors. now designed to be used with scaled data."

