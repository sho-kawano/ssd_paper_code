library(LaplacesDemon); library(BayesLogit); library(Matrix); library(dplyr)

# this is a file with some utility functions I created 
source("utility_functions/mcmc_helper.R")

#' @description: runs one MCMC chain to fit a CAR random effects model
#' @param X the covariates (from model.matrix, include intercept)
#' @param y the response 
#' @param D.i known sample *variances*
#' @param A adjacency matrix 
#' 
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' 
#' @param hyp (Optional) list of the hyperparameters - depends on the prior choice 
#' @param var.prior (Optional) Prior of the *random effects variance*: "half-cauchy" (default) or "inverse-gamma"
#' @param constraint (Optional) whether linear constraints are used for the random effects (TRUE by default - recommended)
#' 
#' @returns list of MCMC containers (matrices) for each set of parameters

car_fit <- function(X, y, D.i, A, ndesired, nburn, nthin, hyp = list(), 
                         var.prior="half-cauchy", constraint=TRUE) {
  print("------------------------------------------------------")
  print("Fitting CAR random effects model")
  # ---- set up containers ----
  nsim <- nthin*ndesired
  n <- nrow(X); j <- ncol(X)
  
  print("------------------------------------------------------")
  # gamma
  Res_gamma <- mcmc_containers(nsim/nthin, j+n, 
                               names=c(colnames(X), paste0("v_", 1:n)))
  Res_beta <- mcmc_containers(nsim/nthin, j, colnames(X))
  Res_v <- mcmc_containers(nsim/nthin, n, names=paste0("v_", 1:n))
  Res_sig.sq_v <- mcmc_containers(nsim/nthin, 1, names="sigma.sq_v")
  Res_rho.v <- mcmc_containers(nsim/nthin, 1, names="rho_v")

  # ---- hyperparameters  ----
  # uses default values unless specified via `hyp`
  # hyper parameters for variances (mostly)
  if(var.prior=="inverse-gamma"){ 
    # c, d are inverse-gamma priors for sig.sq_v
    c <- ifelse(is.null(hyp$c), 0.5, hyp$c)
    d <- ifelse(is.null(hyp$d), 0.5, hyp$d)
    # print out hyperparameters
    print(paste0("Hyperparams sig.sq_v ~ IG(c, d): ", paste(c, d, sep=", ")))
  }else if (var.prior=="half-cauchy"){
    # scale params for half-cauchy
    c.v <- ifelse(is.null(hyp$c.v), 1, hyp$c.v)
    # print out hyperparameters
    print(paste0("Hyperparams: sigma_v ~ Half-Cauchy(c.v): ", paste(c.v, sep=", ")))
  }else{
    print("Invalid variance prior. Implemented: inverse-gamma or half-cauchy")
    return()
  }
  # metropolis-hastings tuning parameters 
  mh_sd.v <- ifelse(is.null(hyp$mh_sd.v), 0.5, hyp$mh_sd.v)
  print(paste0("MH Tuning SD's (rho.v): ", mh_sd.v))
  print("------------------------------------------------------")
  
  # ---- initial values ---- 
  beta <- lm(y~X-1)$coefficients  #L-S Estimate
  v <- rnorm(n, sd=0.01)
  #blocked effects variable
  gamma <- c(beta, v)
  Z <- X %>% cbind(diag(rep(1, n))) 
  sig.sq_v <- 0.01^2
  rho.v <- 0.5
  
  if(var.prior=="half-cauchy"){ u.v <- 100 }
  
  # ---- precision matrices---- 
  D.A <- diag(rowSums(A)) # row stochastic 
  Q_v <- D.A - rho.v * A
  Q_gam <- bdiag(diag(0, nrow=j), Q_v /sig.sq_v) # improper for beta

  # ---- start MCMC chain ----  
  ptm <- start_chain(nsim, nthin, nburn)
  acs <- 0 # acceptance probs 
  for (index in 1:(nsim+nburn)) {
    if((index-nburn)%%(500*nthin)==0 & index > nburn){print(index-nburn)}
    
    # ---- 1. update sigma.sq_v ----
    v <- gamma[(1:n)+j]
    qform_v <- (t(v)%*% Q_v %*%v) %>% as.numeric() 
    if(var.prior=="half-cauchy"){
      sig.sq_v <- 1/rgamma(n=1, shape=(n+1)/2, rate=qform_v/2  + 1/u.v)
      # ** latent: update u.v1 **
      u.v <- 1/rgamma(n=1, shape=1, rate=1/sig.sq_v + 1/c.v^2)
    } else{
      sig.sq_v <- 1/rgamma(1, shape=n/2 + c, rate=qform_v/2  + d)
    }  
    
    # ---- 2. update gamma (block) ----
    m_gam <- t(Z/D.i)%*%y 
    Sig_gam <- solve(t(Z)%*%diag(1/D.i)%*%Z + Q_gam)
    # turn into a PSD covariance matrix, reduce memory load
    # Sig_gam <- as(Sig_gam, "dppMatrix") %>% solve()  
    gamma <- MASS::mvrnorm(1, mu=Sig_gam%*%m_gam, Sigma=Sig_gam)
    # ! sum to zero constraint !
    v <- gamma[(1:n)+j]
    gamma[(1:n)+j] <- v - mean(v) 
    
    # ---- 3. update rho.v (MH Step)----
    rho.prop <- rnorm(1, rho.v, mh_sd.v) # propose 
    rho.prop <- ifelse(rho.prop <= 1, rho.prop, 0.99) 
    rho.prop <- ifelse(rho.prop >= 0, rho.prop, 0.01) 
    
    # calculate acceptance probability
    l.prop <- dmvnp(v, 0, Omega=(D.A - rho.prop * A)/sig.sq_v, log = TRUE)  + 
      dunif(rho.prop, log=TRUE)
    l.current <- dmvnp(v, 0, Omega=(D.A - rho.v * A)/sig.sq_v, log = TRUE) + 
      dunif(rho.v,  log=TRUE)  
    log.q <- l.prop - l.current
    
    # accept/reject 
    if(log(runif(1)) < log.q) { 
      rho.v <- rho.prop
      acs <- acs + 1
      # set new precision matrix
      Q_v <- D.A - rho.v * A
      Q_gam <- bdiag(diag(0, nrow=j), Q_v /sig.sq_v)
    }  
    
    # ---- save results ----
    if (index > nburn && (index-nburn)%%nthin==0){
      # gamma 
      Res_gamma[(index-nburn)/nthin,  ] <- gamma %>% as.numeric()
        Res_beta[(index-nburn)/nthin,  ] <- gamma[1:j]
        Res_v[(index-nburn)/nthin,  ] <- gamma[(1:n)+j]
      Res_sig.sq_v[(index-nburn)/nthin,  ] <- sig.sq_v
      Res_rho.v[(index-nburn)/nthin,] <- rho.v
    }	
  }
  # ---- End MCMC chain ----  
  print(proc.time() - ptm)
  # Save theta (small area means)
  Res_theta <- Res_gamma%*%t(Z)
  
  print("MH Acceptance Probabilities: ")
  print(paste0( (acs /(nsim+nburn)) %>% round(3) , collapse =", "))
  
  return(list(gamma=Res_gamma, beta=Res_beta, v=Res_v, 
              sig.sq_v = Res_sig.sq_v, rho.v = Res_rho.v, theta=Res_theta))
}


CAR_ver <- "v2: untested.. needs some valiation work. could also be sped up and use backsolve"


# ----------------------------- TEST -----------------------------
# 
# simulation_name <- "Illinois"
# setwd(paste0("~/coding/spike_slab_fh/full_simulation_study/", simulation_name, "/load_data"))  
# load('acs_data.RDA')
# D.i <- all_data$povPercSE^2
# covariates <- all_data[, 8:ncol(all_data)] ### CHECK
# X = model.matrix(~., covariates)
# y <- all_data$povPerc; n <- length(y)
# load('adj_mat.RDS')
# 
# test_chain = newCAR_fit(X, y, D.i, A, ndesired=3000, nburn=1000, nthin=1)
# chosen_model = test_chain
# 
# 
# ##  Random Effects Variance 
# plt1 <- chosen_model[[7]] %>%  mcmc_dens()
# plt2 <- chosen_model[[7]] %>%  mcmc_combo(widths=c(2, 1), c("trace", "acf"))
# ggpubr::ggarrange(plt1, plt2, nrow=2)
# 
# 
# ##  Incl. Probs Variance 
# plt1 <- chosen_model[[8]] %>%  mcmc_dens()
# plt2 <- chosen_model[[8]] %>%  mcmc_combo(widths=c(2, 1), c("trace", "acf"))
# ggpubr::ggarrange(plt1, plt2, nrow=2)
# 
# ##  Rho parameter for v
# plt1 <- chosen_model[[9]] %>%  mcmc_dens()
# plt2 <- chosen_model[[9]] %>%  mcmc_combo(widths=c(2, 1), c("trace", "acf"))
# ggpubr::ggarrange(plt1, plt2, nrow=2)
# 
# ##  Rho parameter for p
# plt1 <- chosen_model[[10]] %>%  mcmc_dens()
# plt2 <- chosen_model[[10]] %>%  mcmc_combo(widths=c(2, 1), c("trace", "acf"))
# ggpubr::ggarrange(plt1, plt2, nrow=2)
# 
