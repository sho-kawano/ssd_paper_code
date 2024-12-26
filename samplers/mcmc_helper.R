library(LaplacesDemon); library(dplyr)

#' creates the "containers" to hold MCMC results. Compatible with bayesplot. 
#' @param iters number of desired iterations to store results for
#' @param param number of parameters
#' @param names of the parameters (required for bayesplot)
mcmc_mat <- function(iters, params, names){
  result = matrix(data=NA, nrow=iters, ncol=params)
  dimnames(result) <- list(NULL, names)
  names(dimnames(result)) <- c("iterations", "")
  return(result)
}

#' prints iterations and starts timer for MCMC chain 
#' @param nsim total number of simulations 
#' @param nthin amount of iterations to thin by 
#' @param nburn number of burn-in iterations
#' @returns returns a timer project
start_chain <- function(nsim, nthin, nburn){
  print(paste0("Running: ", nsim, 
               " simulations (excluding burn-in). Thinning by: ", nthin))
  print(paste0("Burn-in iterations: ", nburn ))
  return(proc.time())
}

#' Description: combine results from multiple chains into an array 
#' @param chains a list of posterior samples (a list of MCMC containers) from the same model 
#' @returns array with all of the parameters
combineChains <- function(chains){
  n.chain <- length(chains) 
  n.containers <-  length(chains[[1]])
  n.params <- sapply(1:n.containers, 
                     function(x){ncol(chains[[1]][[x]])}) %>% sum()
  n.iter <- chains[[1]][[1]] %>% nrow() 
  parnames <- sapply(1:n.containers, 
                     function(x){chains[[1]][[x]] %>% colnames()}) %>% unlist()
  
  # create array to store EVERYThing 
  all_chains <- array(NA, dim=c(n.iter, n.chain, n.params), 
                      dimnames = list(iterations=NULL, chains=1:n.chain, parameters=parnames))
  for(k in 1:n.chain){
    one_chain = array(unlist(chains[[k]]), dim=c(n.iter, 1, n.params), 
                      dimnames = list(iterations=NULL, chains=k, 
                                      parameters=parnames))
    all_chains[ , k, ] = one_chain
  }
  return(all_chains)
}

#' Description: transforms and changes the dimnames for a MCMC array 
#' @param samples MCMC samples (matrix or array)
#' @param par_names names of the new transformed parameters 
#' @param dim Optional (default 2) number of dimensions (2 for single chain, 3 for array with multiple chains)
#' @param transf the FUNCTION used to transform the parameters
#' @param nchains Optional (default 4).  Number of chains (only used if dim=3)
#' @returns array with the transformed parameters with dimnames relabeled 
transform_relabel <- function(samples, par_names, dim=2, nchains=4, transf=sqrt){
  samples <- samples %>% transf()
  if(dim==2){
    dimnames(samples) <- list(iterations=NULL, par_names)
  }else{
    dimnames(samples) <- list(iterations=NULL, chains=1:nchains, 
                              parameters=par_names)
  }
  return(samples)
} 



