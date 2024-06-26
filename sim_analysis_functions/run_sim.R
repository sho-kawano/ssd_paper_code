library(parallel)
library(doParallel)

# runs a full simulation with all 4 models for tract-level data
run_sim <- function(load_data_list, home.folder, 
                        y_name, cov_names, transform, 
                        no.cores, nsims, 
                        bym_hyp=list(), new_hyp=list(),
                        sampler.dir="~/coding/spike_slab_fh/sparse_spatial_fh/"){
  # -------- Load Data ----------------
  selected_state = load_data_list[[1]][1]
  selected_county = load_data_list[[1]][2]
  # data needed to fit 
  all_data <- load_data_list[[2]]
  n <- nrow(all_data)
  theta <- all_data[,y_name] # this is treated as 'truth' in this simulation setup
  D <- all_data[,paste0(y_name, "SE")]^2
  X <- model.matrix(~., all_data[, cov_names, drop=F]) # design matrix
  A <- load_data_list[[3]]
  
  # -------- Directory Info ----------------
  if(is.null(selected_county) | is.na(selected_county) ){
    sim.folder = paste0(home.folder, "/sim_results_", selected_state, "/")
  }else{
    sim.folder = paste0(home.folder, "/sim_results_", 
                        selected_state, "_", selected_county, "/")
    sim.folder=str_replace_all(string=sim.folder, pattern=" ", repl="-")
  }
  # -------- Set up parallelized clusters ----------------
  options(mc.cores = no.cores)
  cl <- parallel::makeForkCluster(no.cores, setup_strategy = "sequential",
                                  outfile=paste0(sim.folder, "outfile_fullSim.txt"))
  doParallel::registerDoParallel(cl)
  
  # -------- Messages  ----------------
  print("--------------------------------------------------")
  print(paste0("Starting Sim. for ", selected_state, " ", selected_county))
  print(paste0("Number of areas: ", n))
  start_time <- Sys.time()
  print("------------------")
  print("Simulation Starting at: "); print(start_time)
  print("------------------"); print("Starting Timer...")
  timer <- proc.time()
  # -------- Start Parallelized Loop  ----------------
  foreach(sim=1:nsims) %dopar% {
    source(paste0(sampler.dir, "dm_fit.R"))
    source(paste0(sampler.dir, "fh_fit.R"))  
    source(paste0(sampler.dir, "bym_fit.R"))  
    source(paste0(sampler.dir, "spslBYM_fit.R"))  
    set.seed(sim)
    # ---- Load Simulated Direct Estimates ---- 
    subfolder = paste0(sim.folder, sprintf("%03d", sim))
    load(file=paste0(subfolder, "/y_gen.RDA")) #y_gen is loaded
    y_gen = as.numeric(y_gen) #this is the direct estimate
    
    # ---- Transform Y if desired ---- 
    if(transform=="log"){
      y.star <- log(y_gen)
      D.star <- D / y_gen^2
    } else if (transform=="logit"){
      y.star <- logit(y_gen)
      D.star <- D/(y_gen - y_gen^2)^2
    } else{
      y.star <- y_gen
      D.star <- D
    }
    # ---- Fit FH, DM, BYM ---- 
    fh_res <- fh_fit(X, y.star, D.star, 
                     ndesired=2000, nburn=9000, nthin=1, verbose=F)
    dm_res <- dm_fit(X, y.star, D.star, 
                     ndesired=2000, nburn=9000, nthin=1, verbose=F)
    bym_res <- bym_fit(X, y.star, D.star,  A, 
                       ndesired=2000, nburn=2000, nthin=1, 
                       hyp = bym_hyp, verbose=F)
    
    # ---- Fit new model ---- 
    # additional scaling step for new model 
    center <- mean(y.star)
    scale <- sd(y.star)
    y.scaled <- (y.star-center)/scale
    D.scaled <- D.star/scale^2
    
    new_res <- spslBYM_fit(X, y.scaled, D.scaled, A, 
                           ndesired=2000, nburn=2000, nthin=1,
                           hyp = new_hyp, verbose = F)
    
    # ---- Save Results & Wrap up  ----
    save(fh_res, file=paste0(subfolder, "/fh_res.RDA"))
    save(dm_res, file=paste0(subfolder, "/dm_res.RDA"))
    save(bym_res, file=paste0(subfolder, "/bym_res.RDA"))
    save(new_res, center, scale, file=paste0(subfolder, "/new_res.RDA"))
    print("________________________________________________")
    print(paste0("Simulation ", sim, " is finished."))
    print(Sys.time())
    print("________________________________________________")
  }
  duration <- proc.time() - timer; print(duration)
  
  # -------- Function to Calculate MSE / Coverage Results ----------------
  calc_err <- function(sim.no, sim.folder, truth){
    subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
    # load the results 
    load(file=paste0(subfolder, "/y_gen.RDA")) 
    load(file=paste0(subfolder, "/fh_res.RDA"))
    load(file=paste0(subfolder, "/dm_res.RDA"))
    load(file=paste0(subfolder, "/bym_res.RDA"))
    load(file=paste0(subfolder, "/new_res.RDA"))
    # the variables needed to compute mse 
    d_est = y_gen %>% as.numeric()
    fh_theta = fh_res$theta %>% inv.transform() %>% colMeans()
    dm_theta = dm_res$theta %>% inv.transform() %>% colMeans()
    bym_theta = bym_res$theta %>% inv.transform() %>% colMeans()
    new_theta = (new_res$theta*scale + center) %>% 
      inv.transform() %>% colMeans()
    # put results in df
    result = data.frame(area=1:length(truth),  sim=sim.no,
                        D.Est=truth - d_est,
                        FH=truth - fh_theta, 
                        DM=truth - dm_theta, 
                        BYM=truth - bym_theta,
                        New=truth - new_theta) %>%
      gather(key="model", value="err", -area, -sim)
    row.names(result) <- NULL
    return(result)
  }
  
  calc_CI <- function(sim.no, sim.folder, truth){
    subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
    load(file=paste0(subfolder, "/fh_res.RDA"))
    load(file=paste0(subfolder, "/dm_res.RDA"))
    load(file=paste0(subfolder, "/bym_res.RDA"))
    load(file=paste0(subfolder, "/new_res.RDA"))
    methods = c("FH", "DM", "BYM", "New")
    chains = list(fh_res$theta, dm_res$theta, bym_res$theta, (new_res$theta*scale + center))
    lapply(1:4,
           function(idx){
            ci_df = chains[[idx]] %>% inv.transform() %>%
              apply(2, function(x){quantile(x, c(0.05, 0.95))})
              rownames(ci_df) = c("lwr", "upr")
             ci_df = ci_df %>% t() %>% as.data.frame() %>%
               mutate(area=1:length(theta), sim=sim.no, model=methods[idx]) %>%
               mutate(truth=truth) %>%
               mutate(covered=ifelse(truth>lwr & truth < upr, T, F)) %>%
               select(sim, model, area, covered, truth, lwr, upr) %>%
               as.data.frame()
             rownames(ci_df) <- NULL
             return(ci_df)}
    ) %>% bind_rows()
  }
  # -------- Calculate Results ----------------
  if(transform=="log"){
    inv.transform = exp
  } else if (transform=="logit"){
    inv.transform = invlogit
  } else if (transform =="none") {
    inv.transform = function(x){x} # identity 
  }
  
  mse_df <- parLapply(cl, 1:nsims, calc_err, 
                             sim.folder=sim.folder, truth=theta) %>% 
    bind_rows() %>% 
    group_by(sim, model) %>% reframe(mse=mean(err^2)) %>%
    group_by(model) %>% reframe(tmse=sum(mse)) %>%
    arrange(desc(tmse)) 
  
  d_est.idx = which(mse_df$model=="D.Est")
  fh.idx = which(mse_df$model=="FH")
  
  mse_df = mse_df %>% 
    mutate(pdiff_dest = tmse/mse_df$tmse[d_est.idx] -1 ,
           pdiff_fh = tmse/mse_df$tmse[fh.idx] - 1) %>% 
    mutate(pdiff_dest = pdiff_dest %>% round(3), 
           pdiff_fh=pdiff_fh %>% round(3))
  
  covRate_df <- parLapply(cl, 1:nsims, calc_CI, 
                          sim.folder=sim.folder, truth=theta) %>%
    bind_rows() %>%
    group_by(sim, model) %>% reframe(cov_rate=mean(covered)) %>%
    group_by(model) %>% reframe(cov_rate=mean(cov_rate))
  
  # -------- Finish  ----------------
  stopCluster(cl)
  print("-------------------MSE TABLE ----------------------")
  print(mse_df)
  print("--------------------------------------------------")
  saveRDS(mse_df, file=paste0(sim.folder, "mse_result_full.RDS"))
  print("-------------------Cov.Rate TABLE ----------------------")
  print(covRate_df)
  print("--------------------------------------------------")
  saveRDS(covRate_df, file=paste0(sim.folder, "covRate.RDS"))
  print("----------------- Done! --------------------")
  print("--------------------------------------------------")
}

