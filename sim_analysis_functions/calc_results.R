library(dplyr)

calc_err <- function(sim.no, sim.folder, truth, inv.transform){
  subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
  # load the results 
  load(file=paste0(subfolder, "/y_gen.RDA")) 
  load(file=paste0(subfolder, "/fh_res.RDA"))
  load(file=paste0(subfolder, "/dm_res.RDA"))
  load(file=paste0(subfolder, "/car_res.RDA"))
  load(file=paste0(subfolder, "/bym_res.RDA"))
  load(file=paste0(subfolder, "/new_res.RDA"))
  
  # the variables needed to compute mse 
  d_est = y_gen %>% as.numeric()
  fh_theta = fh_res$theta %>% inv.transform() %>% colMeans()
  dm_theta = dm_res$theta %>% inv.transform() %>% colMeans()
  car_theta = car_res$theta %>% inv.transform() %>% colMeans()
  bym_theta = bym_res$theta %>% inv.transform() %>% colMeans()
  new_theta = (new_res$theta*new_res$scale + new_res$center) %>% 
    inv.transform() %>% colMeans()
  # put results in df
  result = data.frame(area=1:length(truth),  sim=sim.no,
                      D.Est=truth - d_est,
                      FH=truth - fh_theta, 
                      DM=truth - dm_theta,
                      CAR=truth - car_theta, 
                      BYM=truth - bym_theta,
                      New=truth - new_theta) %>%
    gather(key="model", value="err", -area, -sim)
  row.names(result) <- NULL
  return(result)
}
calc_CI <- function(sim.no, sim.folder, theta, inv.transform){
  subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
  load(file=paste0(subfolder, "/fh_res.RDA"))
  load(file=paste0(subfolder, "/dm_res.RDA"))
  load(file=paste0(subfolder, "/car_res.RDA"))
  load(file=paste0(subfolder, "/bym_res.RDA"))
  load(file=paste0(subfolder, "/new_res.RDA"))
  methods = c("FH", "DM", "CAR", "BYM", "New")
  chains = list(fh_res$theta, dm_res$theta, car_res$theta,
                bym_res$theta, (new_res$theta*new_res$scale + new_res$center) )
  lapply(1:length(chains),
         function(idx){
           ci_df = chains[[idx]] %>% inv.transform() %>%
             apply(2, function(x){quantile(x, c(0.05, 0.95))})
           rownames(ci_df) = c("lwr", "upr")
           ci_df = ci_df %>% t() %>% as.data.frame() %>%
             mutate(area=1:n(), sim=sim.no, model=methods[idx]) %>%
             mutate(truth=theta) %>%
             mutate(covered=ifelse(truth>lwr & truth < upr, T, F)) %>%
             select(sim, model, area, covered, truth, lwr, upr) %>%
             as.data.frame()
           rownames(ci_df) <- NULL
           return(ci_df)}
  ) %>% bind_rows()
}


calc_postMeans <- function(sim.no, sim.folder, inv.transform){
  subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
  
  load(file=paste0(subfolder, "/fh_res.RDA")) 
  load(file=paste0(subfolder, "/dm_res.RDA")) 
  load(file=paste0(subfolder, "/car_res.RDA"))
  load(file=paste0(subfolder, "/bym_res.RDA")) 
  load(file=paste0(subfolder, "/new_res.RDA")) 
  
  methods = c("FH", "DM", "CAR", "BYM", "New")
  chains = list(fh_res$theta, dm_res$theta, car_res$theta, 
                bym_res$theta, (new_res$theta*new_res$scale + new_res$center))
  n <- fh_res$theta %>% ncol()
  
  # add Direct Estimates 
  subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
  load(file=paste0(subfolder, "/y_gen.RDA")) 
  dEst_df <- data.frame(sim=sim.no, model="D.Est", area=1:n, 
                        post.mean =as.numeric(y_gen))
  # everything else 
  result <- lapply(1:length(chains), function(idx){
    df = chains[[idx]] %>% inv.transform() %>% colMeans() %>% 
      as.data.frame()
    colnames(df) <- c("post.mean")
    df = df %>% 
      as.data.frame() %>%
      mutate(sim=sim.no, area=1:n, model=methods[idx]) %>% 
      select(sim, model, area, post.mean)
    rownames(df) <- NULL
    return(df)}) %>% bind_rows()
  # attach the direct estimate & we're done
  result %>% rbind(dEst_df) 
}



calc_inclProb <- function(sim.no, sim.folder){
  subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
  load(file=paste0(subfolder, "/dm_res.RDA")) 
  load(file=paste0(subfolder, "/new_res.RDA")) 
  
  methods = c("DM", "New")
  chains = list(dm_res, new_res)
  
  lapply(1:2, function(idx){
    df = chains[[idx]]$p.delta %>% colMeans() %>% as.data.frame()
    colnames(df) <- c("post.mean")
    # result 
    df = df %>% 
      as.data.frame() %>%
      mutate(sim=sim.no, area=1:n(), model=methods[idx]) %>% 
      select(sim, model, area, post.mean)
    rownames(df) <- NULL
    return(df)}
  ) %>% bind_rows()
}


calc_rEffects <- function(sim.no, sim.folder){
  helper <- function(df,  method){
    colnames(df) <- c("post.mean")
    # result 
    df = df %>% 
      as.data.frame() %>%
      mutate(sim=sim.no, area=1:nrow(df), model=method) %>% 
      select(sim, model, area, post.mean)
    rownames(df) <- NULL
    return(df)
  }
  subfolder = paste0(sim.folder, sprintf("%03d", sim.no))
  df_list = list()
  
  load(file=paste0(subfolder, "/fh_res.RDA")) 
  df = fh_res$v %>% colMeans() %>% as.data.frame()
  df = helper(df, method="FH")
  df_list[[1]] <- df
  
  load(file=paste0(subfolder, "/dm_res.RDA")) 
  df = (dm_res$v*dm_res$delta) %>% colMeans() %>% as.data.frame()
  df = helper(df, method="DM")
  df_list[[2]] <- df
  
  load(file=paste0(subfolder, "/car_res.RDA")) 
  df = car_res$v %>% colMeans() %>% as.data.frame()
  df = helper(df, method="CAR")
  df_list[[3]] <- df
  
  load(file=paste0(subfolder, "/bym_res.RDA")) 
  df = (bym_res$v1 + bym_res$v2) %>% colMeans() %>% as.data.frame()
  df = helper(df, method="BYM")
  df_list[[4]] <- df
  
  load(file=paste0(subfolder, "/new_res.RDA")) 
  reffects_scl = new_res$delta*(new_res$v1 + new_res$v2)
  df = (reffects_scl*new_res$scale) %>% colMeans() %>% as.data.frame()
  df = helper(df, method="New")
  df_list[[5]] <- df
  
  df_list %>% bind_rows()
}


pairwiseMSE <- function(errs_all_sims, plot=F){
  mse_tbl = errs_all_sims %>% group_by(sim, model) %>% 
    reframe(mse=mean(err^2)) %>% 
    pivot_wider(names_from = model, values_from = mse) %>% 
    select(-sim) %>% 
    as.matrix()
  
  # result 6 x 6 table 
  result = matrix(NA, nrow=6, ncol=6)
  colnames(result) <- colnames( mse_tbl)
  if(plot){
    rownames(result) <- colnames( mse_tbl)
  }else{
    rownames(result) <- paste0(colnames( mse_tbl), " beats")
  }
  
  # create the rank matrix 
  rank_mat = mse_tbl %>% apply(1, rank) %>% t()
  for (i in 1:6){
    row.i = colSums(rank_mat[,-i] > rank_mat[,i]) / 100
    result[i, -i]  = row.i
  }
  if(plot){
    library(reshape2)
    result %>% melt() %>% 
      rename(`Model A`=Var1, `Model B`=Var2, `A bests B (%)`=value) %>% 
      ggplot(aes(`Model A`, `Model B`)) + 
      geom_tile(aes(fill = `A bests B (%)`)) + 
      geom_text(aes(label = round(`A bests B (%)`, 2))) +
      coord_flip()+ 
      scale_fill_gradient2() + 
      theme_bw()
  }else{
    return(result)
  }
}


pairwiseIntScore <- function(ci_all_sims, plot=F){
  iscore_tbl = ci_all_sims %>% 
    group_by(sim, model) %>%
    reframe(int_score = (upr-lwr) + (2/0.1)*(lwr-truth)*(truth<lwr) +
              (2/0.1)*(truth-upr)*(truth>upr)) %>% 
    group_by(sim, model) %>% 
    reframe(int_score = sum(int_score)) %>% 
    pivot_wider(names_from = model, values_from = int_score) %>% 
    select(-sim) %>% 
    as.matrix()
  
  # result 5 x 5 table (direct estimate not included)
  result = matrix(NA, nrow=5, ncol=5)
  colnames(result) <- colnames(iscore_tbl)
  if(plot){
    rownames(result) <- colnames(iscore_tbl)
  }else{
    rownames(result) <- paste0(colnames(iscore_tbl), " beats")
  }
  
  # create the rank matrix 
  rank_mat = iscore_tbl %>% apply(1, rank) %>% t()
  for (i in 1:5){
    row.i = colSums(rank_mat[,-i] > rank_mat[,i]) / 100
    result[i, -i]  = row.i
  }
  if(plot){
    library(reshape2)
    result %>% melt() %>% 
      rename(`Model A`=Var1, `Model B`=Var2, `A bests B (%)`=value) %>% 
      ggplot(aes(`Model A`, `Model B`)) + 
      geom_tile(aes(fill = `A bests B (%)`)) + 
      geom_text(aes(label = round(`A bests B (%)`, 2))) +
      coord_flip()+ 
      scale_fill_gradient2() + 
      theme_bw()
  }else{
    return(result)
  }
}

