library(tidycensus)
library(LaplacesDemon)
library(tidyverse)

setup <- function(load_data_list, home.folder, y_name, transform, nsims){
  selected_state = load_data_list[[1]][1]
  selected_county = load_data_list[[1]][2]
  all_data = load_data_list[[2]]
  
  if(is.null(selected_county) | is.na(selected_county) ){
    sim.folder = paste0(home.folder, "/sim_results_", selected_state, "/")
  }else{
    sim.folder = paste0(home.folder, "/sim_results_", 
                        selected_state, "_", selected_county, "/")
    sim.folder=str_replace_all(string=sim.folder, pattern=" ", repl="-")
  }
  dir.create(sim.folder)
  curr.wd = getwd()
  setwd(home.folder)
  
  y = all_data[,y_name]
  D = all_data[,paste0(y_name, "SE")]^2
  
  for (sim in 1:nsims){
    newfolder = sprintf("%03d", sim)
    newpath = file.path(paste0(sim.folder, newfolder))
    dir.create(newpath)
    # generate data 
    # everything is transformed back to original units to avoid confusion
    set.seed(sim)
    if(transform=="log"){
      y_gen = rmvn(mu=log(y), Sigma=diag(D/y^2))
      y_gen = exp(y_gen)
    } else if (transform=="logit"){
      y_gen = rmvn(mu=logit(y), Sigma=diag(D / (y-y^2)^2))
      y_gen = invlogit(y_gen) 
    } else if (transform =="none") {
      y_gen = rmvn(mu=y, Sigma=diag(D))
    }
    # save results 
    folder = paste0(sim.folder, newfolder)
    save(y_gen, file=paste0(folder, "/y_gen.RDA"))
  }
  print(paste0("Set up is complete for ", selected_state))
  if(!selected_county %>% is.null()) print(selected_county)
  setwd(curr.wd)
}