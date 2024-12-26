library(tidyverse); library(parallel)
"
Prior to running this code:
1. please run `download_data.R` 
2. please enter the INPUTS below 

Note:
* this code uses FORK clusters which won't work on Windows
* random seed for each simulation (1-100) is simply set to (1-100)
"
# INPUTS 
key <- "9fb9cc8584a00b7006e7550c685129b65240b9aa" # your census API KEY
n.sims <- 20 # for the full simulation, set to 100 
n.cores <- 5 # change as needed 

# ----- Code to recreate South Atlantic Simulation Study ------ 
source("sim_analysis_functions/sim_setup.R")
source("sim_analysis_functions/run_sim.R")

load("data/data_multi-state/data.RDA")
all_data = all_data %>% mutate(state_code=strtrim(all_data$fips, 2)) %>% 
  left_join(fips_codes %>% 
              select(state, state_code, state_name) %>% 
              unique.data.frame(), by="state_code") 

data = list("SA", all_data, A)
covs = c("degree", "povPerc", "black", "state")

data %>% setup(home.folder=getwd(), y_name="rentBurden", 
               nsims=n.sims, transform = "log")

data %>% run_sim(home.folder=getwd(), y_name="rentBurden", cov_names = covs,
                 transform = "log", no.cores=n.cores, nsims=n.sims)

# ----- Post-Processing Simulation Results ------ 
source("sim_analysis_functions/calc_results.R")
cl <- makeForkCluster(n.cores)
doParallel::registerDoParallel(cl)

n <- nrow(all_data)
truth <- all_data$rentBurden
sim_folder = paste0(getwd(), "/sim_results_SA/")

# calculate Errors Across Simulations
errs_all_sims <- parLapply(cl, 1:n.sims, calc_err, sim.folder=sim_folder,
                           truth=truth, inv.transform=exp) %>% bind_rows()

# calculate CIs Across Simulations
ci_all_sims <- parLapply(cl, 1:n.sims, calc_CI, sim.folder=sim_folder,
                         theta=truth, inv.transform=exp) %>% bind_rows()

# calculate Posterior Means Across Simulations for bias calculations
postMean_all_sims <- parLapply(cl, 1:n.sims, calc_postMeans,
                               sim.folder=sim_folder, inv.transform=exp) %>%
  bind_rows()

# calculate Posterior Means of Incl. Prob
inclProb_all_sims <- parLapply(cl, 1:n.sims, calc_inclProb,
                               sim.folder=sim_folder ) %>% bind_rows()

#calculate Posterior Means of Random Effects
rEffects_all_sims <- parLapply(cl, 1:n.sims, calc_rEffects,
                               sim.folder=sim_folder ) %>% bind_rows()

# wrapping things up
stopCluster(cl)
save(postMean_all_sims, ci_all_sims, errs_all_sims, inclProb_all_sims,
     rEffects_all_sims, file=paste0(sim_folder, "/results_by_sim.RDA"))

