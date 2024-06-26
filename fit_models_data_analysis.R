source("samplers/fh_fit.R")
source("samplers/dm_fit.R")
source("samplers/bym_fit.R")
source("samplers/spslBYM_fit.R")

# ---- Loading Data ----
home_folder <- paste0(getwd(), "/data_analysis/")
load("data/data_multi-state/data.RDA")

# ---- Core Inputs ----
covs = c("degree", "assistance", "no_car", "povPerc",  
         "white", "black", "native", "asian", "hispanic")
X <- model.matrix(~., all_data[, covs, drop=F]) 

y <- all_data$rentBurden
D <- all_data$rentBurdenSE^2
y.star <- log(y)
D.star <- D / y^2

# ---- Run Chains ----
start_time <- Sys.time()
print("------------------")
print("Simulation Starting at: "); print(start_time)
timer <- proc.time()
# fit FH
fh_res <- fh_fit(X, y.star, D.star, 
                 ndesired=2500, nburn=10000, nthin=1, verbose=F)
duration <- proc.time() - timer; print(duration)
save(fh_res, file=paste0(home_folder, "fh_res.RDA"))
# fit DM
dm_res <- dm_fit(X, y.star, D.star, 
                 ndesired=2500, nburn=10000, nthin=1, verbose=F)
duration <- proc.time() - timer; print(duration)
save(dm_res, file=paste0(home_folder, "dm_res.RDA"))

# fit BYM
bym_res <- bym_fit(X, y.star, D.star,  A, 
                   ndesired=2500, nburn=1500, nthin=1, verbose=F)
duration <- proc.time() - timer; print(duration)
save(bym_res, file=paste0(home_folder, "bym_res.RDA"))

# fit SSD
center <- mean(y.star)
scale <- sd(y.star)
y.scaled <- (y.star-center)/scale
D.scaled <- D.star/scale^2

new_res <- spslBYM_fit(X, y.scaled, D.scaled, A, 
                       ndesired=2500, nburn=1500, nthin=1, verbose = F)
duration <- proc.time() - timer; print(duration)
save(new_res, center, scale, file=paste0(home_folder, "new_res.RDA"))


