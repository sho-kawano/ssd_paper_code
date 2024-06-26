source("sim_analysis_functions/load_county.R")
key <- "9fb9cc8584a00b7006e7550c685129b65240b9aa" # PUT YOUR CENSUS API HERE
home_folder <- paste0(getwd(), "/data") 

# download NC only 
state <- "NC"
data <- load_county(state, home_folder, key, save_data = T)

# downloading the ACS data with geometry - used for plotting
state_abbrv="NC"
data_geom <- get_acs(geography = "county", variables = "B25071_001", 
               year = 2019, state = state_abbrv, geometry = T, cb=F)
saveRDS(data_geom, file="data/nc_data_w_geom.RDS")

# download south atlantic region 
region <- c("MD", "DE", "DC", "VA", "WV", "NC", "SC", "GA", "FL")
data <- load_county(region, home_folder, key, save_data = T)

data_geom <- get_acs(geography = "county", variables = "B25071_001", 
               year = 2019, state = region, geometry = T, cb=F)
saveRDS(data_geom, file="data/sa_data_w_geom.RDS")