source("sim_analysis_functions/load_county.R")
key <- "" # PUT YOUR CENSUS API HERE
home_folder <- paste0(getwd(), "/data") 

# download NC data used for analysis/simulation study
state <- "NC"
data <- load_county(state, home_folder, key, save_data = T)

# downloading the ACS data with geometry - used for plotting
state_abbrv="NC"
data_geom <- get_acs(geography = "county", variables = "B25071_001", 
               year = 2019, state = state_abbrv, geometry = T, cb=F)
saveRDS(data_geom, file="data/nc_data_w_geom.RDS")

# download south atlantic data used for analysis/simulation study
region <- c("MD", "DE", "DC", "VA", "WV", "NC", "SC", "GA", "FL")
data <- load_county(region, home_folder, key, save_data = T)

# downloading the ACS data with geometry - used for plotting
data_geom <- get_acs(geography = "county", variables = "B25071_001", 
               year = 2019, state = region, geometry = T, cb=F)
saveRDS(data_geom, file="data/sa_data_w_geom.RDS")