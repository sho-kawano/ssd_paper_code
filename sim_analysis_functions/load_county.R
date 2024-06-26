library(tidycensus)
library(tigris)
library(spdep)
library(spatialreg)
library(tidyverse)

# v17 <- load_variables(2019, "acs5", cache = TRUE)
# View(v17)

# downloads ACS / Adjancency Matrix
load_county <- function(selected_state, home.folder, api.key, 
                       selected_yr=2019, save_data=F, geom=F){
  
  census_api_key(api.key)
  
  # --------------- LOAD POPULATION DATA --------------- 
  pop <- get_acs(geography = "county", 
                 variables = "B01003_001", 
                 year = selected_yr,
                 state = selected_state, geometry = geom)
  
  if(!geom){
    # save as dataframe
    main_df <- data.frame(fips=pop$GEOID, geo_name=pop$NAME, pop = pop$estimate)
  }else{
    # save as sf/dataframe 
    main_df <- pop %>% rename(fips=GEOID, geo_name=NAME, pop=estimate) %>% 
      select(fips, geo_name, pop, geometry)
  }
  
  # --------------- LOAD INCOME DATA --------------- 
  income <- get_acs(geography = "county", 
                    variables = "B19013_001", 
                    year = selected_yr,
                    state = selected_state)
  income <- data.frame(fips=income$GEOID, 
                       income=income$estimate, incomeSE=income$moe/1.645)  
  
  # --------------- LOAD POVERTY DATA ---------------  
  pov <- get_acs(geography = "county", 
                 variables = "B17001_002",
                 year = selected_yr,
                 state = selected_state)
  
  poverty <- data.frame(fips=pop$GEOID, pov=pov$estimate, povSE=pov$moe/1.645) %>%
    mutate(povPerc=pov/pop$estimate, povPercSE=povSE/pop$estimate)
  
  # --------------- LOAD MEDIAN RENT AS A %age of Income ---------------  
  rent <- get_acs(geography = "county", variables = "B25071_001",
                  year = selected_yr,
                  state = selected_state)
  
  rent <- data.frame(fips=pop$GEOID,
                     rentBurden=rent$estimate, rentBurdenSE=rent$moe/1.645) %>%
    mutate(rentBurden=rentBurden/100, rentBurdenSE=rentBurdenSE/100) # decimal %age
  
  
  # --------------- LOAD EDUCATION DATA --------------- 
  college = get_acs(geography = "county",
                    variables = c(degree = "B06009_005", total = "B06009_001"),
                    state = selected_state, 
                    year = selected_yr) %>%
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>%
    mutate(degree=degree/total) %>%
    select(-total)
  
  # --------------- LOAD RACE DATA ---------------  
  race = get_acs(geography = "county", 
                 variables = c(total="B02001_001", white="B02001_002", 
                               black="B02001_003", native="B02001_004", 
                               asian="B02001_005", 
                               other1="B02001_006", other2="B02001_007", 
                               mult1="B02001_008", mult2="B02001_009", 
                               mult3="B02001_010"), 
                 state = selected_state, 
                 year = selected_yr)
  
  race = race %>% 
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>% 
    mutate(other=other1+other2+mult1+mult2+mult3) %>% 
    select(-other1, -other2, -mult1, -mult2, -mult3) %>% 
    relocate(geo_name, total, white, native, black, asian, other) %>% 
    mutate(white=white/total, native=native/total, black=black/total, 
           asian=asian/total, other=other/total)
  race = race %>% select(geo_name, white, native, black, asian)
  
  # --------------- LOAD LATINO DATA --------------- 
  latino = get_acs(geography = "county", 
                   variables = c(non_hispanic="B08134_001", hispanic="B08134_010"), 
                   state = selected_state, 
                   year = selected_yr) %>% 
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>% 
    mutate(total = non_hispanic+hispanic) %>% 
    mutate(non_hispanic=non_hispanic/total, hispanic= hispanic/total) %>% 
    select(-total, -non_hispanic)
  
  # ---------------  ASSISTANCE & FOOD STAMPS ---------------  
  #PUBLIC ASSISTANCE INCOME OR FOOD STAMPS/SNAP IN THE PAST 12 MONTHS FOR HOUSEHOLDS
  assistance = get_acs(geography = "county",
                       variables = c(assistance = "B19058_002", total = "B19058_001"),
                       state = selected_state, 
                       year = selected_yr) %>%
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>%
    mutate(assistance=assistance/total) %>%
    select(-total)
  
  # ---------------  Households with Children/Marital Status 
  parents = get_acs(geography = "county", 
                    variables = c(married = "B09002_002", 
                                  male_only = "B09002_009", 
                                  female_only = "B09002_015"), 
                    state = selected_state, 
                    year = selected_yr) %>%
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>%
    mutate(total=female_only + male_only + married) %>% 
    mutate(married=married/total, 
           male_only = male_only/total, 
           female_only=female_only/total) %>% 
    select(-total)
  
  # ---------------  Vehicle Data  
  vehicle = get_acs(geography = "county",
                    variables = c(no_car = "B08014_002", total = "B08014_001"),
                    state = selected_state,
                    year = selected_yr)  %>%
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>%
    mutate(no_car=no_car/total) %>%
    select(-total)
  
  # ---------------  Load Age Data 
  age = get_acs(geography = "county", 
                variables = c(total="B06001_001",
                              age.under_5="B06001_002", age.05_17 ="B06001_003", 
                              age.18_24 ="B06001_004", age.25_34 ="B06001_005", 
                              age.35_44 ="B06001_006", age.45_54 ="B06001_007", 
                              age.55_59 ="B06001_008", age.60_61 ="B06001_009", 
                              age.62_64 ="B06001_010", age.65_74 ="B06001_011",
                              age.75_over ="B06001_012"), 
                state = selected_state, 
                year = selected_yr) %>% 
    rename(geo_name = NAME) %>% 
    select(geo_name, variable, estimate) %>%
    spread(key=variable, value=estimate ) %>%
    mutate(age.60_plus = age.60_61+age.62_64+age.65_74+age.75_over) %>% 
    # leave out 44-59 age group. 
    select(-age.45_54, -age.55_59, -age.60_61, -age.62_64, -age.65_74, -age.75_over) %>% 
    mutate_at(vars(-geo_name, -total), funs(. / total)) %>%
    relocate(age.under_5, age.05_17)%>%
    select(-total)
  
  # ---------------  Combine
  all_data = main_df %>% 
    left_join(income, by = join_by(fips)) %>% 
    left_join(poverty, by = join_by(fips)) %>% 
    left_join(rent, by = join_by(fips)) %>%
    left_join(college, by = join_by(geo_name)) %>% 
    left_join(race, by = join_by(geo_name)) %>% 
    left_join(latino, by = join_by(geo_name)) %>% 
    left_join(assistance, by = join_by(geo_name)) %>% 
    left_join(parents, by = join_by(geo_name)) %>%
    left_join(vehicle, by = join_by(geo_name)) %>%
    left_join(age, by = join_by(geo_name)) %>%
    arrange(fips)
  

  # ---------------  Adjacency Matrix  ---------------
  acs_shape <- counties(cb=F, state=selected_state, year=selected_yr) 
  rownames(acs_shape) <- acs_shape$GEOID
  acs_shape <- acs_shape %>% arrange(GEOID)
  A <- nb2mat(poly2nb(acs_shape), style='B', zero.policy = F)
  
  # --------------- SAVE DATA (if necessary) ---------------
  if(save_data){
    if(length(selected_state)==1){
      save.folder = paste0(home.folder, "/data_", selected_state, "/")
    }else{
      save.folder = paste0(home.folder, "/data_multi-state/")
    }
    dir.create(save.folder)
    save(all_data, A, file=paste0(save.folder, "data.RDA"))
  }
  return(list(selected_state, all_data, A))
}
