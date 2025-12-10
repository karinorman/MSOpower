####################################################################################################################
## Power Analysis Goal: Owl occupancy rates must show a stable or increasing trend after 10 years of monitoring.
## The study design to verify this criterion must have a power of 90% (Type II error rate β = 0.10) to detect a
## 25% decline in occupancy rate over the 10-year period with a Type I error rate (α) of 0.10.
####################################################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(nimble)
library(nimbleEcology)
library(MCMCvis)

# real world vegtypes for each emu
emu_veg <- read.csv(here::here("data/EMU_veg_types.csv")) %>%
  # let's say which we think has high or low occupancy
  mutate(occupancy = case_when(
    veg_type_landfire == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "high",
    veg_type_landfire == "Madrean Upper Montane Conifer-Oak Forest and Woodland" ~ "high",
    veg_type_landfire == "Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Aspen Forest and Woodland" ~ "low",
    veg_type_landfire == "Inter-Mountain Basins Aspen-Mixed Conifer Forest and Woodland" ~ "low",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Savanna" ~ "low",
    veg_type_landfire == "Inter-Mountain Basins Subalpine Limber-Bristlecone Pine Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Bigtooth Maple Ravine Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Riparian Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Lodgepole Pine Forest" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Lodgepole Pine Forest" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "low",
    veg_type_landfire == "Southern Rocky Mountain Mesic Montane Mixed Conifer Forest and Woodland" ~ "high",
    veg_type_landfire == "Madrean Pinyon-Juniper Woodland" ~ "low",
    .default = NA
  ))

emu_ratio <- emu_veg %>%
  group_by(UNIT, occupancy) %>%
  summarize(hex_count = sum(hex_num)) %>%
  mutate(emu = case_when(
    UNIT == "Basin & Range - East" ~ "BRE",
    UNIT == "Basin & Range - West" ~ "BRW",
    UNIT == "Colorado Plateau" ~ "CP",
    UNIT == "Southern Rocky Mountains" ~ "SRM",
    UNIT == "Upper Gila Mountains" ~ "UGM"
  )) %>%
  ungroup()


###########################################
######## Define simulation parameters #####
###########################################

## Fixed study characteristics
nyear = 10
n_sites = sum(emu_veg$hex_num)
n_vis = 2

### Parameters ###

# get dataframe of all possible scenarios
sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = c(0.4, 0.8, .8)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  # and give each unique combination an ID
  mutate(sim_num = row_number(),
         occupancy = ifelse(psi == 0.03, "low", "high")) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  mutate(sd_phi = 0.04, sd_gamma = 0.01, perc_red = 0.25)

# get scenarios, one for each emu
sim_scenarios_emu <- bind_rows(sim_scenarios %>% mutate(emu = "BRE"),
                               sim_scenarios %>% mutate(emu = "BRW"),
                               sim_scenarios %>% mutate(emu = "CP"),
                               sim_scenarios %>% mutate(emu = "SRM"),
                               sim_scenarios %>% mutate(emu = "UGM")) %>%
  # get sample sizes for low and high occupancy for each emu
  left_join(emu_ratio %>%
              select(emu, occupancy, hex_count)) %>%
  # tidyr::pivot_wider(names_from = occupancy, values_from = hex_count) %>%
  # rename(low_n = low, high_n = high)) %>%
  # get the columns in the right order
  select(sim_num, emu, psi, phi, sd_phi, sd_gamma, p, n = hex_count, perc_red) %>%
  filter(!is.na(n))

# Let's get sample size of high quality hexes
# If an emu has enough area, we want the max sample size to be 2500, otherwise max sample is entire high quality area
sample_size_df <- emu_ratio %>%
  filter(occupancy == "high") %>%
  select(emu, hex_count) %>%
  mutate(log_max_samp = ifelse(hex_count > 2000, log(2000), log(hex_count))) %>%
  rowwise() %>%
  mutate(log_samp = list(c(seq(2.3, log_max_samp, by = 0.5), log_max_samp))) %>%
  unnest(log_samp) %>%
  mutate(samp_size = round(exp(log_samp)))

emu_sample_sizes <- emu_ratio %>%
  filter(occupancy == "high") %>%
  select(emu, hex_count) %>%
  left_join(sample_size_df %>% select(emu, n_samp = samp_size)) %>%
  group_by(emu, n_samp) %>%
  # expand again to get a row for each simulation scenario
  slice(rep(row_number() , n_distinct(sim_scenarios$sim_num))) %>%
  mutate(sim_num = 1:n_distinct(sim_scenarios$sim_num))

sim_scenarios_emu <- sim_scenarios_emu %>%
  left_join(emu_sample_sizes %>% select(sim_num, emu, n_samp)) %>%
  group_by(emu) %>%
  mutate(sim_num = row_number()) %>%
  unite("sim_id", emu, sim_num)


###########################################
########### Generate data sets ############
###########################################

source(here::here("R/sim_dataset.R"))
simn <- 100

#single_rep <- purrr::pmap(sim_scenarios_emu %>% select(-sim_id), sim_dataset, nyear = nyear, n_vis = 2) %>% set_names(sim_scenarios_emu$sim_id)

plan(multisession, workers = 70)
sim_list <- map(1:simn, ~furrr::future_pmap(sim_scenarios_emu %>%
                                              select(-sim_id, -n_samp), sim_dataset, nyear = nyear, n_vis = 2,
                                            .options=furrr_options(seed = TRUE)) %>%
                  set_names(sim_scenarios_emu$sim_id)
) %>%
  set_names(paste0("rep", 1:simn))

usethis::use_data(sim_list)

######### Format simulated data for model fitting ############

# reorder so top level of nested list is a sim scenario
sim_list_emu <- map(sim_scenarios_emu$sim_id, function(emu) {
  map(1:simn, ~pluck(sim_list, .x, emu)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios_emu$sim_id)

# map high occupancy sims to their low occupancy counterpart
sim_map_names <- sim_scenarios_emu %>%
  select(sim_id, psi, phi, p, n_samp) %>%
  filter(psi != 0.03) %>%
  separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
  select(-psi) %>%
  rename(high_name = sim_id) %>%
  left_join(sim_scenarios_emu %>%
              select(sim_id, psi, phi, p, n_samp) %>%
              filter(psi == 0.03) %>%
              separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
              select(-psi) %>%
              rename(low_name = sim_id)) %>%
  select(high_name, low_name, sample_size = n_samp)


# helper function to grab and sample observed occurrence for a rep
sample_data <- function(name, rep, sample_size){
  data <- pluck(sim_list_emu, name, rep, "obs_occ")
  data_ids <- unique(data$site_id)
  samp_data <- data %>%
    filter(site_id %in% sample(data_ids, sample_size, replace = FALSE)) %>%
    mutate(name = name, rep = rep)
}

# map low and high_name to the same scenario
sim_map <- sim_map_names %>%
  separate(high_name, c("emu", "scenario_num"), remove = FALSE) %>%
  select(-scenario_num) %>%
  group_by(high_name, sample_size) %>%
  nest() %>%
  mutate(rep = map(data, ~1:simn)) %>%
  unnest(cols = c("data", "rep")) %>%
  rename(high_n = sample_size) %>%
  mutate(low_n = round(high_n*(1/3)), total_n = (high_n + low_n)) %>%
  #mutate(temp_id = paste(high_name, rep, sep = "_")) %>%
  ungroup()

sim_map_data <- bind_rows(sim_map %>% rename(sample_size = high_n) %>%
                            mutate(name = high_name) %>%
                            select(-c(low_name, low_n)),
                          sim_map %>% rename(name = low_name, sample_size = low_n) %>%
                            select(-c(high_n))) %>%
  nest(data = c(name, rep, sample_size)) %>%
  rowwise() %>%
  mutate(samp_data = list(pmap(data, sample_data) %>% bind_rows())) %>%
  select(-data) %>%
  unnest(samp_data) %>%
  nest(data = -c(emu, total_n)) %>%
  ungroup()

sim_map_data_test <- bind_rows(sim_map %>% rename(sample_size = high_n) %>%
                                 mutate(name = high_name) %>%
                                 select(-c(low_name, low_n)),
                               sim_map %>% rename(name = low_name, sample_size = low_n) %>%
                                 select(-c(high_n))) %>%
  filter(rep < 11) %>%
  nest(data = c(name, rep, sample_size)) %>%
  rowwise() %>%
  mutate(samp_data = list(pmap(data, sample_data) %>% bind_rows())) %>%
  select(-data) %>%
  unnest(samp_data) %>%
  nest(data = -c(emu, total_n)) %>%
  ungroup() %>%
  filter(total_n < 15)

##############################################
########### Model and Power check ############
##############################################

# function to initialize a model object for a given sample size
init_model <- function(n, year, visit, model_obj){

  nsite <- n
  nseason <- year
  nrep <- visit
  
  print(c(nsite, "before model fit"))

  # make observed occurrence an array with site x year(season) x visit(rep)
  obs_occ_init <- array(sample(c(0,1), (nsite * nseason * nrep), replace = TRUE), c(nsite, nseason, nrep))

  # Build the model
  mod <- nimbleModel(
    code = model_obj,
    constants = list(nsite = n, nrep = nrep, nseason = nseason,
                     start_indexes = rep(1, nseason),
                     end_indexes = rep(nrep, nseason)
    ),
    data = list(y = obs_occ_init),
    inits = list(
      colonize = 0.5,
      init_occ = 0.5,
      detect = 0.5,
      persist = rep(0.5, (nseason-1))
    )
  )

  print(c(nsite, "after model fit"))
  
  # shouldn't NA, infinite, or positive (that's a dist issue)
  #Non-NA means we're fully initialized
  if (!is.finite(mod$calculate())){
    stop("Model did not initialize properly.")
  }

  # Build an MCMC
  conf <- configureMCMC(mod)
  conf$addMonitors(c("psi", "perc_change"))
  mcmc <- buildMCMC(conf)

  # Compile
  complist <- compileNimble(mod, mcmc)

  print(c(nsite, "after compile"))
  return(complist)
}


# let's start at the bottom
# sampled data needs to be all sampled data for that EMU and sample size (multiple scenarios)
#fit_model_reps <- function(emu, total_n, n_year, n_visit, data){
fit_model_reps <- function(df){  
  browser()
  print(emu)
  # Model code for single EMU year estimate
  dynoccmod_code <- nimbleCode({
    
    # The whole likelihood for the dynamic occupancy model is contained inside
    # dDynOcc_sss. The suffix _sss indicates that persistence, colonization, and
    # detection are provided as scalars (one value for the whole site's
    # detection history). Other variants exist with suffixes like _svm (which
    # would mean that persistence is (s)calar, colonization is a (v)ector
    # varying with season, and detection is a (m)atrix varying with season and
    # with replicate)
    
    for (i in 1:nsite) {
      y[i, 1:nseason, 1:nrep] ~ dDynOcc_vss(probPersist = persist[1:(nseason-1)],
                                            probColonize = colonize,
                                            init = init_occ,
                                            p = detect,
                                            start = start_indexes[1:nseason], # Start and end arguments allow you to provide ragged mtx data
                                            end = end_indexes[1:nseason])
      
      
    }
    
    # Define Priors
    for (i in 1:(nseason-1)){
      persist[i] ~ dunif(0,1)
    }
    
    colonize ~ dunif(0,1)
    init_occ ~ dunif(0,1)
    detect ~ dunif(0,1)
    
    # Derive posterior for year
    psi[1] <-  init_occ
    for (i in 2:nseason){
      psi[i] <- psi[i-1]*(persist[i-1]) + (1-psi[i-1])*colonize
    }
    perc_change <- (psi[10] - psi[1])/psi[1]
  })
  
  # create template model that can be updated with data
  compile_model <- init_model(n = total_n, year = n_year, visit = n_visit, model_obj = dynoccmod_code)

  # create dataframe with all scenarios and reps, maybe from sampled data?
  map_df <- data %>%
    select(scenario = high_name, rep_id = rep) %>%
    distinct()

  ## Map across scenarios for EMU and reps (multiple scenarios with the same sample size for each EMU)
  purrr::pmap(map_df, function(scenario, rep_id, data){

    # get only reps for the sample size that we have an initialized model for
    rep_data <- data %>%
      filter(high_name == scenario, rep == rep_id) %>%
      arrange(visit)

    obs_occ_array <- rep_data %>%
      dplyr::select(-c(high_name, site_id, name, rep)) %>%
      split(rep_data$visit) %>%
      purrr::map( ~ .x |> dplyr::select(-visit) |> as.matrix()) %>%
      simplify2array()

    # list of new initialized variables
    new_inits <-  list(
      colonize = 0.5,
      init_occ = 0.5,
      detect = 0.5,
      # number of years hard coded here!!
      persist = rep(0.5, 9))

    # update model with data
    compile_model$mod$y <- obs_occ_array
    #compile_model$setData("y")
    fit <- nimble::runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                           samplesAsCodaMCMC = TRUE,
                           inits = new_inits)

    summary <- MCMCvis::MCMCsummary(fit, probs = c(0.025, 0.5, 0.95, 0.975)) |>
      dplyr::mutate(high_name = scenario, rep = rep_id) |>
      tibble::rownames_to_column(var = "parameter")

    readr::write_csv(summary, paste0(here::here("data/nimble/emu_summary/"), "/", scenario, "_", rep_id, "_summary.csv"))
    saveRDS(fit, paste0(here::here("data/nimble/emu_posterior/"), "/", scenario, "_", rep_id, "_posterior.rds"))

  }, data = data)

}

plan(multisession, workers = 5)
furrr::future_pmap(sim_map_data_test, fit_model_reps, n_year = 10, n_visit = 2,
                   .options=furrr_options(seed = TRUE, packages = c("nimble"))
                   )


## Trying lapply

# test <- split(sim_map_data_test, 1:nrow(sim_map_data_test)) 
# 
# lapply(test, fit_model_reps)

###########################################################
################## Processing runs ########################
###########################################################

## This stuff we should already have, but we can read it in too
# ## Fixed study characteristics
# nyear = 10
# #n_sites = sum(emu_veg$hex_num)
# n_vis = 2
#
# # read in data we need
# true_occ_paired <- read.csv(here::here("data/true_occ_paired.csv"))
# testrds <- readRDS(here::here("data/nimble/emu_posterior/BRE_101_1_posterior.rds"))
# sim_map_metadata <- read.csv(here::here("data/sim_map_metadata.csv"))

# read in estimates
nimble_output <- purrr::map_dfr(list.files(here::here("data/nimble/emu_summary/"), full.names = TRUE), ~read.csv(.x) %>% mutate(parameter = colnames(testrds$chain1)))

# perc_change power checks
perc_change_check <- nimble_output %>%
  filter(parameter == "perc_change") %>%
  select(mean, ci025 = X2.5., ci97.5 = X97.5., high_name, rep) %>%
  left_join(true_occ_paired %>% select(high_name, true_perc_change = perc_change, rep)) %>%
  # pivot_longer(starts_with("ci"), names_to = "ci_type", values_to = "ci_value") %>%
  # mutate(ci_low = (mean - abs(ci_value)), ci_high = (mean + abs(ci_value))) %>%
  rowwise() %>%
  mutate(ci_two_tail = between(true_perc_change, ci025, ci97.5) & !between(0,  ci025, ci97.5)) %>%
  group_by(high_name) %>%
  summarize(ci_two_tail = sum(ci_two_tail)/n(),
            rep_count = n()) %>%
  left_join(sim_map_metadata %>% select(high_name, total_n) %>% distinct()) %>%
  left_join(sim_scenarios_emu %>% select(high_name = sim_id, psi, phi, p)) %>%
  group_by(psi, p, phi) %>%
  mutate(line_id = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

readr::write_csv(perc_change_check, here::here("data/nimble_power_check.csv"))

