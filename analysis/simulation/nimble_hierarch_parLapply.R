############################################
########### Hierarchical Power #############
############################################

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(spOccupancy)

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

high_hex_count <- emu_ratio %>% filter(occupancy == "high") %>% pull(hex_count) %>% sum()
#sample_sizes <- seq(log(100), log(3000), by = 0.3) %>% exp() %>% round()
sample_sizes <- seq(log(100), log(high_hex_count/2), by = 0.5) %>% exp() %>% round()

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
  filter(!is.na(n)) %>%
  group_by(emu, sim_num) %>%
  slice(rep(row_number(), length(sample_sizes))) %>% mutate(total_samp = sample_sizes) %>%
  group_by(emu) %>%
  mutate(sim_num = row_number()) %>%
  unite("sim_id", emu, sim_num)

# let's figure out relative area of each emu and total area
# emu_ratio %>%
#   filter(occupancy == "high") %>%
#   select(emu, hex_count) %>%
#   mutate(total_hex = sum(hex_count)) %>%
#   mutate(proportion = hex_count/total_hex)

# map high occupancy sims to their low occupancy counterpart
# set seed for random seed generator
set.seed(524879)
sim_map <- sim_scenarios_emu %>%
  #select(sim_id, psi, phi, p, total_samp) %>%
  filter(psi != 0.03) %>%
  separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
  rename(high_name = sim_id, high_hex_count = n) %>%
  left_join(sim_scenarios_emu %>%
              select(sim_id, psi, phi, p, total_samp,low_hex_count = n) %>%
              filter(psi == 0.03) %>%
              separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
              select(-psi) %>%
              rename(low_name = sim_id), by = c("emu", "total_samp", "phi", "p")) %>%
  group_by(phi, p, psi, total_samp) %>%
  mutate(scenario_id = cur_group_id(),
         seed = (1 + floor(runif(1) * 100000))) %>%
  ungroup() %>%
  mutate(high_n = round(total_samp*0.75), low_n = round(total_samp*0.25),
         low_psi = 0.03)

##############################################
########### Model and Power check ############
##############################################

source(here::here("R/sim_dataset.R"))

# function to initialize a model object for a given sample size
init_model <- function(n, year, visit, emu_vec, model_obj){

  nsite <- n
  nseason <- year
  nrep <- visit

  # make observed occurrence an array with site x year(season) x visit(rep)
  obs_occ_init <- array(sample(c(0,1), (nsite * nseason * nrep), replace = TRUE), c(nsite, nseason, nrep))

  # Build the model
  mod <- nimble::nimbleModel(
    code = model_obj,
    constants = list(nsite = n, nrep = nrep, nseason = nseason,
                     start_indexes = rep(1, nseason),
                     end_indexes = rep(nrep, nseason),
                     EMU = sample(emu_vec, nsite), num_EMU = n_distinct(emu_vec)
    ),
    data = list(y = obs_occ_init),
    inits = list(
      colonize = 0.5,
      init_occ = 0.5,
      detect = 0.5,
      persist_int = rep(0.5, (nseason-1)),
      ranef = rep(0, n_distinct(emu_vec)),
      sigma_ranef = 1
    )
  )

  # shouldn't NA, infinite, or positive (that's a dist issue)
  #Non-NA means we're fully initialized
  if (!is.finite(mod$calculate())){
    stop("Model did not initialize properly.")
  }

  # Build an MCMC
  conf <- nimble::configureMCMC(mod)
  conf$addMonitors(c("psi", "perc_change", "ranef"))
  mcmc <- nimble::buildMCMC(conf)

  # Compile
  complist <- nimble::compileNimble(mod, mcmc)

  return(complist)
}

sample_data <- function(data, sample_size){
  data_ids <- unique(data$site_id)
  samp_data <- data %>%
    filter(site_id %in% sample(data_ids, sample_size, replace = FALSE))
}


# let's start at the bottom
# sampled data needs to be all sampled data for that EMU and sample size (multiple scenarios)
fit_model_reps <- function(chunk, reps, n_year, n_visit){

  # get dataframe of scenarios
  map_df <- sim_map %>%
    dplyr::filter(scenario_id == chunk) %>%
    select(-scenario_id)

  sample_size <- unique(map_df$total_samp)
  high_n <- unique(map_df$high_n)
  low_n <- unique(map_df$low_n)

  #### Set up nimble model for that sample size ####

  # Model code for single EMU year estimate
  dynoccmod_code <- nimble::nimbleCode({

    # The whole likelihood for the dynamic occupancy model is contained inside
    # dDynOcc_sss. The suffix _sss indicates that persistence, colonization, and
    # detection are provided as scalars (one value for the whole site's
    # detection history). Other variants exist with suffixes like _svm (which
    # would mean that persistence is (s)calar, colonization is a (v)ector
    # varying with season, and detection is a (m)atrix varying with season and
    # with replicate)

    for (i in 1:nsite) {
      y[i, 1:nseason, 1:nrep] ~ dDynOcc_vss(probPersist = persist[i, 1:(nseason-1)],
                                            probColonize = colonize,
                                            init = init_occ,
                                            p = detect,
                                            start = start_indexes[1:nseason], # Start and end arguments allow you to provide ragged mtx data
                                            end = end_indexes[1:nseason])


    }

    # Define Priors
    for (i in 1:(nseason-1)){
      persist_int[i] ~ dunif(0,1)
      for (j in 1:nsite){
        logit(persist[j, i]) <- logit(persist_int[i]) + ranef[EMU[j]]
      }
    }

    for (r in 1:num_EMU) {
      # do sd = so life isn't ruined (might think its precision)
      ranef[r] ~ dnorm(0, sd = sigma_ranef)
    }

    colonize ~ dunif(0,1)
    init_occ ~ dunif(0,1)
    detect ~ dunif(0,1)
    sigma_ranef ~ dunif(0, 10)

    # Derive posterior for year
    psi[1] <-  init_occ
    for (i in 2:nseason){
      # gives the estimate for year based on mean persistance (not a level of random effect)
      psi[i] <- psi[i-1]*(persist_int[i-1]) + (1-psi[i-1])*colonize
    }
    perc_change <- (psi[10] - psi[1])/psi[1]
  })

  # create template model that can be updated with data
  compile_model <- init_model(n = sample_size, year = n_year, visit = n_visit, emu_vec = 1:n_distinct(map_df$emu), model_obj = dynoccmod_code)

  set.seed(unique(map_df$seed))
  # for reach replicate, generate high and low data, and fit model
  for (i in 1:reps){
    # get high data
    high_data <- purrr::pmap_dfr(map_df %>% select(-c(high_n, low_n)), function(high_name, emu, psi, phi, sd_phi, sd_gamma, p, high_hex_count, perc_red, total_samp, low_name,
                                 low_hex_count, low_psi, seed, year, visit){

      high_data <- sim_dataset(psi = psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                               n_sites = high_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
        append(c("sim_id" = high_name, "rep" = i))

      saveRDS(high_data, paste0(here::here("data/nimble/hier_simulated_data/"), "/", high_name, "_", i, "_simdata.rds"))

      return(high_data$obs_occ %>% mutate(emu = emu, landtype = "high") %>%
               # get unique site id across emu's
               mutate(site_id = paste0(emu, site_id)))
    }, year = n_year, visit = n_visit)

    # get low data
    low_data <- purrr::pmap_dfr(map_df %>% select(-c(high_n, low_n)), function(high_name, emu, psi, phi, sd_phi, sd_gamma, p, high_hex_count, perc_red, total_samp, low_name,
                                 low_hex_count, low_psi, seed, year, visit){

      low_data <- sim_dataset(psi = low_psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                              n_sites = low_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
        append(c("sim_id" = low_name, "rep" = i))

      saveRDS(low_data, paste0(here::here("data/nimble/emu_simulated_data/"), "/", low_name, "_", i, "_simdata.rds"))

      return(low_data$obs_occ %>% mutate(emu = emu, landtype = "low") %>%
               # get unique site id across emu's
               mutate(site_id = paste0(emu, site_id)))
    }, year = n_year, visit = n_visit)

    # sample observed occupancy
    occ_data <- bind_rows(sample_data(high_data, high_n),
                          sample_data(low_data, low_n)) %>%
      dplyr::arrange(visit)

    occ_array <- occ_data %>%
      dplyr::select(-c(site_id, emu, landtype)) %>%
      split(occ_data$visit) %>%
      purrr::map( ~ .x |> dplyr::select(-visit) |> as.matrix()) %>%
      simplify2array()

    emu_var <-  occ_data %>%
      filter(visit == 1) %>%
      mutate(emu_num = as.numeric(as.factor(emu))) %>%
      pull(emu_num)

    # list of new initialized variables
    new_inits <-  list(
      colonize = 0.5,
      init_occ = 0.5,
      detect = 0.5,
      persist_int = rep(0.5, (n_year-1)),
      ranef = rep(0, n_distinct(emu_vec)),
      sigma_ranef = 1)

    # update model with data
    compile_model$mod$y <- occ_array
    #compile_model$setData("y")
    fit <- nimble::runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                           samplesAsCodaMCMC = TRUE,
                           inits = new_inits)

    summary <- MCMCvis::MCMCsummary(fit, probs = c(0.025, 0.5, 0.95, 0.975)) |>
      dplyr::mutate(high_name = high_name, rep = i) |>
      tibble::rownames_to_column(var = "parameter")

    readr::write_csv(summary, paste0(here::here("data/nimble/emu_summary/"), "/", high_name, "_", i, "_summary.csv"))
    saveRDS(fit, paste0(here::here("data/nimble/emu_posterior/"), "/", high_name, "_", i, "_posterior.rds"))
  }
}

# initialize save out directories
dir.create(here::here("data/nimble/hier_summary"), recursive = TRUE)
dir.create(here::here("data/nimble/hier_posterior"))
dir.create(here::here("data/nimble/hier_simulated_data"))

ncores <- 50
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map', 'sim_dataset', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})

simn = 100

chunk_list <- unique(sim_map$scenario_id)
results <- parLapply(cl, chunk_list, fit_model_reps,
                     reps = simn, n_year = 10, n_visit = 2)

lapply(chunk_list, fit_model_reps,
       reps = 100, n_year = 10, n_visit = 2)
