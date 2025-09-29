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
library(parallel)

## Fixed study characteristics
nyear = 10
n_sites = sum(emu_veg$hex_num)
n_vis = 2

load(here::here("data/sim_map_names.rda"))

##############################################
########### Model and Power check ############
##############################################

source(here::here("R/sim_dataset.R"))

# function to initialize a model object for a given sample size
init_model <- function(n, year, visit, model_obj){

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
                     end_indexes = rep(nrep, nseason)
    ),
    data = list(y = obs_occ_init),
    inits = list(
      colonize = 0.5,
      init_occ = 0.5,
      detect = 0.5,
      persist_intercept = rep(0.5, (nseason-1)),
      beta = 0.5
    )
  )

  # shouldn't NA, infinite, or positive (that's a dist issue)
  #Non-NA means we're fully initialized
  if (!is.finite(mod$calculate())){
    stop("Model did not initialize properly.")
  }

  # Build an MCMC
  conf <- nimble::configureMCMC(mod)
  conf$addMonitors(c("psi", "perc_change"))
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
fit_model_reps <- function(chunk, reps, n_year, n_visit, data){

  # get dataframe of scenarios
  map_df <- data %>%
    dplyr::filter(chunk_num == chunk) %>%
    select(-chunk_num)

  sample_size <- unique(map_df$total_n)

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
      y[i, 1:nseason, 1:nrep] ~ dDynOcc_vss(probPersist = persist[1:(nseason-1)],
                                            probColonize = colonize,
                                            init = init_occ,
                                            p = detect,
                                            start = start_indexes[1:nseason], # Start and end arguments allow you to provide ragged mtx data
                                            end = end_indexes[1:nseason])


    }

    # Define Priors
    for (i in 1:(nseason-1)){
      persist_intercept[i] ~ dunif(0,1)
      logit(persist[i]) <- logit(persist_intercept[i]) + logit(beta) * i
    }

    # priors
    beta ~ dunif(0,1)
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
  compile_model <- init_model(n = sample_size, year = n_year, visit = n_visit, model_obj = dynoccmod_code)

  ## Map across scenarios for EMU and reps (multiple scenarios with the same sample size for each EMU)
  purrr::pmap(map_df, function(high_name, psi, phi, sd_phi, sd_gamma, p, high_hex_count, perc_red,
                               high_n, low_n, total_n, low_name, low_hex_count, low_psi, scenario_id, seed, year, visit){

    set.seed(seed)

    browser()
    for (i in 1:reps){

      #simulate high occupancy
      high_data <- sim_dataset(psi = psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                               n_sites = high_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
      append(c("sim_id" = high_name, "rep" = i))

      #simulate low occupancy
      low_data <- sim_dataset(psi = low_psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                              n_sites = low_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
        append(c("sim_id" = low_name, "rep" = i))

      # save data out
      saveRDS(high_data, paste0(here::here("data/nimble/emu_simulated_data/"), "/", high_name, "_", i, "_", scenario_id, "_simdata.rds"))
      saveRDS(low_data, paste0(here::here("data/nimble/emu_simulated_data/"), "/", low_name, "_", i, "_", scenario_id, "_simdata.rds"))

      #sample observed occupancy as model input
      sample_occ <- dplyr::bind_rows(sample_data(high_data$obs_occ, high_n),
                              sample_data(low_data$obs_occ, low_n)) %>%
        dplyr::arrange(visit)

      occ_array <- sample_occ %>%
        dplyr::select(-site_id) %>%
        split(sample_occ$visit) %>%
        purrr::map( ~ .x |> dplyr::select(-visit) |> as.matrix()) %>%
        simplify2array()

      # list of new initialized variables
      new_inits <-  list(
        colonize = 0.5,
        init_occ = 0.5,
        detect = 0.5,
        persist_intercept = rep(0.5, (year-1)),
        beta = 0.5)

      # update model with data
      compile_model$mod$y <- occ_array
      #compile_model$setData("y")
      fit <- nimble::runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                             samplesAsCodaMCMC = TRUE,
                             inits = new_inits)

      summary <- MCMCvis::MCMCsummary(fit, probs = c(0.025, 0.5, 0.95, 0.975)) |>
        dplyr::mutate(high_name = high_name, rep = i) |>
        tibble::rownames_to_column(var = "parameter")

      readr::write_csv(summary, paste0(here::here("data/nimble/emu_summary/"), "/", high_name, "_", i, "_", scenario_id, "_summary.csv"))
      saveRDS(fit, paste0(here::here("data/nimble/emu_posterior/"), "/", high_name, "_", i, "_", scenario_id, "_posterior.rds"))
    }
  }, year = n_year, visit = n_visit)
}

# initialize save out directories
dir.create(here::here("data/nimble/emu_summary"), recursive = TRUE)
dir.create(here::here("data/nimble/emu_posterior"))
dir.create(here::here("data/nimble/emu_simulated_data"))

ncores <- 11
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map_names', 'sim_dataset', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})

simn = 200

chunk_list <- unique(sim_map_names$chunk_num)
results <- parLapply(cl, chunk_list, fit_model_reps,
       reps = simn, n_year = 10, n_visit = 2, data = sim_map_names)


lapply(chunk_list, fit_model_reps,
       reps = 100, n_year = 10, n_visit = 2, data = sim_map_names)

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
# sim_map_metadata <- read.csv(here::here("data/sim_map_metadata.csv"))

# read in estimates
nimble_output <- purrr::map_dfr(list.files(here::here("data/nimble/emu_summary/"), full.names = TRUE), ~read.csv(.x) %>%
                                  filter(parameter == "perc_change") %>%
                                  select(mean, ci025 = X2.5., ci97.5 = X97.5., high_name, rep))

# get dataframe of sims and reps we've already done
sim_data_files <- data.frame(files = list.files(here::here("data/nimble/emu_simulated_data/")),
                         file_paths = list.files(here::here("data/nimble/emu_simulated_data/"), full.names = TRUE)) %>%
                         separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
    unite("high_name", c("EMU", "sim_id")) %>%
  mutate(scenario_id = as.integer(scenario_id)) %>%
  left_join(sim_map_names %>% select(high_name, low_name, scenario_id)) %>%
  filter(!is.na(low_name)) %>%
  select(-file_paths)

# get true occupancy
true_occ <- pmap_dfr(sim_data_files %>% select(high_name, rep, low_name, scenario_id), function(high_name, rep, low_name, scenario_id){
    #browser()
    readRDS(paste0(here::here("data/nimble/emu_simulated_data/"), "/", high_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ %>%
      bind_rows(readRDS(paste0(here::here("data/nimble/emu_simulated_data/"), "/", low_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ) %>%
      select(-site_id) %>%
      ungroup() %>%
      summarize(across(everything(), mean)) %>%
      mutate(rep = rep, high_name = high_name, low_name = low_name, scenario_id = scenario_id)
}) %>%
  mutate(true_perc_change = (t10-t1)/t1) %>%
  select(high_name, low_name, rep, true_perc_change)

# this works if things are no longer in progress
# true_occ <- pmap(sim_map_names %>% select(high_name, low_name), function(high_name, low_name){
#   map_dfr(1:simn, function(high_name, low_name){
#
#     readRDS(paste0(here::here("data/nimble/emu_simulated_data/"), "/", high_name, "_", .x, "_simdata.rds"))$true_occ %>%
#             bind_rows(readRDS(paste0(here::here("data/nimble/emu_simulated_data/"), "/", low_name, "_", .x, "_simdata.rds"))$true_occ) %>%
#             select(-site_id) %>%
#             ungroup() %>%
#             summarize(across(everything(), mean)) %>%
#             mutate(rep = .x)
#           }, high_name = high_name, low_name = low_name) %>%
#     mutate(high_name = high_name, low_name = low_name)
# })


# perc_change power checks
perc_change_check <- nimble_output %>%
  left_join(true_occ %>% select(high_name, true_perc_change, rep) %>%
              mutate(rep = as.integer(rep))) %>%
  # pivot_longer(starts_with("ci"), names_to = "ci_type", values_to = "ci_value") %>%
  # mutate(ci_low = (mean - abs(ci_value)), ci_high = (mean + abs(ci_value))) %>%
  rowwise() %>%
  mutate(ci_two_tail = between(true_perc_change, ci025, ci97.5) & !between(0,  ci025, ci97.5)) %>%
  group_by(high_name) %>%
  summarize(ci_two_tail = sum(ci_two_tail)/n(),
            rep_count = n()) %>%
  left_join(sim_map_names %>% select(high_name, total_n, psi, phi, p) %>% distinct()) %>%
  group_by(psi, p, phi) %>%
  mutate(line_id = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

readr::write_csv(perc_change_check, here::here("data/nimble_power_check.csv"))


### Let's look at how the estimates converge on the true mean trend across sample sizes
# get mean true occurrence across replicates for each scenario/emu/sample size
library(ggplot2)

mean_true_trend <- true_occ %>%
  left_join(sim_map_names %>%
              select(high_name, total_n, psi, p, phi) %>% distinct()) %>%
  group_by(high_name, low_name, total_n, psi, p, phi) %>%
  summarize(mean_true_trend = mean(true_perc_change)) %>%
  group_by(psi, p, phi) %>%
  mutate(scenario = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

est_trend_reps <-  nimble_output %>%
  select(mean, high_name, rep) %>%
  left_join(sim_map_names %>%
              select(high_name, low_name, total_n, psi, p, phi) %>% distinct()) %>%
  group_by(psi, p, phi) %>%
  mutate(scenario = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

mean_est_trend <- est_trend_reps %>%
  group_by(high_name, low_name, total_n, psi, p, phi, scenario, emu) %>%
  summarize(mean_est_trend = mean(mean))

mean_precision <- nimble_output %>%
  mutate(width = ci97.5 - ci025) %>%
  group_by(high_name) %>%
  summarize(mean_width = mean(width)) %>%
  left_join(sim_map_names %>%
              select(high_name, low_name, total_n, psi, p, phi) %>% distinct()) %>%
  group_by(psi, p, phi) %>%
  mutate(scenario = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

est_trend_reps %>%
  filter(scenario == 7, emu == "BRE") %>%
  ggplot() +
  geom_point(aes(x = total_n, y = mean)) +
  ylim(c(-2, 2)) +
  geom_line(data = mean_true_trend %>% filter(scenario == 7, emu == "BRE"),
            aes(x = total_n, y = mean_true_trend), color = "red") +
  geom_line(data = mean_est_trend %>% filter(scenario == 7, emu == "BRE"),
            aes(x = total_n, y = mean_est_trend), color = "blue") +
  theme_classic()

# # for when the whole thing doesn't run in one go
# missing_runs <- perc_change_check %>%
#   filter(rep_count == 100) %>%
#   select(-c(sim_num, ci_two_tail, rep_count, line_id, emu)) %>%
#   mutate(finished = TRUE) %>%
#   right_join(sim_map_names) %>%
#   filter(is.na(finished)) %>%
#   select(-finished)
#
# ncores <- 50
# cl <- makeCluster(ncores, type = "PSOCK")
# clusterExport(cl, c('init_model', 'missing_runs', 'sim_dataset', 'sample_data'))
# capture <- clusterEvalQ(cl, {
#   library(nimbleEcology)
#   library(magrittr)
#   library(purrr)
#   library(dplyr)
# })
#
# simn = 100
#
# # this uses the scenario as the level of parallelization, not chunk, so have to change
# # the dataframe filter statement for it to work
# chunk_list <- unique(missing_runs$high_name)
# results <- parLapply(cl, chunk_list, fit_model_reps,
#                      reps = simn, n_year = 10, n_visit = 2,
#                      data = missing_runs)
#
