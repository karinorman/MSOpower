library(dplyr)
library(nimble)
library(nimbleEcology)
library(MCMCvis)

# model object, which can be made outside of iterating and passed in

init_model <- function(n, year, visit, model_obj){

  nsite <- n
  nseason <- year
  nrep <- visit

  # make observed occurrence an array with site x year(season) x visit(rep)
  obs_occ_init <- array(sample(c(0,1), (nsite * nseason * nrep), replace = TRUE), c(nsite, nseason, nrep))

  # Build the model
  mod <- nimbleModel(
    code = model_obj,
    constants = list(nsite = nsite, nrep = nrep, nseason = nseason,
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

  return(complist)
}

model_check_nimble <- function(high_name, low_name, high_n, low_n, sample_size, repn, high_occ, low_occ, compile_model) {

  print(c(high_name, repn))

  # get appropriate sample size from the obs_occ
  if(low_n > n_distinct(low_occ$site_id)){
    low_n = n_distinct(low_occ$site_id)
  }

  obs_occ_df <- bind_rows(low_occ %>% filter(site_id %in% sample(unique(low_occ$site_id), low_n, replace = FALSE)),
                          high_occ %>% filter(site_id %in% sample(unique(high_occ$site_id), high_n, replace = FALSE))) %>%
    arrange(visit)

  obs_occ_array <- obs_occ_df %>%
    select(-c(site_id, landtype)) %>%
    split(obs_occ_df$visit) %>%
    map(., ~ .x %>% select(-visit) %>% as.matrix()) %>%
    simplify2array()

  # list of new initialized variables
  new_inits <-  list(
    colonize = 0.5,
    init_occ = 0.5,
    detect = 0.5,
    persist = rep(0.5, (nseason-1)))

  # update model with data
  compile_model$mod$y <- obs_occ_array
  #compile_model$setData("y")
  fit <- runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                 samplesAsCodaMCMC = TRUE,
                 inits = new_inits)

  # Run the MCMC
  samples <- runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                     samplesAsCodaMCMC = TRUE)

  summary <- MCMCsummary(samples, probs = c(0.025, 0.5, 0.95, 0.975)) %>%
    mutate(high_name = high_name, low_name = low_name, rep = repn)

  readr::write_csv(summary, paste0(here::here("data/nimble/emu_summary/"), high_name, "_", repn, "_summary.csv"))
  saveRDS(samples, paste0(here::here("data/nimble/emu_posterior/"), high_name, "_", repn, "_posterior.rds"))
}

# initialize save out directories
dir.create(here::here("data/nimble/emu_summary"), recursive = TRUE)
dir.create(here::here("data/nimble/emu_posterior"))


#1: create model object
# this model has survival that can change between years and constant colonization
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

# for each sample size:
# 2: initialize and compile model
purrr::map(unique(BRE_sim_occ$total_n), function(samp_size, n_year, n_visit){
  browser()

  comp_model <- init_model(n = samp_size, year = n_year, visit = n_visit, model_obj = dynoccmod_code)

  # get only reps for the sample size that we have an initialized model for
  rep_df <- BRE_sim_occ %>%
    select(high_name, low_name, high_n, low_n, sample_size = total_n, repn = rep, high_occ, low_occ) %>%
    filter(sample_size == samp_size)

  # for each rep in each sample size
  # 3: initialize and update model, get posteriors
  plan(multisession, workers = 70)
  power_check_list <- furrr::future_pmap(rep_df, model_check_nimble, compile_model = comp_model,
                                         .options=furrr_options(seed = TRUE))

}, n_year = nseason, n_visit = n_vis # the arguments that are constant
)

# let's get a single EMU example
BRE_scenarios <- sim_map_occ %>% select(high_name) %>%
  tidyr::separate(high_name, into = c("EMU", "scenario"), sep = "_", remove = FALSE) %>%
  filter(EMU == "BRE") %>%
  pull(high_name)

BRE_sim_occ <- sim_map_occ %>% filter(high_name %in% BRE_scenarios)






