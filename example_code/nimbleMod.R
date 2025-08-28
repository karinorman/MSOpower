#### Nimble Stuff ####

library(dplyr)
library(nimble)

load(here::here("data/sim_list.rda"))


sim_ex <- pluck(sim_list, "rep1", "BRE_1")
usethis::use_data(sim_ex)

rm(sim_list)


true_occ <- sim_ex$true_occ
obs_occ <- sim_ex$obs_occ


# parameters
nsite <- n_distinct(true_occ$site_id)
nseason <- 10
nrep <- 2

# make observed occurrence an array with site x year(season) x visit(rep)
obs_occ_array <- obs_occ %>%
  select(-site_id) %>%
  split(obs_occ$visit) %>%
  map(., ~ .x %>% select(-visit) %>% as.matrix()) %>%
  simplify2array()

#### NIMBLE model ####

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

# Build the model
mod <- nimbleModel(
  code = dynoccmod_code,
  constants = list(nsite = nsite, nrep = nrep, nseason = nseason,
                   start_indexes = rep(1, nseason),
                   end_indexes = rep(nrep, nseason)
                   ),
  data = list(y = obs_occ_array),
  inits = list(
    colonize = 0.5,
    init_occ = 0.5,
    detect = 0.5,
    persist = rep(0.5, (nseason-1))
  )
)

# shouldn't NA, infinite, or positive (that's a dist issue)
mod$calculate() #Non-NA means we're fully initialized


# Build an MCMC
conf <- configureMCMC(mod)
conf$addMonitors(c("psi", "perc_change"))
mcmc <- buildMCMC(conf)

# Compile
complist <- compileNimble(mod, mcmc)

# Run the MCMC
samples <- runMCMC(complist$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                   samplesAsCodaMCMC = TRUE)

# for updating
# update data and reinitialize
new_inits <-  list(
  colonize = 0.5,
  init_occ = 0.5,
  detect = 0.5,
  persist = rep(0.5, (nseason-1))

complist$mod$y <- new_y
complist$setData("y")
fit <- runMCMC(complist$mcmc, niter = 1000, nchains = 2, nburnin = 500,
        samplesAsCodaMCMC = TRUE,
        inits = new_inits)


# Get a summary df of posterior samples
summary <- MCMCsummary(samples)
