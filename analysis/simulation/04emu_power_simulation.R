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

source(here::here("R/sim_dataset.R"))
source(here::here("R/fit_model_reps.R"))

load(here::here("data/sim_map_names.rda"))

# Get scenario mapping for all EMU's
sim_map_names <- sim_map_names 

## Fixed study characteristics
nyear = 10
n_vis = 2
simn = 200

# directory to save outputs to
path <- here::here("data/nimble/emu_simulations")

###############################################################
############ Power simulation for all EMU's ###################
################ using recursive approach #####################
###############################################################

ncores <- 45
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map_names', 'sim_dataset', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- unique(sim_map_names$chunk_num)
results <- parLapply(cl, chunk_list, fit_model_reps,
                     reps = simn, n_year = 10, n_visit = 2, data = sim_map_names,
                     method = "recursive", path = path, save_ending = "")

