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

# Let's start with just one EMU
sim_map_names <- sim_map_names %>%
  filter(stringr::str_detect(high_name, "^SRM"))

## Fixed study characteristics
nyear = 10
n_vis = 2
simn = 300

###############################################################
################## Recursive Trend Sim ########################
###############################################################

# initialize save out directories
dir.create(here::here("data/nimble/emu_summary"), recursive = TRUE)
dir.create(here::here("data/nimble/emu_posterior"))
dir.create(here::here("data/nimble/emu_simulated_data"))

ncores <- 13
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
       method = "recursive", save_ending = "")


# lapply(chunk_list, fit_model_reps,
#        reps = simn, n_year = 10, n_visit = 2, data = sim_map_names, method = "recursive", save_ending = "")

#################################################################
################## Equilibrium Trend Sim ########################
#################################################################

path_ending <- "_equil"

# initialize save out directories
dir.create(paste0(here::here("data/nimble/emu_summary"), path_ending), recursive = TRUE)
dir.create(paste0(here::here("data/nimble/emu_posterior"), path_ending))
dir.create(paste0(here::here("data/nimble/emu_simulated_data"), path_ending))

ncores <- 13
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map_names', 'sim_dataset_equil', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- unique(sim_map_names$chunk_num)
results <- parLapply(cl, chunk_list, fit_model_reps,
                     reps = simn, n_year = 10, n_visit = 2, data = sim_map_names, 
                     method = "equilibrium", save_ending = path_ending)


# lapply(chunk_list, fit_model_reps,
#        reps = simn, n_year = 10, n_visit = 2, data = sim_map_names, 
#        method = "equilibrium", save_ending = path_ending)

#######################################################################
################## Constant Survival Trend Sim ########################
#######################################################################

path_ending <- "_constphi"

# initialize save out directories
dir.create(paste0(here::here("data/nimble/emu_summary"), path_ending), recursive = TRUE)
dir.create(paste0(here::here("data/nimble/emu_posterior"), path_ending))
dir.create(paste0(here::here("data/nimble/emu_simulated_data"), path_ending))

sim_map_constphi <- sim_map_names %>%
  mutate(low_gamma = ((1-perc_red)*low_psi*(1-phi))/(1 - ((1-perc_red)*low_psi)),
                high_gamma = ((1-perc_red)*psi*(1-phi))/(1 - ((1-perc_red)*psi)))

ncores <- 13
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map_constphi', 'sim_dataset_constphi', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- unique(sim_map_constphi$chunk_num)
results <- parLapply(cl, chunk_list, fit_model_reps,
                     reps = simn, n_year = 10, n_visit = 2, data = sim_map_constphi, 
                     method = "const_phi", save_ending = path_ending)


# lapply(chunk_list, fit_model_reps,
#        reps = simn, n_year = 10, n_visit = 2, data = sim_map_constphi, 
#        method = "const_phi", save_ending = path_ending)

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
