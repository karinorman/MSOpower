############################################
########### Hierarchical Power #############
############################################

library(dplyr)
library(tidyr)
library(nimble)
library(nimbleEcology)
library(MCMCvis)
library(parallel)
library(purrr)

###########################################
######## Define simulation parameters #####
###########################################

## Fixed study characteristics
nyear = 10
n_vis = 2
simn = 200

load("data/sim_map_hier.rda")

# directory to save outputs to
path <- here::here("data/nimble/hierarchical_simulations")

##############################################
########### Model and Power check ############
##############################################

source(here::here("R/sim_dataset.R"))
source(here::here("R/fit_model_reps_hier.R"))

ncores <- 55
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map_names', 'sim_dataset', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- unique(sim_map_hier$scenario_id)
results <- parLapply(cl, chunk_list, fit_model_reps,
                     reps = simn, n_year = 10, n_visit = 2,
                     data = sim_map_hier)

lapply(chunk_list, fit_model_reps_hier,
       reps = 100, n_year = 10, n_visit = 2,
       data = sim_map_hier)

###########################################################
################## Processing runs ########################
###########################################################

# read in estimates
output_files <- data.frame(files = list.files(here::here("data/nimble/hier_summary/")),
                           file_paths = list.files(here::here("data/nimble/hier_summary/"), full.names = TRUE)) %>%
  separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
  select(EMU, file_paths)


nimble_output <- purrr::pmap_dfr(output_files, function(EMU, file_paths) {
  read.csv(file_paths) %>%
    filter(parameter == "perc_change") %>%
    select(mean, ci025 = X2.5., ci97.5 = X97.5., scenario_id, rep) %>%
    mutate(emu = EMU)
  })

# get dataframe of sims and reps we've already done
high_data_files <- data.frame(files = list.files(here::here("data/nimble/hier_simulated_data/")),
                             file_paths = list.files(here::here("data/nimble/hier_simulated_data/"), full.names = TRUE)) %>%
  separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
  unite("high_name", c("EMU", "sim_id")) %>%
  mutate(scenario_id = as.integer(scenario_id)) %>%
  left_join(sim_map %>% select(scenario_id, high_name, low_name)) %>%
  #left_join(sim_map %>% select(scenario_id, low_name), by = c("sim_id" = "low_name")) %>%
  filter(!is.na(low_name)) %>%
  select(-file_paths)

low_data_files <- data.frame(files = list.files(here::here("data/nimble/hier_simulated_data/")),
                              file_paths = list.files(here::here("data/nimble/hier_simulated_data/"), full.names = TRUE)) %>%
  separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
  unite("low_name", c("EMU", "sim_id")) %>%
  mutate(scenario_id = as.integer(scenario_id)) %>%
  #left_join(sim_map %>% select(scenario_id, high_name), by = c("sim_id" = "high_name", "scenario_id")) %>%
  left_join(sim_map %>% select(scenario_id, low_name, high_name)) %>%
  filter(!is.na(high_name)) %>%
  select(-file_paths)

sim_data_files <- full_join(high_data_files, low_data_files) %>%
  filter(!is.na(high_name), !is.na(low_name))

# get true occupancy
true_occ <- pmap_dfr(sim_data_files %>% select(high_name, rep, low_name, scenario_id), function(high_name, rep, low_name, scenario_id){

  readRDS(paste0(here::here("data/nimble/hier_simulated_data/"), "/", high_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ %>%
    bind_rows(readRDS(paste0(here::here("data/nimble/hier_simulated_data/"), "/", low_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ) %>%
    select(-site_id) %>%
    ungroup() %>%
    summarize(across(everything(), mean)) %>%
    mutate(rep = rep, high_name = high_name, low_name = low_name)
}) %>%
  mutate(true_perc_change = (t10-t1)/t1) %>%
  select(high_name, low_name, rep, true_perc_change)

perc_change_check <- nimble_output %>%
  left_join(true_occ %>% select(scenario_id, true_perc_change, rep) %>%
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
