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
simn = 200

# directory to save outputs to
path <- here::here("data/nimble/srm_simulations")

###############################################################
################## Recursive Trend Sim ########################
###############################################################

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
       method = "recursive", path = path, save_ending = "")


# lapply(chunk_list, fit_model_reps,
#        reps = simn, n_year = 10, n_visit = 2, data = sim_map_names, method = "recursive", save_ending = "")

#################################################################
################## Equilibrium Trend Sim ########################
#################################################################

path_ending <- "_equil"

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
                     method = "equilibrium", path = path, save_ending = path_ending)


# lapply(chunk_list, fit_model_reps,
#        reps = simn, n_year = 10, n_visit = 2, data = sim_map_names,
#        method = "equilibrium", save_ending = path_ending)

#######################################################################
################## Constant Survival Trend Sim ########################
#######################################################################

path_ending <- "_constphi"

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
                     method = "const_phi", path = path, save_ending = path_ending)


# lapply(chunk_list, fit_model_reps,
#        reps = simn, n_year = 10, n_visit = 2, data = sim_map_constphi,
#        method = "const_phi", save_ending = path_ending)

###########################################################
################## Processing runs ########################
###########################################################

# read in estimates
nimble_output <- purrr::map_dfr(list.files(paste0(path, "/emu_summary/"), full.names = TRUE), ~read.csv(.x) %>%
                                  filter(parameter == "perc_change") %>%
                                  select(mean, ci025 = X2.5., ci97.5 = X97.5., high_name, rep)) %>%
  mutate(type = "recursive") %>%
  bind_rows(purrr::map_dfr(list.files(paste0(path, "/emu_summary_equil/"), full.names = TRUE), ~read.csv(.x) %>%
                             filter(parameter == "perc_change") %>%
                             select(mean, ci025 = X2.5., ci97.5 = X97.5., high_name, rep)) %>%
              mutate(type = "equilibrium"),
            purrr::map_dfr(list.files(paste0(path, "/emu_summary_constphi/"), full.names = TRUE), ~read.csv(.x) %>%
                             filter(parameter == "perc_change") %>%
                             select(mean, ci025 = X2.5., ci97.5 = X97.5., high_name, rep)) %>%
              mutate(type = "constant_phi"))

# fully completed scenarios
comp_scenario <- nimble_output %>% 
  group_by(high_name, type) %>% 
  summarize(rep_count = n_distinct(rep)) %>% 
  filter(rep_count == 200)

# get dataframe of sims and reps we've already done
sim_data_files <- data.frame(files = list.files(paste0(path, "/emu_simulated_data/")),
                         file_paths = list.files(paste0(path, "/emu_simulated_data/"), full.names = TRUE)) %>%
  mutate(type = "recursive", ending = "") %>%
  bind_rows(data.frame(files = list.files(paste0(path, "/emu_simulated_data_equil/")),
                       file_paths = list.files(paste0(path, "/emu_simulated_data_equil/"), full.names = TRUE)) %>%
              mutate(type = "equilibrium", ending = "_equil")) %>%
  bind_rows(data.frame(files = list.files(paste0(path, "/emu_simulated_data_constphi/")),
                       file_paths = list.files(paste0(path, "/emu_simulated_data_constphi/"), full.names = TRUE)) %>%
              mutate(type = "constant_phi", ending = "_constphi")) %>%
  separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
  unite("high_name", c("EMU", "sim_id")) %>%
  mutate(scenario_id = as.integer(scenario_id)) %>%
  left_join(sim_map_names %>% select(high_name, low_name, scenario_id)) %>%
  filter(!is.na(low_name)) %>%
  select(-file_paths)

# get true occupancy
# for each simulation setting
true_occ_indiv <- pmap_dfr(sim_data_files %>%
                              select(high_name, rep, low_name, scenario_id, type, ending) %>%
                              pivot_longer(cols = c(low_name, high_name), names_to = "initial_occ", values_to = "name") %>%
                              select(name, rep, scenario_id, type, ending, initial_occ),
                            function(name, rep, scenario_id, type, ending, initial_occ){

                              #browser()
                              readRDS(paste0(path, "/emu_simulated_data", ending, "/", name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ %>%
                                select(-site_id) %>%
                                ungroup() %>%
                                summarize(across(everything(), mean)) %>%
                                mutate(rep = rep, name = name, scenario_id = scenario_id, type = type, initial_occ = initial_occ)
                            }) %>%
  mutate(true_perc_change = (t10-t1)/t1) %>%
  left_join(sim_map_names %>%
              select(name = high_name, psi, phi, scenario_id) %>%
              bind_rows(sim_map_names %>%
                          select(name = low_name, psi = low_psi, phi, scenario_id)))

readr::write_csv(true_occ_indiv, here::here("data/true_occ_indiv_SRM.csv"))

# and paired for each scenario
true_occ_paired <- pmap_dfr(sim_data_files %>% select(high_name, rep, low_name, scenario_id, type, ending),
                                                                         function(high_name, rep, low_name, scenario_id, type, ending){

                                                                           #browser()
    readRDS(paste0(path, "/emu_simulated_data", ending, "/", high_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ %>%
      bind_rows(readRDS(paste0(path, "/emu_simulated_data",  ending, "/", low_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ) %>%
      select(-site_id) %>%
      ungroup() %>%
      summarize(across(everything(), mean)) %>%
      mutate(rep = rep, high_name = high_name, low_name = low_name, scenario_id = scenario_id, type = type)
}) %>%
  mutate(true_perc_change = (t10-t1)/t1)

readr::write_csv(true_occ_paired, here::here("data/true_occ_paired_SRM.csv"))

# perc_change power checks
perc_change_check <- nimble_output %>%
  left_join(true_occ_paired %>% select(high_name, type, true_perc_change, rep) %>%
              mutate(rep = as.integer(rep))) %>%
  # pivot_longer(starts_with("ci"), names_to = "ci_type", values_to = "ci_value") %>%
  # mutate(ci_low = (mean - abs(ci_value)), ci_high = (mean + abs(ci_value))) %>%
  rowwise() %>%
  mutate(ci_two_tail = between(true_perc_change, ci025, ci97.5) & !between(0,  ci025, ci97.5),
         bias = true_perc_change - mean) %>%
  group_by(high_name, type) %>%
  summarize(ci_two_tail = sum(ci_two_tail)/n(),
            bias = mean(bias),
            rep_count = n()) %>%
  left_join(sim_map_names %>% select(high_name, total_n, psi, phi, p) %>% distinct()) %>%
  group_by(psi, p, phi) %>%
  mutate(line_id = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

readr::write_csv(perc_change_check, here::here("data/nimble_power_check_SRM.csv"))


### When simulation gets interrupted, restart here
#################################################################
################## Equilibrium Trend Sim ########################
#################################################################

comp_equil_scenario <- comp_scenario %>% filter(type == "equilibrium") %>% pull(high_name)

sim_map_restart <- sim_map_names %>%
  filter(!high_name %in% comp_equil_scenario)

path_ending <- "_equil"

ncores <- 2
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_map_restart', 'sim_dataset_equil', 'sample_data'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- unique(sim_map_restart$chunk_num)
results <- parLapply(cl, chunk_list, fit_model_reps,
                     reps = simn, n_year = 10, n_visit = 2, data = sim_map_restart,
                     method = "equilibrium", path = path, save_ending = path_ending)

#######################################################################
################## Constant Survival Trend Sim ########################
#######################################################################

path_ending <- "_constphi"

comp_constphi_scenario <- comp_scenario %>% filter(type == "constant_phi") %>% pull(high_name)

sim_map_constphi <- sim_map_names %>%
  mutate(low_gamma = ((1-perc_red)*low_psi*(1-phi))/(1 - ((1-perc_red)*low_psi)),
         high_gamma = ((1-perc_red)*psi*(1-phi))/(1 - ((1-perc_red)*psi))) %>%
  filter(!high_name %in% comp_constphi_scenario)

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
                     method = "const_phi", path = path, save_ending = path_ending)




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
