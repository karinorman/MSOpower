#################################################
##### Check power to detect common estimates ####
#################################################

library(dplyr)
library(tidyr)
library(purrr)

load(here::here("data/sim_map_names.rda"))

ex_est <- read.csv(here::here("data/nimble/emu_summary/UGM_119_80_375_summary.csv"))


# get estimates of parameters of interest
estimates <- purrr::map_dfr(list.files(here::here("data/nimble/emu_summary/"), full.names = TRUE), ~read.csv(.x) %>% 
                 filter(parameter %in% c("init_occ", "detect", "colonize", "psi[10]")) %>%
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

# get true estimates of 
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
