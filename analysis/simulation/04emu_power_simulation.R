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


###########################################################
################## Processing runs ########################
###########################################################

nimble_output <- purrr::map_dfr(list.files(paste0(path, "/emu_summary/"), full.names = TRUE), ~read.csv(.x) %>%
                                  filter(parameter == "perc_change") %>%
                                  select(mean, ci025 = X2.5., ci97.5 = X97.5., high_name, rep))

# get dataframe of sims and reps we've already done
sim_data_files <- data.frame(files = list.files(paste0(path, "/emu_simulated_data/")),
                             file_paths = list.files(paste0(path, "/emu_simulated_data/"), full.names = TRUE)) %>%
  separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
  unite("high_name", c("EMU", "sim_id")) %>%
  mutate(scenario_id = as.integer(scenario_id)) %>%
  left_join(sim_map_names %>% select(high_name, low_name, scenario_id)) %>%
  filter(!is.na(low_name)) %>%
  select(-file_paths)

# are any scenarios missing?
# get scenarios for which we didn't get all the reps
missing_scenarios <- sim_data_files %>%
  group_by(high_name, scenario_id) %>%
  summarize(reps = n_distinct(rep)) %>%
  select(-scenario_id) %>%
  filter(reps == 200) %>%
  mutate(complete = "yes") %>%
  right_join(sim_map_names) %>%
  filter(is.na(complete)) %>%
  select(-c(complete, reps))

# Let's do the missing ones
# want a worker for each row instead of using the chunk approach
# missing_scenarios_expt <- missing_scenarios %>%
#   ungroup() %>%
#   mutate(chunk_num = row_number())
#
# ncores <- 17
# cl <- makeCluster(ncores, type = "PSOCK")
# clusterExport(cl, c('init_model', 'missing_scenarios_expt', 'sim_dataset', 'sample_data'))
# capture <- clusterEvalQ(cl, {
#   library(nimbleEcology)
#   library(magrittr)
#   library(purrr)
#   library(dplyr)
# })
#
#
# chunk_list <- unique(missing_scenarios_expt$chunk_num)
# results <- parLapply(cl, chunk_list, fit_model_reps,
#                      reps = simn, n_year = 10, n_visit = 2, data = missing_scenarios_expt,
#                      method = "recursive", path = path, save_ending = "")



# get true occurrence at the scenario level (paired high and low)
true_occ_paired <- pmap_dfr(sim_data_files %>% select(high_name, rep, low_name, scenario_id),
                            function(high_name, rep, low_name, scenario_id){

                              #browser()
                              readRDS(paste0(path, "/emu_simulated_data", "/", high_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ %>%
                                bind_rows(readRDS(paste0(path, "/emu_simulated_data", "/", low_name, "_", rep, "_", scenario_id, "_simdata.rds"))$true_occ) %>%
                                select(-site_id) %>%
                                ungroup() %>%
                                summarize(across(everything(), mean)) %>%
                                mutate(rep = rep, high_name = high_name, low_name = low_name, scenario_id = scenario_id)
                            }) %>%
  mutate(true_perc_change = (t10-t1)/t1)

readr::write_csv(true_occ_paired, here::here("data/true_occ_paired_emu.csv"))

# perc_change power checks
perc_change_check <- nimble_output %>%
  left_join(true_occ_paired %>% select(high_name, true_perc_change, rep) %>%
              mutate(rep = as.integer(rep))) %>%
  rowwise() %>%
  mutate(ci_two_tail = between(true_perc_change, ci025, ci97.5) & !between(0,  ci025, ci97.5),
         bias = true_perc_change - mean) %>%
  group_by(high_name) %>%
  summarize(ci_two_tail = sum(ci_two_tail)/n(),
            bias = mean(bias),
            rep_count = n()) %>%
  left_join(sim_map_names %>% select(high_name, total_n, psi, phi, p) %>% distinct()) %>%
  group_by(psi, p, phi) %>%
  mutate(line_id = cur_group_id()) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE)

readr::write_csv(perc_change_check, here::here("data/nimble_power_check_emu.csv"))

### Export dataframe of sample sizes and power values for before and after the threshold is crossed

pre_df <- perc_change_check %>%
  select(emu, ci_two_tail, total_n, psi, phi, p) %>%
  mutate(distance = ci_two_tail - 0.90) %>%
  filter(distance < 0) %>%
  group_by(emu, psi, phi, p) %>%
  filter(abs(distance) == min(abs(distance))) %>%
  rename(pre = ci_two_tail, pre_n = total_n) %>%
  select(-distance)

# get before and after crossing the power threshold
post_df <- perc_change_check %>%
  select(emu, line_id, ci_two_tail, total_n, psi, phi, p) %>%
  group_by(emu, psi, phi, p) %>%
  mutate(past_threshold = (ci_two_tail > 0.9)) %>%
  filter(past_threshold == TRUE) %>%
  filter(total_n == min(total_n)) %>%
  rename(post = ci_two_tail, post_n = total_n) %>%
  select(-past_threshold) %>%
  ungroup()

# get data frame with pre and post power threshold points
# calculate the exact threshold
threshold_df <- perc_change_check %>%
  select(emu, line_id, ci_two_tail, total_n, psi, phi, p) %>%
  rename(pre = ci_two_tail, pre_n = total_n) %>%
  left_join(post_df) %>%
  mutate(distance = post_n - pre_n) %>%
  filter(distance > 0) %>%
  group_by(emu, psi, phi, p) %>%
  filter(distance == min(distance)) %>%
  select(-distance) %>%
  mutate(slope = (post - pre)/(post_n - pre_n),
         intercept = post - (slope*post_n),
         threshold = (.9 - intercept)/slope)

readr::write_csv(threshold_df, here::here("data/emu_power_thresholds.csv"))
