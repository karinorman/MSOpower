library(dplyr)
library(tidyr)
library(nimble)
library(nimbleEcology)
library(MCMCvis)
library(parallel)
library(purrr)

## Fixed study characteristics
nyear = 10
n_vis = 2
simn = 100

load("data/sim_map_hier.rda")

# now instead of generating each EMU separately, let's do the whole landscape in one go,
# maintaining the high/low ratio
hex_count <- sim_map_hier %>%
  select(emu, high_hex_count, low_hex_count) %>%
  distinct()

landscape_hex_count <- hex_count %>%
  select(-emu) %>%
  summarise(across(everything(), ~ sum(., na.rm = TRUE)))

# let's get all the rows to fill out a full curve
set.seed(5242342)

sim_df <-  sim_map_hier %>%
  select(chunk_num = scenario_id, total_n, high_n, low_n, psi, phi, sd_phi, sd_gamma, p, perc_red, low_psi) %>%
  # this one is particularly biased
  #filter(psi == 0.6, phi == 0.8, p == 0.4, total_n < 3500) %>%
  distinct() %>%
  mutate(high_hex_count = landscape_hex_count$high_hex_count,
         low_hex_count = landscape_hex_count$low_hex_count,
         seed = 1 + floor(runif(n()) * 100000))

###############################################################
######### Whole landscape, no hierarchical structure ##########
###############################################################

source(here::here("R/sim_dataset.R"))
source(here::here("R/fit_model_reps_working.R"))

path <- here::here("data/nimble/whole_landscape")

ncores <- 7
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_df', 'sim_dataset', 'sample_data', 'path'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- unique(sim_df$chunk_num)
results <- parLapply(cl, chunk_list, fit_model_reps_working,
                     reps = simn, n_year = nyear, n_visit = n_vis, data = sim_df,
                     method = "recursive", path = path, save_ending = "")


results <- lapply(chunk_list, fit_model_reps_working,
                     reps = simn, n_year = nyear, n_visit = n_vis, data = sim_df,
                     method = "recursive", path = path, save_ending = "")


###############################################################
############## Now with hierarchical structure ################
###############################################################
source(here::here("R/fit_model_reps_hier_working.R"))

path <- here::here("data/nimble/whole_landscape_hier")

ncores <- 13
cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, c('init_model', 'sim_df', 'sim_dataset', 'sample_data_hier', 'path', 'hex_count'))
capture <- clusterEvalQ(cl, {
  library(nimbleEcology)
  library(magrittr)
  library(purrr)
  library(dplyr)
})


chunk_list <- sort(unique(sim_df$chunk_num), decreasing = TRUE)
results <- parLapply(cl, chunk_list, fit_model_reps_hier_working,
                     reps = simn, n_year = nyear, n_visit = n_vis, data = sim_df,
                     path = path, save_ending = "")

# results <- lapply(chunk_list, fit_model_reps_hier_working,
#                   reps = simn, n_year = nyear, n_visit = n_vis, data = sim_df,
#                   path = path, save_ending = "")

##############################
###### Read in results #######
##############################

nimble_output <- purrr::map_dfr(list.files(paste0(path, "/emu_summary/"), full.names = TRUE), ~read.csv(.x) %>%
                                  filter(parameter == "perc_change") %>%
                                  select(mean, ci025 = X2.5., ci97.5 = X97.5., chunk_num, rep))

# get dataframe of sims and reps we've already done
sim_data_files <- data.frame(files = list.files(paste0(path, "/emu_simulated_data/"))) %>%
  separate(files, c("occ_type", "chunk_num", "rep")) %>%
  mutate(across(c(chunk_num, rep), as.numeric)) #%>%
#pivot_wider(names_from = occ_type, values_from = )

# get true occurrence for each scenario and rep
true_occ_paired <- pmap_dfr(sim_data_files %>% select(chunk_num, rep) %>% distinct(),
                            function(chunk_num, rep){

                              #browser()
                              readRDS(paste0(path, "/emu_simulated_data", "/", "high", "_", chunk_num, "_", rep, "_simdata.rds"))$true_occ %>%
                                bind_rows(readRDS(paste0(path, "/emu_simulated_data", "/", "low", "_", chunk_num, "_", rep, "_simdata.rds"))$true_occ) %>%
                                select(-site_id) %>%
                                ungroup() %>%
                                summarize(across(everything(), mean)) %>%
                                mutate(rep = rep, chunk_num = chunk_num)
                            }) %>%
  mutate(true_perc_change = (t10-t1)/t1)

perc_change_check <- nimble_output %>%
  left_join(true_occ_paired %>% select(chunk_num, true_perc_change, rep) %>%
              mutate(across(c(chunk_num, rep), as.numeric))) %>%
  rowwise() %>%
  mutate(ci_two_tail = between(true_perc_change, ci025, ci97.5) & !between(0,  ci025, ci97.5),
         bias = true_perc_change - mean) %>%
  group_by(chunk_num) %>%
  summarize(ci_two_tail = sum(ci_two_tail)/n(),
            bias = mean(bias),
            rep_count = n()) %>%
  left_join(sim_df %>% select(chunk_num, total_n, psi, phi, p) %>% distinct())

readr::write_csv(perc_change_check, here::here("data/nimble_power_check_hier_landscape.csv"))

## RUN MISSING ONES ##
# comp_chunk <- perc_change_check %>%
#   filter(rep_count == 100) %>%
#   pull(chunk_num)
#
# missing_chunks <- chunk_list[!chunk_list %in% comp_chunk]

# get before and after crossing the power threshold
pre_df <- perc_change_check %>%
  select(ci_two_tail, total_n, psi, phi, p) %>%
  mutate(distance = ci_two_tail - 0.90) %>%
  filter(distance < 0) %>%
  group_by(psi, phi, p) %>%
  filter(abs(distance) == min(abs(distance))) %>%
  rename(pre = ci_two_tail, pre_n = total_n) %>%
  select(-distance)

post_df <- perc_change_check %>%
  select(ci_two_tail, total_n, psi, phi, p) %>%
  mutate(distance = ci_two_tail - 0.90) %>%
  filter(distance > 0) %>%
  group_by(psi, phi, p) %>%
  filter(distance == min(distance)) %>%
  rename(post = ci_two_tail, post_n = total_n) %>%
  select(-distance) %>%
  # for some scenarios more than one sample size has the same power, need to filter
  group_by(psi, phi, p) %>%
  filter(post_n == min(post_n))

# get data frame with pre and post power threshold points
# calculate the exact threshold
threshold_df <- left_join(pre_df, post_df) %>%
  mutate(slope = (post - pre)/(post_n - pre_n),
         intercept = post - (slope*post_n),
         threshold = (.9 - intercept)/slope)

readr::write_csv(threshold_df, here::here("data/hier_power_thresholds.csv"))

