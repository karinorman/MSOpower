library(dplyr)
library(tidyr)
library(purrr)
library(furrr)

source(here::here("R/sim_dataset.R"))
source(here::here("R/sim_data_equilib.R"))

# We want to compare the two approaches for simulating data by visualizing the true occupancy trend
sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = c(0.8, 0.8, .8)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  # and give each unique combination an ID
  mutate(sim_num = paste0("sim_", row_number()),
         occupancy = ifelse(psi == 0.03, "low", "high"),
         # Use site numbers from the Upper Gila Mountains EMU
         n_sites = ifelse(psi == 0.03, 3553, 5016)) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  mutate(sd_phi = 0.04, sd_gamma = 0.01, perc_red = 0.25) %>%
  select(psi, phi, sd_phi, sd_gamma, p, n_sites, perc_red, sim_num, occupancy)

simn <- 100

########################################
######### Simulate using final #########
########### trend approach #############
########################################

set.seed(42)
plan(multisession, workers = 14)
sim_list_recurs <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-sim_num, -occupancy),
                                                   sim_dataset, nyear = 10, n_vis = 2,
                                                   .options=furrr_options(seed = TRUE)) %>%
                  set_names(sim_scenarios$sim_num)) %>%
  set_names(paste0("rep", 1:simn))


# reorder so top level of nested list is a sim scenario
sim_list_scenario_recurs <- map(sim_scenarios$sim_num, function(scenario) {
  map(1:simn, ~pluck(sim_list_recurs, .x, scenario)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios$sim_num)

#get true occurrence for each rep and sim
true_occ_recurs <- map_dfr(sim_scenarios$sim_num, function(scenario){
  map_dfr(1:simn, ~pluck(sim_list_scenario_recurs, scenario, .x, "true_occ") %>%
            select(-site_id) %>%
            ungroup() %>%
            summarize(across(everything(), mean)) %>%
            mutate(rep = .x)) %>%
    mutate(sim_num = scenario)
})

true_occ_stats_recurs <- true_occ_recurs %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(sim_num, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(sim_scenarios)

true_occ_stats_recurs %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(color = sim_num)) +
  theme_classic() +
  #facet_wrap(~emu, scales = "free") +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted") +
  geom_hline(yintercept = 0.45, linetype = "dotted") +
  geom_hline(yintercept = 0.6, linetype = "dotted") +
  geom_hline(yintercept = 0.43, linetype = "dotted") +
  geom_hline(yintercept = 0.3225, linetype = "dotted")


##############################################
######### Simulate using equilibrium #########
############## trend approach ################
##############################################

set.seed(42)
plan(multisession, workers = 14)
sim_list_equil <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-sim_num, -occupancy),
                                                   sim_dataset_equil, nyear = 10, n_vis = 2,
                                                   .options=furrr_options(seed = TRUE)) %>%
                         set_names(sim_scenarios$sim_num)) %>%
  set_names(paste0("rep", 1:simn))

# reorder so top level of nested list is a sim scenario
sim_list_scenario_equil <- map(sim_scenarios$sim_num, function(scenario) {
  map(1:simn, ~pluck(sim_list_equil, .x, scenario)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios$sim_num)

#get true occurrence for each rep and sim
true_occ_equil <- map_dfr(sim_scenarios$sim_num, function(scenario){
  map_dfr(1:simn, ~pluck(sim_list_scenario_equil, scenario, .x, "true_occ") %>%
            select(-site_id) %>%
            ungroup() %>%
            summarize(across(everything(), mean)) %>%
            mutate(rep = .x)) %>%
    mutate(sim_num = scenario)
})

true_occ_stats_equil <- true_occ_equil %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(sim_num, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(sim_scenarios)

true_occ_stats_equil %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(color = sim_num)) +
  theme_classic() +
  #facet_wrap(~emu, scales = "free") +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted") +
  geom_hline(yintercept = 0.45, linetype = "dotted") +
  geom_hline(yintercept = 0.6, linetype = "dotted") +
  geom_hline(yintercept = 0.43, linetype = "dotted") +
  geom_hline(yintercept = 0.3225, linetype = "dotted")

true_occ_stats_equil %>%
  filter(psi == 0.03) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(color = sim_num)) +
  theme_classic() +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted")
