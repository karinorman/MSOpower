#### Let's check if a constant annual decline will match the model expectation ###

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(ggplot2)
library(patchwork)

source(here::here("R/sim_data_constphi.R"))

simn <- 100

sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.6, 0.6, 0.6), phi = c(.6, 0.3, 0.4, 0.8), p = rep(0.8, 4)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p)%>%
  mutate(sd_phi = 0.04, sd_gamma = 0.01) %>%
  mutate(occupancy = ifelse(psi == 0.03, "low", "high"),
         # Use site numbers from Southern Rocky Mountains EMU
         n_sites = ifelse(psi == 0.03, 11400, 7577),
         gamma = ifelse(psi == 0.03, 0.02, 0.03)) %>%
  group_by(phi, p) %>%
  mutate(simulation_id = cur_group_id()) %>%
  ungroup() %>%
  relocate(c(p, gamma, n_sites), .after = last_col()) %>%
  mutate(name = row_number())

sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = rep(0.8, 3)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  # and give each unique combination an ID
  mutate(occupancy = case_when(
    psi == 0.03 ~ "low",
    psi == 0.43 ~ "med",
    psi == 0.6 ~ "high"),
         # Use site numbers from Southern Rocky Mountains EMU
         n_sites = ifelse(psi == 0.03, 11400, 7577),
         gamma = (0.75*psi*(1-phi))/(1 - (0.75*psi)),
         sd_phi = 0.04, sd_gamma = 0.01) %>%
  group_by(phi, p) %>%
  mutate(simulation_id = paste0(occupancy, "_", cur_group_id())) %>%
  ungroup()

########################################
######### Simulate using final #########
########### trend approach #############
########################################

plan(multisession, workers = 6)
sim_list_recurs <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(psi, phi, sd_phi, sd_gamma, p, gamma, n_sites),
                                                   sim_dataset_constphi, nyear = 10, n_vis = 2,
                                                   .options=furrr_options(seed = TRUE)) %>%
                         set_names(sim_scenarios$simulation_id)) %>%
  set_names(paste0("rep", 1:simn))

map_sim <- sim_scenarios %>%
  filter(occupancy %in% c("med", "high")) %>%
  rename(high_name = simulation_id) %>%
  select(-occupancy) %>%
  left_join(sim_scenarios %>%
              filter(occupancy == "low") %>%
              select(simulation_id, low_name = simulation_id, p, phi))

# true_occ_comb <- pmap_dfr(map_sim %>% select(high_name, low_name), function(high_name, low_name){
#
#   map(1:simn, function(rep, sim_id, high, low){
#
#     bind_rows(pluck(sim_list_recurs, rep, high)$true_occ,
#               pluck(sim_list_recurs, rep, low)$true_occ
#     ) %>%
#       select(-site_id) %>%
#       ungroup() %>%
#       summarize(across(everything(), mean)) %>%
#       mutate(rep = rep) %>%
#     mutate(high_name = high)
#   }, high = high_name, low = low_name)
# })
#
# true_occ_stats_comb <- true_occ_comb %>%
#   pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
#   select(-rep) %>%
#   group_by(high_name, time) %>%
#   summarize(mean = mean(occ),
#             lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
#             upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
#   ungroup() %>%
#   mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
#   left_join(map_sim)

true_occ_sep <-  map_dfr(c(map_sim$high_name, map_sim$low_name), function(name){

  map(1:simn, function(rep, id_name){

    pluck(sim_list_recurs, rep, id_name, "true_occ") %>%
      select(-site_id) %>%
      ungroup() %>%
      summarize(across(everything(), mean)) %>%
      mutate(rep = rep) %>%
      mutate(simulation_id = id_name)
  }, id_name = name)
})



true_occ_stats_sep <- true_occ_sep %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(simulation_id, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(sim_scenarios)

true_occ_stats_sep %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = simulation_id), alpha = 0.3) +
  geom_line(aes(color = as.factor(simulation_id))) +
  #facet_wrap(~sd_phi, nrow = 1) +
  theme_classic() +
  #facet_wrap(~emu, scales = "free") +
  #scale_color_discrete(name = "Sim Scenario") +
 # scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted") +
  geom_hline(yintercept = 0.43, linetype = "dotted") +
  geom_hline(yintercept = 0.3225, linetype = "dotted") +
  geom_hline(yintercept = 0.45, linetype = "dotted") +
  geom_hline(yintercept = 0.6, linetype = "dotted")


###############################
## OK values from the Viorel ##
###############################

sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = rep(0.6, 5), phi = c(1, 0.2, 0.5, 0.8, 0.95), p = rep(0.8, 5)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p)%>%
  mutate(sd_phi = 0.04, sd_gamma = 0.01, gamma = 0.4) %>%
  mutate(#occupancy = ifelse(psi == 0.03, "low", "high"),
    # Use site numbers from Southern Rocky Mountains EMU
    #n_sites = ifelse(psi == 0.03, 11400, 7577)#,
    n_sites = 7577
    #gamma = ifelse(psi == 0.03, 0.02, 0.03)
  ) %>%
  group_by(phi, p) %>%
  mutate(simulation_id = cur_group_id()) %>%
  ungroup() %>%
  relocate(c(p, gamma, n_sites), .after = last_col()) %>%
  mutate(name = row_number())

plan(multisession, workers = 5)
sim_list_recurs <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-simulation_id, -name),
                                                   sim_dataset_constphi, nyear = 10, n_vis = 2,
                                                   .options=furrr_options(seed = TRUE)) %>%
                         set_names(sim_scenarios$name)) %>%
  set_names(paste0("rep", 1:simn))


true_occ_sim <- map_dfr(sim_scenarios$simulation_id, function(simulation_id){

  map(1:simn, function(rep, sim_id){

    pluck(sim_list_recurs, rep, sim_id, "true_occ") %>%
      select(-site_id) %>%
      ungroup() %>%
      summarize(across(everything(), mean)) %>%
      mutate(rep = rep) %>%
      mutate(simulation_id = sim_id)
  }, sim_id = simulation_id)
})

mean_true_occ <- true_occ_sim %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(time, simulation_id) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t")))

mean_true_occ %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = simulation_id), alpha = 0.3) +
  geom_line(aes(color = as.factor(simulation_id))) +
  theme_classic() +
  #facet_wrap(~emu, scales = "free") +
  #scale_color_discrete(name = "Sim Scenario") +
  # scale_fill_discrete(name = "Sim Scenario") +
  #geom_hline(yintercept = 0.03, linetype = "dotted") +
  #geom_hline(yintercept = 0.0225, linetype = "dotted")# +
  geom_hline(yintercept = 0.45, linetype = "dotted") +
  geom_hline(yintercept = 0.6, linetype = "dotted")
