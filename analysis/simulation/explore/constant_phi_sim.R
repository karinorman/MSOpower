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
  

########################################
######### Simulate using final #########
########### trend approach #############
########################################

plan(multisession, workers = 8)
sim_list_recurs <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-simulation_id, -occupancy, -name),
                                                   sim_dataset_constphi, nyear = 10, n_vis = 2,
                                                   .options=furrr_options(seed = TRUE)) %>%
                         set_names(sim_scenarios$name)) %>%
  set_names(paste0("rep", 1:simn))

map_sim <- sim_scenarios %>% 
  filter(occupancy == "high") %>% 
  rename(high_name = name) %>% 
  left_join(sim_scenarios %>% 
              filter(occupancy == "low") %>% 
              select(simulation_id, low_name = name))

true_occ_sim <- pmap_dfr(map_sim %>% select(simulation_id, high_name, low_name), function(simulation_id, high_name, low_name){

  map(1:simn, function(rep, sim_id, high, low){

    bind_rows(pluck(sim_list_recurs, rep, high)$true_occ,
              pluck(sim_list_recurs, rep, low)$true_occ
    ) %>%
      select(-site_id) %>%
      ungroup() %>%
      summarize(across(everything(), mean)) %>%
      mutate(rep = rep) %>%
    mutate(simulation_id = sim_id)
  }, sim_id = simulation_id, high = high_name, low = low_name)
})

true_occ_high <-  map_dfr(map_sim$high_name, function(high_name){
  
  map(1:simn, function(rep, high){
    
    pluck(sim_list_recurs, rep, high)$true_occ %>%
      select(-site_id) %>%
      ungroup() %>%
      summarize(across(everything(), mean)) %>%
      mutate(rep = rep) %>%
      mutate(high_name = high)
  }, high = high_name)
})


true_occ_stats_mean <- true_occ_sim %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(simulation_id, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(map_sim)

true_occ_stats_high <- true_occ_high %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(high_name, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(map_sim)

true_occ_stats_high %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = simulation_id), alpha = 0.3) +
  geom_line(aes(color = as.factor(simulation_id))) +
  facet_wrap(~sd_phi, nrow = 1) +
  theme_classic() +
  #facet_wrap(~emu, scales = "free") +
  #scale_color_discrete(name = "Sim Scenario") +
 # scale_fill_discrete(name = "Sim Scenario") +
  #geom_hline(yintercept = 0.03, linetype = "dotted") +
  #geom_hline(yintercept = 0.0225, linetype = "dotted")# +
  #geom_hline(yintercept = 0.45, linetype = "dotted") +
  geom_hline(yintercept = 0.6, linetype = "dotted") 
