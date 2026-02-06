library(dplyr)

# get dataframe of all possible scenarios
sim_scenarios_table <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = c(0.4, 0.8, .8)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  # and give each unique combination an ID
  mutate(sim_num = row_number(),
         occupancy = ifelse(psi == 0.03, "low", "high")) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  #mutate(sd_phi = 0.04, sd_gamma = 0.01) %>%
  select(-sim_num) %>%
  #rowwise() %>%
  mutate(gamma = round((psi * (1-phi))/(1-psi), 2)) %>%
  #ungroup() %>%
  #group_by(phi, p, sd_phi, sd_gamma) %>%
  pivot_wider(names_from = "occupancy", values_from = c("gamma", "psi"), values_fn = list) %>%
  unnest(c(gamma_high, psi_high)) %>%
  arrange(psi_high, phi) %>%
  mutate(simulation_scenario = row_number()) %>%
  select(simulation_scenario,everything())

# Collate threshold dataframes to output for table 2 in the manuscript

emu_threshold <- read.csv(here::here("data/emu_power_thresholds.csv")) %>%
  select(emu, psi, phi, p, threshold) %>%
  pivot_wider(names_from = "emu", values_from = "threshold")

hier_threshold <- read.csv(here::here("data/hier_power_thresholds.csv")) %>%
  select(psi, phi, p, hierarchical = threshold)

trend_threshold <- read.csv(here::here("data/srm_trendgen_power_thresholds.csv")) %>%
  select(type, psi, phi, p, threshold) %>%
  pivot_wider(names_from = "type", values_from = "threshold")

thresholds <- emu_threshold %>%
  left_join(hier_threshold) %>%
  left_join(trend_threshold) %>%
  mutate(across(-c(psi, phi, p), round)) %>%
  left_join(sim_scenarios_table %>% select(simulation_scenario, phi, p, psi = psi_high)) %>%
  ungroup() %>%
  arrange(simulation_scenario)

readr::write_csv(thresholds, here::here("data/thresholds.csv"))

thresholds_long <- thresholds %>%
  pivot_longer(cols = -c(psi, phi, p, simulation_scenario), names_to = "simulation_type", values_to = "sample_size") %>%
  select(simulation_scenario, everything())

readr::write_csv(thresholds_long, here::here("data/thresholds_long.csv"))
