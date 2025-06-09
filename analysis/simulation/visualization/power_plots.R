library(dplyr)
library(tidyr)
library(ggplot2)

## Fixed study characteristics
nyear = 10
#n_sites = sum(emu_veg$hex_num)
n_vis = 2

### Parameters ###

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

hier_plot_df <- read.csv(here::here("data/hier_plot_df.csv"))

hier_plot_df <- hier_plot_df %>%
  select(-scenario_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high"))

power_plt <- power_plot_df %>%
  mutate(sim_type = paste0("Psi = ", psi, ", Phi = ", phi)) %>%
  ggplot(aes(x = total_samp, y = success)) +
  geom_line(aes(color = as.factor(sim_type), linetype = forcats::fct_rev(as.factor(p)))) +
  #facet_wrap(~emu, scales = "free") +
  theme_classic() +
  scale_color_discrete(name = "Scenario") +
  scale_linetype_discrete(name = "Detection") +
  geom_hline(yintercept = 0.9, color = "darkgrey")#, linetype = "dotted")

hier_plot_df %>%
  arrange(psi) %>%
  mutate(sim_type = paste0(simulation_scenario, ": psi = ", psi, ", phi = ", phi, ", p = ", p)) %>%
  ggplot(aes(x = total_samp, y = success)) +
  geom_line(aes(color = sim_type, linetype = sim_type)) +
  theme_classic() +
  scale_colour_discrete("") +
  scale_linetype_manual("", values=c(1,2,1,2,1,2)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dotted")


ggsave(here::here("figures/hier_power_plot.jpeg"), power_plt)
