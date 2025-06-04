library(dplyr)
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


### Hierarchical Power ###
hier_plot_df <- read.csv(here::here("data/hier_plot_df. csv"))

hier_plot_df <- hier_plot_df %>%
  select(-scenario_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high"))


hier_power_plt <- hier_plot_df %>%
  arrange(psi) %>%
  mutate(sim_type = paste0(simulation_scenario, ": psi = ", psi, ", phi = ", phi, ", p = ", p)) %>%
  ggplot(aes(x = total_samp, y = success)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 1) +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#76BED0"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 2000, color = "darkgrey", linewidth = 1) +
  theme(legend.position = "none",
        text=element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  ylab("Percent Success") +
  xlab("Sample Size")

ggsave(here::here("figures/hier_power_plot.jpeg"), hier_power_plt)


### EMU Power ###

power_plot_df <- read.csv(here::here("data/power_plot_df.csv"))

power_plt <- power_plot_df %>%
  select(-sim_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ": psi = ", psi, ", phi = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  ggplot(aes(x = total_n, y = success)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 1) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#76BED0"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.25),
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        legend.key.spacing.y = unit(0.5, 'cm'),
        panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
        ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1))# byrow = TRUE),
         #color = guide_legend(byrow = TRUE)
         ) +
  ylab("Percent Success") +
  xlab("Sample Size")

ggsave(here::here("figures/power_plot.jpeg"), power_plt)

library(patchwork)

power_join <- power_plt + plot_spacer() +
  (plot_spacer() + hier_power_plt + plot_spacer() + plot_layout(ncol = 1, heights = c(0.5,2,0.5))) +
  plot_layout(nrow = 1, widths = c(2, 0.15, 1))

ggsave(here::here("figures/power_plot_join.jpeg"), power_join, width = 23, height = 11.5)

