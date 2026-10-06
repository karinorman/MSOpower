library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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
hier_plot_df <- read.csv(here::here("data/nimble_power_check_hier.csv")) %>%
  mutate(scenario_id = chunk_num)

hier_plot_df <- hier_plot_df %>%
  select(-scenario_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high"))

# ci_plotting_df <- tibble(ci = c(0.95, 0.90, 0.95, 0.90),
#                          y = c("ci_check", "ci_check_left_tail", "ci_check_any_decline", "post_check_any_decline"))

hier_power_plt <- #purrr::pmap(ci_plotting_df, ~hier_plot_df %>%
                                #filter(CI_type == .x) %>%
  hier_plot_df %>%
  arrange(psi) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  ggplot(aes(x = total_n, y = ci_two_tail)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  #geom_vline(xintercept = 2000, color = "darkgrey", linewidth = 1) +
  theme(legend.position = "none",
        text=element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  ylab("Percent Success") +
  xlab("Sample Size")
#)

#ggsave(here::here("figures/hier_power_plot.jpeg"), hier_power_plt)

hier_power_simp <- hier_plot_df %>%
  arrange(psi) %>%
  filter(phi == 0.6) %>%
  mutate(sim_type = paste0(":  \u03A8 = ", psi, ", p = ", p)) %>%
  ggplot(aes(x = total_n, y = ci_two_tail)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 1.25) +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  #geom_vline(xintercept = 2000, color = "darkgrey", linewidth = 1) +
  theme(#legend.position = "none",
        text=element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent', color = NA)) +
  ylab("Percent Success") +
  xlab("Sample Size")
#)

ggsave(here::here("figures/hier_power_plot_simplified.png"), hier_power_simp, bg = "transparent", height = 8, width = 9)

### EMU Power ###


power_eval <- read.csv(here::here("data/nimble_power_check_emu.csv"))

power_plt <- power_eval %>%
  #filter(total_n < 1200) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  ggplot(aes(x = total_n, y = ci_two_tail)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  # geom_vline(xintercept = 300) +
  # geom_vline(xintercept = 420) +
  # geom_vline(xintercept = 560) +
  # geom_vline(xintercept = 700) +
  # geom_vline(xintercept = 930) +
  # geom_vline(xintercept = 1000) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.25),
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        legend.key.spacing.y = unit(0.3, 'cm'),
        panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent', color = NA,),
        strip.background =element_rect(fill="transparent")
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1))# byrow = TRUE),
         #color = guide_legend(byrow = TRUE)
  ) +
  ylab("Percent Success") +
  xlab("Sample Size")

ggsave(here::here("figures/power_plot.png"), power_plt, bg = "transparent")

power_join <- power_plt + plot_spacer() +
  (plot_spacer() + hier_power_plt + plot_spacer() + plot_layout(ncol = 1, heights = c(0.5,2,0.5))) +
  plot_layout(nrow = 1, widths = c(2, 0.15, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(here::here("figures/power_plot_join.jpeg"), power_join, width = 18, height = 9, bg = "white")

# # get plots for different kinds of power checks
# ggsave(here::here("figures/power_check_twotail.png"), power_plt[[1]], width = 16.5, height = 9.87)
# ggsave(here::here("figures/power_check_lefttail.png"), power_plt[[2]], width = 16.5, height = 9.87)
# ggsave(here::here("figures/power_check_negative.png"), power_plt[[3]], width = 16.5, height = 9.87)


################################
#### Separate by detection #####
################################

hier_power_p4 <- hier_plot_df %>%
  arrange(psi) %>%
  filter(p == 0.4) %>%
  mutate(sim_type = paste0("\u03A8 = ", psi, ", \u03C6 = ", phi)) %>%
  ggplot(aes(x = total_samp, y = success)) +
  geom_line(aes(color = sim_type), linewidth = 0.75) +
  theme_classic() +
  scale_colour_discrete("", type = c("#8A6240", "#87A96B", "#28587B", "#c9673a")) +
  #scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 2000, color = "darkgrey", linewidth = 1) +
  theme(#legend.position = "none",
        text=element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  ylab("Percent Success") +
  xlab("Sample Size")

hier_power_p8 <- hier_plot_df %>%
  arrange(psi) %>%
  filter(p == 0.8) %>%
  mutate(sim_type = paste0("\u03A8 = ", psi, ", \u03C6 = ", phi)) %>%
  ggplot(aes(x = total_samp, y = success)) +
  geom_line(aes(color = sim_type), linewidth = 0.75) +
  theme_classic() +
  scale_colour_discrete("", type = c("#8A6240", "#87A96B", "#28587B", "#c9673a")) +
  #scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 2000, color = "darkgrey", linewidth = 1) +
  theme(#legend.position = "none",
        text=element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  ylab("Percent Success") +
  xlab("Sample Size")

power_p4 <- power_plot_df %>%
  select(-sim_id) %>%
  filter(p == 0.4) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ": \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  ggplot(aes(x = total_n, y = success)) +
  geom_line(aes(color = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = c("#8A6240", "#87A96B", "#28587B", "#c9673a")) +
  #scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(legend.position = "none",
        #legend.position = "inside", legend.position.inside = c(0.85, 0.25),
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        legend.key.spacing.y = unit(0.5, 'cm'),
        panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1))# byrow = TRUE),
         #color = guide_legend(byrow = TRUE)
  ) +
  ylab("Percent Success")

power_p8 <- power_plot_df %>%
  select(-sim_id) %>%
  filter(p == 0.8) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ": \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  ggplot(aes(x = total_n, y = success)) +
  geom_line(aes(color = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = c("#8A6240", "#87A96B", "#28587B", "#c9673a")) +
  #scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(legend.position = "none",
        #legend.position = "inside", legend.position.inside = c(0.85, 0.25),
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


(power_p4 + plot_spacer() +
  (plot_spacer() + hier_power_p4 + plot_spacer() + plot_layout(ncol = 1, heights = c(0.5,2,0.5))) +
  plot_layout(nrow = 1, widths = c(2, 0.15, 1))) /
  plot_spacer() /
  (power_p8 + plot_spacer() +
     (plot_spacer() + hier_power_p8 + plot_spacer() + plot_layout(ncol = 1, heights = c(0.5,2,0.5))) +
     plot_layout(nrow = 1, widths = c(2, 0.15, 1))) +
  plot_layout(heights = c(2,.3, 2))

power_by_detection <- (power_p4 + plot_spacer() + power_p8 + plot_layout(axis_titles = "collect", heights = c(2,.3, 2))) /
  plot_spacer() /
  (plot_spacer() + hier_power_p4 + guide_area() + hier_power_p8 + plot_spacer() + plot_layout(guides = "collect", ncol = 1, heights = c(0.15, 1, 0.75, 1, 0.15))) +
  plot_layout(ncol = 3, widths = c(2, 0.15, 1))

ggsave(here::here("figures/power_plot_detection.jpeg"), power_by_detection, width = 20, height = 15)



##### Bias #####

## Absolute Bias
bias_plt <- power_eval %>%
  #filter(total_n > 250) %>%
  select(-sim_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  filter(CI_type == 0.90) %>%
  ggplot(aes(x = total_n, y = bias)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 1) +
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
  ylab("Mean absolute bias") +
  xlab("Sample Size")

ggsave(here::here("figures/bias_plot.png"), bias_plt, width = 19.6, height = 12)

## Relative Bias
rel_bias_plt <- power_eval %>%
  #filter(total_n > 250) %>%
  select(-sim_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  filter(CI_type == 0.90) %>%
  ggplot(aes(x = total_n, y = relative_bias)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 1) +
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
  ylab("Mean absolute bias") +
  xlab("Sample Size")

ggsave(here::here("figures/relative_bias_plot.png"), rel_bias_plt, width = 19.6, height = 12)

## For a two-tailed, alpha = 0.05 test, what percentage of reps is the true trend less than the estimated trend
### This plot shows that we're slightly more likely to underestimate the trend (true trend is less than estimated trend) across all replicates
power_eval %>%
  select(-sim_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  filter(CI_type == 0.90) %>%
  ggplot(aes(x = total_n, y = percent_trend_lower)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.5, color = "darkgrey", linetype = "dashed", linewidth = 1) +
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
  ylab("Percent of Reps where True Trend is less than estimated trend") +
  xlab("Sample Size")



## For a two-tailed, alpha = 0.05 test, what percentage of reps is the true trend less than the estimated trend for only reps where
## estimated trend is not included in CI
### This plot shows that for reps not in the confidence interval, we're way more likely to overestimate the trend (true trend is greater than estimated trend)
bias_ci_in_plt <- power_eval %>%
  select(-sim_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  filter(CI_type == 0.90) %>%
  ggplot(aes(x = total_n, y = percent_exclude_trend_lower)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.5, color = "darkgrey", linetype = "dashed", linewidth = 1) +
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
  ylab("Percent of Reps where True Trend not in CI and is less than estimated trend") +
  xlab("Sample Size")

#ggsave(here::here("figures/bias_plot.png"), bias_ci_in_plt, width = 19.6, height = 12)

## Precision figure ##

precision_plt <- power_eval %>%
  select(-sim_id) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  mutate(unit = case_when(
    emu == "BRE" ~ "Basin & Range - East",
    emu == "BRW" ~ "Basin & Range - West",
    emu == "CP" ~ "Colorado Plateau",
    emu == "SRM" ~ "Southern Rocky Mountains",
    emu == "UGM" ~ "Upper Gila Mountains"
  )) %>%
  filter(CI_type == 0.95) %>%
  ggplot(aes(x = total_n, y = width)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~unit, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  #geom_hline(yintercept = 0.5, color = "darkgrey", linetype = "dashed", linewidth = 1) +
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
  ylab("Average CI Width") +
  xlab("Sample Size")

ggsave(here::here("figures/precision_plot.png"), precision_plt, width = 16.2, height = 12)
