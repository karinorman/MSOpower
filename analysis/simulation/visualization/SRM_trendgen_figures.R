library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)

# trend gen comparison data for SRM EMU
srm_power_check <- read.csv(here::here("data/nimble_power_check_SRM.csv"))
# true occurrence for SRM comparision
srm_true_occ <- read.csv(here::here("data/true_occ_SRM.csv"))
# true occurrence by individual sim_id rather than paired
indv_srm_true_occ <- read.csv(here::here("data/true_occ_indiv_SRM.csv"))

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

# Figure comparing trend generation and modeling options for SRM EMU
nimble_srm_plt <- srm_power_check %>%
  #filter(total_n < 1501) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p),
         type_label = case_when(
           type == "constant_phi" ~ "Constant Survival",
           type == "equilibrium" ~ "Equilibrium ",
           type == "recursive" ~ "Recursive"
         )) %>%
  #filter(CI_type == .x) %>%
  ggplot(aes(x = total_n, y = ci_two_tail)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  # geom_vline(xintercept = 300) +
  # geom_vline(xintercept = 420) +
  # geom_vline(xintercept = 560) +
  # geom_vline(xintercept = 700) +
  # geom_vline(xintercept = 1000) +
  # geom_vline(xintercept = 1500) +
  facet_wrap(~type_label, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(#legend.position = "inside", legend.position.inside = c(0.9, 0.3),
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        legend.key.spacing.y = unit(0.3, 'cm'),
        panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        strip.background = element_blank(),
        strip.text = element_text(size = 13)
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1))# byrow = TRUE),
         #color = guide_legend(byrow = TRUE)
  ) +
  ylab("Percent Success") +
  xlab("Sample Size")

ggsave(here::here("figures/SRM_trendgen_power.png"), nimble_srm_plt, width = 17, height = 8)

bias_srm_plt <- srm_power_check %>%
  filter(total_n < 1501) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  #filter(CI_type == .x) %>%
  ggplot(aes(x = total_n, y = bias)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~type, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.3),
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
  ylab("mean(true trend - estimate)") +
  xlab("Sample Size")

ggsave(here::here("figures/SRM_trendgen_bias.png"), bias_srm_plt, width = 17, height = 8)

bias_srm_zoom <- srm_power_check %>%
  filter(total_n < 1501, total_n > 150) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
  #filter(CI_type == .x) %>%
  ggplot(aes(x = total_n, y = bias)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_wrap(~type, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.3),
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
  ylab("mean(true trend - estimate)") +
  xlab("Sample Size")

ggsave(here::here("figures/SRM_trendgen_bias_zoom.png"), bias_srm_zoom, width = 17, height = 8)

###########################################
####### True trend generated by sim #######
###########################################

# Look at mean trend generated with error bars for each of the generation options
trend_gen_df <- indv_srm_true_occ %>%
  pivot_longer(cols = paste0("t", 1:10), names_to = "time", values_to = "occupancy") %>%
  group_by(type, psi, phi, time) %>%
  summarize(mean_occ = mean(occupancy), mean_perc_change = mean(true_perc_change),
            lower = mean(occupancy) - qt(1- 0.05/2, (n() - 1))*sd(occupancy)/sqrt(n()),
            upper = mean(occupancy) + qt(1- 0.05/2, (n() - 1))*sd(occupancy)/sqrt(n())) %>%
  mutate(time = as.integer(unlist(stringr::str_extract_all(time, "[0-9]+"))),
         type = ifelse(type == "constant_phi", "constant survival", type))



pal <- c("#392759", "#EC9A29", "#A8201A")

high_occ_plot <- trend_gen_df %>%
  filter(psi == 0.6) %>%
  mutate(sim_type_name = paste0(type, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean_occ, linetype = sim_type_name, color = sim_type_name)) +
  #geom_ribbon(aes(group = type, ymin = lower, ymax = upper, fill = type), alpha = 0.3) +
  geom_line(linewidth = 0.75) +
  theme_classic() +
  #facet_wrap(~sd_phi, nrow = 1) +
  ylim(c(0.4, 0.65)) +
  ylab("Occupancy, \u03A8") +
  scale_colour_manual("", values = rep(pal, each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1)) +
  geom_segment(aes(x = -0.5, xend = 10.25, y = 0.45, yend = 0.45), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = -0.5, xend = 10.25, y = 0.6, yend = 0.6), linetype = "dashed", color = "grey") +
  annotate("text", x = I(1), y = 0.45, hjust = 0, label = deparse(bquote("0.75\u03A8" [i])),
           color = "darkgrey", parse = TRUE, size = 5) +
  annotate("text", x = I(1), y = 0.6, hjust = 0, label = deparse(bquote("\u03A8" [i])),
           color = "darkgrey", parse = TRUE, size = 5) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  coord_cartesian(xlim = c(0, 10), clip = 'off') +
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.key.size = unit(1,"cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        # plot.margin = margin(1,6,1,1, "cm"),
        # legend.position = c(1.35, .5)
  ) +
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))

med_occ_plot <- trend_gen_df %>%
  filter(psi == 0.43) %>%
  mutate(sim_type_name = paste0(type, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean_occ, linetype = sim_type_name, color = sim_type_name)) +
  #geom_ribbon(aes(group = type, ymin = lower, ymax = upper, fill = type), alpha = 0.3) +
  geom_line(linewidth = 0.75) +
  theme_classic() +
  #facet_wrap(~sd_phi, nrow = 1) +
  ylim(c(0.3, 0.457)) +
  ylab("Occupancy, \u03A8") +
  scale_colour_manual("", values = rep(pal, each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1)) +
  geom_segment(aes(x = -0.5, xend = 10.25, y = 0.3225, yend = 0.3225), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = -0.5, xend = 10.25, y = 0.43, yend = 0.43), linetype = "dashed", color = "grey") +
  annotate("text", x = I(1), y = 0.3225, hjust = 0, label = deparse(bquote("0.75\u03A8" [i])),
           color = "darkgrey", parse = TRUE, size = 5) +
  annotate("text", x = I(1), y = 0.43, hjust = 0, label = deparse(bquote("\u03A8" [i])),
           color = "darkgrey", parse = TRUE, size = 5) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  coord_cartesian(xlim = c(0, 10), clip = 'off') +
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.key.size = unit(1,"cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        # plot.margin = margin(1,6,1,1, "cm"),
        # legend.position = c(1.35, .5)
  ) +
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))

low_occ_plot <- trend_gen_df %>%
  filter(psi == 0.03) %>%
  mutate(sim_type_name = paste0(type, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean_occ, linetype = sim_type_name, color = sim_type_name)) +
  #geom_ribbon(aes(group = type, ymin = lower, ymax = upper, fill = type), alpha = 0.3) +
  geom_line(linewidth = 0.75) +
  theme_classic() +
  #facet_wrap(~sd_phi, nrow = 1) +
  ylim(c(0.3, 0.457)) +
  ylab("Occupancy, \u03A8") +
  scale_colour_manual("", values = rep(pal, each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1)) +
  geom_segment(aes(x = -0.5, xend = 10.25, y = 0.03, yend = 0.03), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = -0.5, xend = 10.25, y = 0.0225, yend = 0.0225), linetype = "dashed", color = "grey") +
  annotate("text", x = I(1), y = 0.0225, hjust = 0, label = deparse(bquote("0.75\u03A8" [i])),
           color = "darkgrey", parse = TRUE, size = 5) +
  annotate("text", x = I(1), y = 0.03, hjust = 0, label = deparse(bquote("\u03A8" [i])),
           color = "darkgrey", parse = TRUE, size = 5) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = c(0.02, 0.025, 0.03), limits = c(0.02, 0.032)) +
  coord_cartesian(xlim = c(0, 10), clip = 'off') +
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.key.size = unit(1,"cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        # plot.margin = margin(1,6,1,1, "cm"),
        # legend.position = c(1.35, .5)
  ) +
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))

# it fills rows first!!! this is very annoying!!
occ_trend_plot <- high_occ_plot + plot_spacer() + plot_spacer() +
  plot_spacer() + plot_spacer() + plot_spacer() +
  med_occ_plot + plot_spacer() + guide_area() +
  plot_spacer() + plot_spacer() + plot_spacer() +
  low_occ_plot + plot_spacer() + plot_spacer() +
  plot_layout(ncol = 3, nrow = 5, heights = c(2,.1, 2, .1, 2), widths = c(1, 0.1, 1), axis_titles = "collect", guides = "collect")

ggsave(here::here("figures/occ_trend_plot_srm.jpg"), occ_trend_plot, width = 8, height = 11)

###########################################
####### comparison of power checks  #######
###########################################

power_crit_plt <- srm_power_check %>%
  #filter(type == "recursive") %>%
  pivot_longer(c(ci_two_tail, ci_any_decline, ci_left_tail), names_to = "ci_type", values_to = "power") %>%
  #filter(total_n < 1501) %>%
  left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
  mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p),
         type_label = case_when(
           ci_type == "ci_two_tail" ~ "Two Tail",
           ci_type == "ci_any_decline" ~ "Any Decline",
           ci_type == "ci_left_tail" ~ "Left Tail",
         ),
         gen_label = case_when(
           type == "constant_phi" ~ "Constant Survival",
           type == "equilibrium" ~ "Equilibrium",
           type == "recursive" ~ "Recursive")) %>%
  #filter(CI_type == .x) %>%
  ggplot(aes(x = total_n, y = power)) +
  geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
  facet_grid(type_label ~ gen_label, scales = "free_x") +
  theme_classic() +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
  geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(#legend.position = "inside", legend.position.inside = c(0.9, 0.3),
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        legend.key.spacing.y = unit(0.3, 'cm'),
        panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        strip.background = element_blank(),
        strip.text = element_text(size = 15)
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1))# byrow = TRUE),
         #color = guide_legend(byrow = TRUE)
  ) +
  ylab("Percent Success") +
  xlab("Sample Size")

ggsave(here::here("figures/SRM_trendgen_power.png"), power_crit_plt, width = 15, height = 12)


###########################################
######## comparison of thresholds  ########
###########################################

threshold_df <- read.csv(here::here("data/srm_trendgen_power_thresholds.csv")) %>%
  select(-line_id) %>%
  left_join(sim_scenarios_table %>% select(simulation_scenario, phi, p, psi = psi_high))

bracket_ends <- threshold_df %>%
  mutate(threshold = round(threshold)) %>%
  group_by(simulation_scenario) %>%
  mutate(difference = max(threshold) - min(threshold),
         type_rank = case_when(
           threshold == max(threshold) ~ "max",
           threshold == min(threshold) ~ "min"
         ),
         xval = simulation_scenario + .25,
         halfway = min(threshold) + (difference/2)) %>%
  filter(!is.na(type_rank)) %>%
  ungroup()

bracket_df <- bracket_ends %>%
  select(simulation_scenario, threshold, type_rank, xval, halfway, difference) %>%
  pivot_wider(names_from = "type_rank", values_from = "threshold")

threshold_plot <-
  threshold_df %>%
  mutate(shape_var = paste0(type, "_", simulation_scenario)) %>%
  ggplot() +
  geom_point(data = bracket_ends, aes(x = xval, y = threshold), shape = 95, size = 5) +
    geom_segment(data = bracket_df,
                 aes(x = xval, y = max, xend = xval, yend = min)) +
    geom_text(data = bracket_df, aes(x = (xval + .2), y = halfway, label = difference)) +
  geom_point(aes(x = simulation_scenario, y = threshold, color = as.factor(simulation_scenario), shape = shape_var), size = 4) +
  scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
  scale_shape_manual("", values= rep(c(1, 19, 2, 17, 0, 15), 4)) +
  ylab("Threshold Sample Size") +
  xlab("Scenario") +
  theme_classic() +
  theme(legend.position = "none",
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        #legend.key.spacing.y = unit(0.5, 'cm'),
        panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  scale_x_continuous(breaks = 1:8)

legend_plot <- threshold_df %>%
  mutate(type_label = case_when(
    type == "constant_phi" ~ "Constant Survival",
    type == "equilibrium" ~ "Equilibrium",
    type == "recursive" ~ "Recursive"
  )) %>%
  #mutate(shape_var = paste0(type, "_", line_id)) %>%
  ggplot() +
  geom_point(aes(x = as.factor(simulation_scenario), y = threshold, shape = type_label), color = "grey", size = 3) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text=element_text(size=12),
        legend.key.width = unit(1,"cm"),
        text=element_text(size=14),
        legend.key.spacing.y = unit(0.5, 'cm'),
        #panel.spacing = unit(30, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

legend <- get_legend(legend_plot)

bottom_row <- plot_grid(threshold_plot, legend, nrow = 1, rel_widths = c(1, .2))

trend_eval_plt <- plot_grid(nimble_srm_plt, bottom_row, nrow = 2, labels = "AUTO")

save_plot(here::here("figures/trend_eval.png"), trend_eval_plt, nrow = 2, ncol = 2, bg = 'white', base_asp = 1.5, base_height = 5)


# # Power check figure for single trend generation option for all EMU's
# nimble_emu_plt <- nimble_power_check %>%
#   #filter(total_n < 1500) %>%
# left_join(sim_scenarios_table, by = c("phi", "p", "psi" = "psi_high")) %>%
#   mutate(sim_type = paste0(simulation_scenario, ":  \u03A8 = ", psi, ", \u03C6 = ", phi, ", p = ", p)) %>%
#   mutate(unit = case_when(
#     emu == "BRE" ~ "Basin & Range - East",
#     emu == "BRW" ~ "Basin & Range - West",
#     emu == "CP" ~ "Colorado Plateau",
#     emu == "SRM" ~ "Southern Rocky Mountains",
#     emu == "UGM" ~ "Upper Gila Mountains"
#   )) %>%
#   #filter(CI_type == .x) %>%
#   ggplot(aes(x = total_n, y = ci_two_tail)) +
#   geom_line(aes(color = sim_type, linetype = sim_type), linewidth = 0.75) +
#   facet_wrap(~unit, scales = "free_x") +
#   theme_classic() +
#   scale_colour_discrete("", type = rep(c("#8A6240", "#87A96B", "#28587B", "#c9673a"), each = 2)) +
#   scale_linetype_manual("", values=c(2,1,2,1,2,1,2,1)) +
#   geom_hline(yintercept = 0.9, color = "darkgrey", linetype = "dashed", linewidth = 1) +
#   theme(legend.position = "inside", legend.position.inside = c(0.85, 0.25),
#         legend.text=element_text(size=12),
#         legend.key.width = unit(1,"cm"),
#         text=element_text(size=14),
#         legend.key.spacing.y = unit(0.5, 'cm'),
#         panel.spacing = unit(30, "pt"),
#         axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
#         axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
#   ) +
#   guides(linetype = guide_legend(override.aes = list(linewidth = 1))# byrow = TRUE),
#          #color = guide_legend(byrow = TRUE)
#   ) +
#   ylab("Percent Success") +
#   xlab("Sample Size")
#
# ggsave(here::here("figures/emu_power_nimble.png"), nimble_emu_plt, height = 12, width = 15)

