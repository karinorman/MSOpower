library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(ggplot2)
library(patchwork)

source(here::here("R/sim_dataset.R"))
source(here::here("R/sim_data_equilib.R"))
source(here::here("R/sim_data_constphi.R"))

# We want to compare the two approaches for simulating data by visualizing the true occupancy trend
sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = c(0.8, 0.8, .8)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  mutate(sim_type = paste0("sim_", row_number())) %>%
  tidyr::crossing(sd_phi = c(0.04, 0.1)) %>%
  mutate(sd_gamma = ifelse(sd_phi == 0.04, 0.01, 0.06)) %>%
  # and give each unique combination an ID
  mutate(sim_num = paste0("sim_", row_number()),
         occupancy = ifelse(psi == 0.03, "low", "high"),
         # Use site numbers from the Upper Gila Mountains EMU
         n_sites = ifelse(psi == 0.03, 3553, 5016)) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  mutate(perc_red = 0.25) %>%
  select(psi, phi, sd_phi, sd_gamma, p, n_sites, perc_red, sim_num, sim_type, occupancy)

simn <- 100

########################################
######### Simulate using final #########
########### trend approach #############
########################################

plan(multisession, workers = 14)
sim_list_recurs <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-sim_num, -sim_type, -occupancy),
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

# get variance in the true trend generated for each scenario
true_occ_recurs %>%
  mutate(perc_change = (t10 - t1)/t1) %>%
  select(sim_num, perc_change) %>%
  group_by(sim_num) %>%
  summarize(mean = mean(perc_change),
            variance = var(perc_change)) %>%
  mutate(bias = mean + 0.25)

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
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_type), alpha = 0.3) +
  geom_line(aes(color = sim_type)) +
  facet_wrap(~sd_phi, nrow = 1) +
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

#usethis::use_data(true_occ_stats_recurs)

##############################################
######### Simulate using equilibrium #########
############## trend approach ################
##############################################

plan(multisession, workers = 14)
sim_list_equil <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-sim_num, -sim_type, -occupancy),
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

# get variance in the true trend generated for each scenario
true_occ_equil %>%
  mutate(perc_change = (t10 - t1)/t1) %>%
  select(sim_num, perc_change) %>%
  group_by(sim_num) %>%
  summarize(mean = mean(perc_change),
            variance = var(perc_change)) %>%
  mutate(bias = mean + 0.25)

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

##############################################
######### Simulate using const phi ###########
############## trend approach ################
##############################################

plan(multisession, workers = 14)
sim_list_constphi <- map(1:simn, ~furrr::future_pmap(sim_scenarios %>% select(-sim_num, -sim_type, -occupancy, -perc_red) %>%
                                                    mutate(gamma = (0.75*psi*(1-phi))/(1 - (0.75*psi))) %>%
                                                    relocate(c(gamma, n_sites), .after = last_col()),
                                                  sim_dataset_constphi, nyear = 10, n_vis = 2,
                                                  .options=furrr_options(seed = TRUE)) %>%
                        set_names(sim_scenarios$sim_num)) %>%
  set_names(paste0("rep", 1:simn))

# reorder so top level of nested list is a sim scenario
sim_list_scenario_constphi <- map(sim_scenarios$sim_num, function(scenario) {
  map(1:simn, ~pluck(sim_list_constphi, .x, scenario)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios$sim_num)

#get true occurrence for each rep and sim
true_occ_constphi <- map_dfr(sim_scenarios$sim_num, function(scenario){
  map_dfr(1:simn, ~pluck(sim_list_scenario_constphi, scenario, .x, "true_occ") %>%
            select(-site_id) %>%
            ungroup() %>%
            summarize(across(everything(), mean)) %>%
            mutate(rep = .x)) %>%
    mutate(sim_num = scenario)
})

# get variance in the true trend generated for each scenario
true_occ_constphi %>%
  mutate(perc_change = (t10 - t1)/t1) %>%
  select(sim_num, perc_change) %>%
  group_by(sim_num) %>%
  summarize(mean = mean(perc_change),
            variance = var(perc_change)) %>%
  mutate(bias = mean + 0.25)

true_occ_stats_constphi <- true_occ_constphi %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(sim_num, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(sim_scenarios)


######## Plotting ########

pal <- c("#8A6240", "#87A96B", "#28587B", "#c9673a")

true_occ_stats_equil %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_type), alpha = 0.3) +
  geom_line(aes(color = sim_type)) +
  theme_classic() +
  facet_wrap(~sd_phi, nrow = 1) +
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
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_type), alpha = 0.3) +
  geom_line(aes(color = sim_type)) +
  theme_classic() +
  facet_wrap(~sd_phi, nrow = 1) +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted")

#usethis::use_data(true_occ_stats_equil)

####### Joined Plot #######

# This works for above
# join_df <- bind_rows(true_occ_stats_recurs %>%
#             mutate(sim_num = paste0(sim_type, "_recurs"),
#                    algo = "recursive"),
#           true_occ_stats_equil %>%
#             mutate(sim_num = paste0(sim_type, "_equil"),
#                    algo = "equilibrium")) #%>%
#   #filter(sd_phi == 0.04)

# This works for the saved datasets that produce the final figure
join_df <- bind_rows(true_occ_stats_recurs %>%
                       mutate(sim_num = paste0(sim_num, "_recurs"),
                              algo = "recursive"),
                     true_occ_stats_equil %>%
                       mutate(sim_num = paste0(sim_num, "_equil"),
                              algo = "equilibrium"))

high_occ_plot <- join_df %>%
  filter(psi == 0.6) %>%
  mutate(sim_type_name = paste0(algo, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(group = sim_num, linetype = sim_type_name, color = sim_type_name), linewidth = 0.75) +
  theme_classic() +
  #facet_wrap(~sd_phi, nrow = 1) +
  ylim(c(0.4, 0.65)) +
  ylab("Occupancy, \u03A8") +
  scale_colour_discrete("", type = rep(c("#c9673a", "#28587B"), 2)) +
  scale_linetype_manual("", values=c(2,2,1,1)) +
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


med_occ_plot <- join_df %>%
  filter(psi == 0.43) %>%
  mutate(sim_type_name = paste0(algo, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(group = sim_num, linetype = sim_type_name, color = sim_type_name), linewidth = 0.75) +
  theme_classic() +
  #facet_wrap(~sd_phi, nrow = 1) +
  ylim(c(0.3, 0.457)) +
  ylab("Occupancy, \u03A8") +
  scale_colour_discrete("", type = rep(c("#c9673a", "#28587B"), 2)) +
  scale_linetype_manual("", values=c(2,2,1,1)) +
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
        )+
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))


low_occ_plot <- join_df %>%
  filter(psi == 0.03) %>%
  mutate(sim_type_name = paste0(algo, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(group = sim_num, linetype = sim_type_name, color = sim_type_name), linewidth = 0.75) +
  theme_classic() +
  #facet_wrap(~sd_phi, nrow = 1) +
  #ylim(c(0.02, 0.032)) +
  ylab("Occupancy, \u03A8") +
  xlab("Time") +
  scale_colour_discrete("", type = rep(c("#c9673a", "#28587B"), 2)) +
  scale_linetype_manual("", values=c(2,2,1,1)) +
  # geom_hline(yintercept = 0.03, linetype = "dashed", color = "grey") +
  # geom_hline(yintercept = 0.0225, linetype = "dashed", color = "grey") +
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


ggsave(here::here("figures/occ_trend_plot.jpeg"), occ_trend_plot, width = 8, height = 11)
