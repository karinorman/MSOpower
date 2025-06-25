library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(ggplot2)
library(patchwork)

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

usethis::use_data(true_occ_stats_recurs)

##############################################
######### Simulate using equilibrium #########
############## trend approach ################
##############################################

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

pal <- c("#8A6240", "#87A96B", "#28587B", "#c9673a")

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

usethis::use_data(true_occ_stats_equil)

####### Joined Plot #######

join_df <- bind_rows(true_occ_stats_recurs %>%
            mutate(sim_num = paste0(sim_num, "_recurs"),
                   algo = "recursive"),
          true_occ_stats_equil %>%
            mutate(sim_num = paste0(sim_num, "_equil"),
                   algo = "equilibrium"))

high_occ_plot <- join_df %>%
  filter(psi == 0.6) %>%
  mutate(sim_type = paste0(algo, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(group = sim_num, linetype = sim_type, color = sim_type), linewidth = 0.75) +
  theme_classic() +
  ylim(c(0.4, 0.65)) +
  ylab("Occupancy, \u03A8") +
  #facet_wrap(~emu, scales = "free") +
  scale_colour_discrete("", type = rep(c("#c9673a", "#28587B"), 2)) +
  scale_linetype_manual("", values=c(2,2,1,1)) +
  geom_hline(yintercept = 0.45, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "grey") +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.key.size = unit(1,"cm"))+
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))


med_occ_plot <- join_df %>%
  filter(psi == 0.43) %>%
  mutate(sim_type = paste0(algo, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(group = sim_num, linetype = sim_type, color = sim_type), linewidth = 0.75) +
  theme_classic() +
  ylim(c(0.3, 0.457)) +
  ylab("Occupancy, \u03A8") +
  #facet_wrap(~emu, scales = "free") +
  scale_colour_discrete("", type = rep(c("#c9673a", "#28587B"), 2)) +
  scale_linetype_manual("", values=c(2,2,1,1)) +
  geom_hline(yintercept = 0.43, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.3225, linetype = "dashed", color = "grey") +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.key.size = unit(1,"cm"))+
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))


low_occ_plot <- join_df %>%
  filter(psi == 0.03) %>%
  mutate(sim_type = paste0(algo, ", ", "\u03C6 = ", phi)) %>%
  ggplot(aes(x = time, y = mean)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sim_num), alpha = 0.3) +
  geom_line(aes(group = sim_num, linetype = sim_type, color = sim_type), linewidth = 0.75) +
  theme_classic() +
  #ylim(c(0.02, 0.032)) +
  ylab("Occupancy, \u03A8") +
  xlab("Time") +
  #facet_wrap(~emu, scales = "free") +
  scale_colour_discrete("", type = rep(c("#c9673a", "#28587B"), 2)) +
  scale_linetype_manual("", values=c(2,2,1,1)) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.0225, linetype = "dashed", color = "grey") +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = c(0.02, 0.025, 0.03), limits = c(0.02, 0.032)) +
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.key.size = unit(1,"cm")) +
  guides(linetype = guide_legend(override.aes = list(linesize = 1)))

occ_trend_plot <- high_occ_plot + plot_spacer() + med_occ_plot + plot_spacer() + low_occ_plot +
  plot_layout(ncol = 1, guides = "collect",
              heights = c(2,.1, 2, .1, 2), axis_titles =  "collect")

ggsave(here::here("figures/occ_trend_plot.jpeg"), occ_trend_plot, width = 7.6, height = 11)
