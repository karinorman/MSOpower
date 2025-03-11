library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(spOccupancy)

# real world vegtypes for each emu
emu_veg <- read.csv(here::here("data/EMU_veg_types.csv")) %>%
  # let's say which we think has high or low occupancy
  mutate(occupancy = case_when(
    veg_type_landfire == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "low",
    veg_type_landfire == "Madrean Upper Montane Conifer-Oak Forest and Woodland" ~ "high",
    veg_type_landfire == "Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Aspen Forest and Woodland" ~ "low",
    veg_type_landfire == "Inter-Mountain Basins Aspen-Mixed Conifer Forest and Woodland" ~ "low",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Savanna" ~ "low",
    veg_type_landfire == "Inter-Mountain Basins Subalpine Limber-Bristlecone Pine Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Bigtooth Maple Ravine Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Riparian Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Lodgepole Pine Forest" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "low",
    veg_type_landfire == "Rocky Mountain Lodgepole Pine Forest" ~ "low",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "low",
    veg_type_landfire == "Southern Rocky Mountain Mesic Montane Mixed Conifer Forest and Woodland" ~ "high",
    .default = NA
  ))

emu_ratio <- emu_veg %>%
  group_by(UNIT, occupancy) %>%
  summarize(hex_count = sum(hex_num)) %>%
  mutate(emu = case_when(
    UNIT == "Basin & Range - East" ~ "BRE",
    UNIT == "Basin & Range - West" ~ "BRW",
    UNIT == "Colorado Plateau" ~ "CP",
    UNIT == "Southern Rocky Mountains" ~ "SRM",
    UNIT == "Upper Gila Mountains" ~ "UGM"
  )) %>%
  ungroup()


###########################################
######## Define simulation parameters #####
###########################################

## Fixed study characteristics
nyear = 10
n_sites = sum(emu_veg$hex_num)
n_vis = 2

### Parameters ###

# get dataframe of all possible scenarios
sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi1_high = c(0.2, 0.43), phi = c(.6, .8), p = c(0.4, 0.8)) %>%
  # get all possible combinations
  tidyr::expand(psi1_high, phi,p) %>%
  # and give each unique combination an ID
  mutate(sim_id = row_number()) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  mutate(sd_phi = 0.04, sd_gamma = 0.01, psi1_low = 0.03, perc_red = 0.25)

# get scenarios, one for each emu
sim_scenarios_emu <- bind_rows(sim_scenarios %>% mutate(emu = "BRE"),
                               sim_scenarios %>% mutate(emu = "BRW"),
                               sim_scenarios %>% mutate(emu = "CP"),
                               sim_scenarios %>% mutate(emu = "SRM"),
                               sim_scenarios %>% mutate(emu = "UGM")) %>%
  # get sample sizes for low and high occupancy for each emu
  left_join(emu_ratio %>%
              select(emu, occupancy, hex_count) %>%
              tidyr::pivot_wider(names_from = occupancy, values_from = hex_count) %>%
              rename(low_n = low, high_n = high)) %>%
  unite("sim_id", emu, sim_id, sep = "_") %>%
  # get the columns in the right order
  select(sim_id, psi1_low, psi1_high, phi, sd_phi, sd_gamma, p, low_n, high_n, perc_red)

###########################################
########### Generate data sets ############
###########################################

simn <- 100

#single_rep <- purrr::pmap(sim_scenarios_emu %>% select(-sim_id), sim_dataset, nyear = nyear, n_vis = 2) %>% set_names(sim_scenarios_emu$sim_id)

plan(multisession, workers = 15)
sim_list <- furrr::future_map(1:simn, ~purrr::pmap(sim_scenarios_emu %>%
                                                     select(-sim_id), sim_dataset, nyear = nyear, n_vis = 2) %>%
                                set_names(sim_scenarios_emu$sim_id),
                              seed = TRUE) %>%
  set_names(paste0("rep", 1:simn))

# reorder so top level of nested list is a sim scenario
sim_list_emu <- map(sim_scenarios_emu$sim_id, function(emu) {
  map(1:simn, ~pluck(sim_list, .x, emu)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios_emu$sim_id)

#get true occurrence for each rep and sim
true_occ <- map_dfr(sim_scenarios_emu$sim_id, function(emu){
  map_dfr(1:simn, ~pluck(sim_list_emu, emu, .x, "true_occ") %>%
            group_by(cat_var) %>%
            select(-site_id) %>%
            summarize(across(everything(), mean)) %>%
            mutate(rep = .x)) %>%
    mutate(sim_id = emu)
})

true_occ_stats <- true_occ %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(sim_id, cat_var, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  separate(sim_id, c("emu", "sim_num"), sep = "_", remove = FALSE) %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t")))

#This returns giant dataframe, hasn't been processed into encounter histories yet
# obs_occ <- map_dfr(sim_scenarios_emu$sim_id, function(emu){
#   map_dfr(1:3, ~pluck(sim_list_emu, emu, .x, "obs_occ") %>% mutate(rep = .x)) %>%
#     mutate(sim_id = emu)
# })

###########################################
########## Check Realized Trend ###########
###########################################
library(lme4)
library(broom.mixed)

true_occ_model_df <- true_occ %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t")))


model_fit_df <- true_occ_model_df %>%
  group_by(sim_id, cat_var) %>%
  nest() %>%
  # fit model for each sim_id and cat variable
  mutate(model = map(data, ~lmer(occ ~ time + (1|rep), data = .x) %>% broom.mixed::tidy())) %>%
  select(-data) %>%
  unnest(model) %>%
  # get slope and intercept for mean effect
  filter(term %in% c("time", "(Intercept)")) %>%
  select(-std.error, -statistic, -group, -effect) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`) %>%
  mutate(t10 = intercept + (time * 10),
         t1 =  intercept + time,
         percent_change = ((t10 - t1)/abs(t1))) %>%
  separate(sim_id, c("emu", "sim_num"), sep = "_", remove = FALSE) %>%
  left_join(sim_scenarios %>% mutate(sim_id = as.character(sim_id)), by = c("sim_num" = "sim_id"))



###########################################
########### Visualize True Occ ############
###########################################
library(ggplot2)

true_occ_stats %>%
  mutate(sim_num_cat = paste0(cat_var, sim_num)) %>%
  filter(emu == "BRE") %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, line = sim_num_cat, fill = as.factor(sim_num)), alpha = 0.3) +
  geom_line(aes(line = sim_num_cat, color = as.factor(sim_num))) +
  theme_classic() +
  facet_wrap(~cat_var, scales = "free") +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario")


###########################################
########### Sampling Protocol ############
###########################################

n_samp <- seq.int(100, 1000, by = 100)

# 75% of samples go in high occupancy, 25% go in low occupancy, re: recovery plan


