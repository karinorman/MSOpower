####################################################################################################################
## Power Analysis Goal: Owl occupancy rates must show a stable or increasing trend after 10 years of monitoring.
## The study design to verify this criterion must have a power of 90% (Type II error rate β = 0.10) to detect a
## 25% decline in occupancy rate over the 10-year period with a Type I error rate (α) of 0.10.
####################################################################################################################

library(dplyr)
library(tidyr)

# real world vegtypes for each emu
emu_veg <- read.csv(here::here("data/EMU_veg_types.csv")) %>%
  # let's say which we think has high or low occupancy
  mutate(occupancy = case_when(
    veg_type_landfire == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "high",
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
    veg_type_landfire == "Madrean Pinyon-Juniper Woodland" ~ "low",
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
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = c(0.4, 0.8, .8)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  # and give each unique combination an ID
  mutate(occupancy = ifelse(psi == 0.03, "low", "high")) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  mutate(sd_phi = 0.04, sd_gamma = 0.01, perc_red = 0.25)

# get scenarios, one for each emu
sim_scenarios_emu <- bind_rows(sim_scenarios %>% mutate(emu = "BRE"),
                               sim_scenarios %>% mutate(emu = "BRW"),
                               sim_scenarios %>% mutate(emu = "CP"),
                               sim_scenarios %>% mutate(emu = "SRM"),
                               sim_scenarios %>% mutate(emu = "UGM")) %>%
  # get sample sizes for low and high occupancy for each emu
  left_join(emu_ratio %>%
              select(emu, occupancy)) %>%
  # tidyr::pivot_wider(names_from = occupancy, values_from = hex_count) %>%
  # rename(low_n = low, high_n = high)) %>%
  # get the columns in the right order
  select(emu, psi, phi, sd_phi, sd_gamma, p, perc_red) 

samp_size_df <- data.frame(total_n = c(20, 40, 60, 100, 160, 200, 300, 420, 560, 700, 1000, 1500, 2000)) %>%
  mutate(high_n = 0.75*total_n, low_n = 0.25*total_n) %>%
  group_by(total_n) %>%
  slice(rep(row_number(), n_distinct(emu_ratio$emu))) %>%
  mutate(emu = unique(emu_ratio$emu)) %>%
  # make sure we have enough hexes of each category to do that sample size
  left_join(emu_ratio %>% filter(occupancy == "high") %>% select(emu, high_hex_count = hex_count) %>% distinct()) %>%
  left_join(emu_ratio %>% filter(occupancy == "low") %>% select(emu, low_hex_count = hex_count) %>% distinct()) %>%
  filter(high_hex_count > high_n, low_hex_count > low_n) %>%
  # add back in highest sample size allowed by available high occupancy hex area for BRE and BRW %>%
  bind_rows(data.frame(emu = c("BRE", "BRW"), 
             high_n = c(744, 370)
             )) %>%
  ungroup() %>%
  mutate(low_n = ifelse(is.na(low_n), round(high_n * .25), low_n),
         total_n = ifelse(is.na(total_n), low_n + high_n, total_n)) %>%
  arrange(emu) %>%
  fill(high_hex_count, low_hex_count, .direction = "down")


emu_scenarios_samp <- sim_scenarios_emu %>%
  left_join(samp_size_df) %>%
  group_by(emu) %>%
  mutate(sim_num = row_number()) %>%
  unite("sim_id", emu, sim_num, remove = FALSE) %>%
  select(-sim_num)

# generate master scenario table with random seeds
# set seed for random seed generator
set.seed(524876)
sim_map_names <- emu_scenarios_samp %>%
  filter(psi != 0.03) %>%
  rename(high_name = sim_id) %>%
  left_join(emu_scenarios_samp %>%
              select(low_name = sim_id, psi, phi, p, low_n, total_n, emu) %>%
              filter(psi == 0.03) %>%
              select(-psi)) %>%
  mutate(low_psi = 0.03,
         seed = 1 + floor(runif(n()) * 100000)) %>%
  group_by(emu, total_n) %>%
  mutate(chunk_num = cur_group_id()) %>%
  ungroup() %>%
  select(-emu) %>%
  group_by(low_name, high_name) %>%
  mutate(scenario_id = cur_group_id()) %>%
  ungroup()

usethis::use_data(sim_map_names)
