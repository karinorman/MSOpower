library(dplyr)
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
psi1 = list("low" = 0.03, "high" = 0.43) #initial occupancy, from Wood 2019
p = c(0.4, 0.8) #detection probability, from Wood 2019
epsilon1 = .2 #extinction probability, Wood 2019, based on biology?
phi = 1 - epsilon1 # survival
sd.phi = 0.04 # also 0.01 in Woods 2019

perc_red = 0.25 # 25% decline in occupancy

# derive from psi1 = gamma/(gamma + epsilon); occupancy at equilibrium
gamma = (unlist(psi1) * epsilon1)/(1-unlist(psi1))
sd.gamma = 0.01 # woods 2019

# defining categorical variable impacting psi
# named list of number of sites in each category
# must sum to the number of sites
psi_cat_var <- list("low" = n_sites/2, "high" = n_sites/2)

# we need all of these parameters for each scenario, the rest are fixed
# perc_red, psi1, phi, sd.phi, gamma1, sd.gamma, p, cat_var_list
sim_scenarios <- data.frame(psi1_high = c(0.2, 0.43), phi = c(.6, .8), p = c(0.4, 0.8)) %>%
  tidyr::expand(psi1_high, phi,p) %>%
  mutate(sim_id = row_number()) %>%
  # these are the same for all scenarios right now
  mutate(sd_phi = 0.04, sd_gamma = 0.01, psi1_low = 0.03, perc_red = 0.25)

sim_scenarios_emu <- bind_rows(sim_scenarios %>% mutate(emu = "BRE"),
                               sim_scenarios %>% mutate(emu = "BRW"),
                               sim_scenarios %>% mutate(emu = "CP"),
                               sim_scenarios %>% mutate(emu = "SRM"),
                               sim_scenarios %>% mutate(emu = "UGM")) %>%
  left_join(emu_ratio %>%
              select(emu, occupancy, hex_count) %>%
              tidyr::pivot_wider(names_from = occupancy, values_from = hex_count) %>%
              rename(low_n = low, high_n = high)) %>%
  select(sim_id, emu, psi1_low, psi1_high, phi, sd_phi, sd_gamma, p, low_n, high_n, perc_red)

###########################################
########### Generate data sets ############
###########################################

single_rep <- purrr::pmap(sim_scenarios_emu[1,] %>% select(-sim_id, -emu), sim_dataset, nyear = nyear, n_vis = 2)
