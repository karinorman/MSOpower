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
  summarize(hex_count = sum(hex_num))


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
sim_scenarios <- data.frame(psi1_low = rep(0.03, 2), psi1_high = rep(0.43, 2), phi = rep(0.8, 2),
                               sd_phi = rep(0.04, 2), sd_gamma = rep(0.01, 2), p = c(0.4, 0.8),
                            low_n = rep(15602, 2), high_n = rep(15602, 2), perc_red = rep(0.25, 2))

###########################################
########### Generate data sets ############
###########################################

single_rep <- purrr::pmap(sim_scenarios, sim_dataset, nyear = nyear, n_vis = 2)
