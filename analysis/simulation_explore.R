library(spOccupancy)

### Fixed Parameters ###
## effect_size = .25 decline
## study_duration = 10
## single deployment, two sample periods

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

### Parameters ###
psi1 = 0.43 #initial occupancy, from Wood 2019
p = c(0.4, 0.8) #detection probability, from Wood 2019
epsilon1 = .2 #extinction probability, Wood 2019, based on biology?

# derive from psi1 = gamma/(gamma + epsilon); occupancy at equilibrium
gamma1 = (psi1 * epsilon1)/(1-psi1)

# info about the sample frame
n_pot_tot <- sum(emu_veg$hex_num)

### Design scenarios ###
n_sites
n_recorders


### Simulate decline by multiplicative decline in survival

# single time period, for the whole grid
gen_dim = ceiling(sqrt(n_pot_tot))
test_data <- simOcc(J.x = gen_dim, J.y = gen_dim,
          n.rep = rep(2, gen_dim*gen_dim), beta = c(psi1, 0.2), alpha = c(p[1], 0.2))

# add trend
trend_data <- simTIntOcc()
