library(dplyr)
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

emu_ratio <- emu_veg %>%
  group_by(UNIT, occupancy) %>%
  summarize(hex_count = sum(hex_num))

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
          n.rep = rep(2, gen_dim*gen_dim),
          # categorical variable controlling occupancy, so intercept is lower reference category
          # and slope is reference category + amount to bring up to high occupancy
          beta = c(psi1 - 0.2, psi1),
          # intercept with no covariates
          alpha = c(p[1], 0))

# Detection-nondetection data
y <- test_data$y
# Occurrence design matrix for fixed effects
X <- test_data$X
# Detection design matrix for fixed effets
X.p <- test_data$X.p
# Occurrence values
psi <- test_data$psi
# Spatial coordinates
coords <- test_data$coords

# Package all data into a list
# Occurrence covariates consists of the fixed effects and random effect
occ.covs <- as.matrix(X[, 2])
colnames(occ.covs) <- c('occ.cov.1')

# Detection covariates consists of the fixed effects and random effect
det.covs <- list(det.cov.1 = X.p[, , 2])
                 #det.factor.1 = X.p.re[, , 1])
# Package into a list for spOccupancy
data.list <- list(y = y,
                  occ.covs = occ.covs,
                  det.covs = det.covs,
                  coords = coords)

inits <- list(alpha = 0, beta = 0, sigma.sq.psi = 0.5,
              sigma.sq.p = 0.5, z = apply(test_data$y, 1, max, na.rm = TRUE),
              sigma.sq = 1, phi = 3 / 0.5)
out.full <- PGOcc(occ.formula = ~ occ.cov.1,
                    det.formula = ~ det.cov.1,
                    data = data.list,
                    inits = inits,
                    n.samples = 5000,
                    n.report = 100,
                    n.burn = 1000,
                    n.chains = 1)


# add trend
trend_data <- simTIntOcc()
