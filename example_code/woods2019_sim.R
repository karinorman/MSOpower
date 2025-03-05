
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
gamma1 = (unlist(psi1) * epsilon1)/(1-unlist(psi1))
sd.gamma = 0.01 # woods 2019

# defining categorical variable impacting psi
# named list of number of sites in each category
# must sum to the number of sites
psi_cat_var <- list("low" = nsites/2, "high" = nsites/2)


########################################
############ Simulation ################
########################################

### Define matrix tocc (true territory occupancy state)
tocc <- matrix(rep(0, nyear*n_sites), ncol = nyear)
tocc <- as.data.frame(tocc)
names(tocc) <- paste("t", 1:nyear, sep="")

### Define the matrix obsocc that contains the observations. The first j
### columns will contain observed occupancy states for the first year on
### the 1...j visits, the j + 1 column will contain the observations of
### the first visit in the second year, and so forth.
obsocc <- matrix(rep(0, nyear*n_sites*n_vis), ncol=(nyear*n_vis),
                 dimnames = list(paste("site", 1:n_sites, sep=""),
                                 paste("Yr", rep(1:nyear, each = n_vis), "v", 1:n_vis, sep="")))
obsocc <- as.data.frame(obsocc)

# create a variable assigning sites to categorical variable
sum(unlist(cat_var_names)) == n_sites #check that it will create a vector of the correct length

cat_var <-purrr::imap(psi_cat_var, ~rep(.y, .x)) %>% unlist() %>% unname()

### simulate seasonal effects on local colonization(using an
### additive term on the logit scale drawn from N(0,sd.gamma)

### year effects on colonization without trend
logit.r.mean <- qlogis(gamma1)
year.effect.r <- rnorm((nyear-1), 0, sd.gamma)

for(i in 1:n_sites){
  r_year <- ifelse(cat_var[i] == "low", plogis(logit.r.mean["low"] + year.effect.r), plogis(logit.r.mean["high"] + year.effect.r))
}

#r_year <- plogis(logit.r.mean + year.effect.r)

### year effects on survival without trend
# logit.phi.mean <- qlogis(phi)
# year.effect.phi <- rnorm((nyear-1), 0, sd.phi)
# phi_year <- plogis(logit.phi.mean + year.effect.phi)

### year effects on survival with simulate yearly decreases
year_phi = rep(0,nyear-1)
year_phi[1] = phi

# make survival in a given year a function of previous year and the reduction
year_perc_red <- perc_red/(nyear - 1 )

for (j in 2:(nyear-1)) {
  #year_phi[j] = 1 - ((1-year_phi[j-1]) + (1-year_phi[j-1])*perc_red)
  year_phi[j] = year_phi[j-1] * (1-year_perc_red)
}

# add the noise
year.effect.phi = rep(0,nyear-1)
for (i in 1:(nyear-1)) {
  year.effect.phi[i] < - rnorm(1, year_phi[i], sd.phi)
  }


### Simulate tocc and obsocc from t = 1 to t = nyear
### First use the Control parameter values for all sites
for(i in 1:n_sites) {
  # get initial occupancy and detection histories for year 1
  # generate initial occupancy for each site
  tocc[i,1] = ifelse(cat_var[i] == "low", rbinom(1, 1, psi1$low), rbinom(1, 1, psi1$high))
  # then was it detected? Doesn't allow for false positives
  obsocc[i,1:n_vis] = rbinom(n_vis, 1, tocc[i,1]*p)

  ### the first transition (yr1 to yr2) is parameterized with the initial phi and gamma
  ### now we get the rest of the years
  for(t in 2:nyear) {
    # if the site was occupied in the previous time step
    if (tocc[i,t-1] == 1) {
      # does it persist? based on survival with noise
      tocc[i,t] <- rbinom(1,1, year.effect.phi[t-1])
    # if the site wasn't occupied
    } else if (tocc[i,t-1] == 0){
      # is it colonized, based on colonization rate with noise
      tocc[i,t] <- rbinom(1,1,r_year[t-1])
    }
    # for each visit within a year/season
    for(j in 1:n_vis) {
      # figure out whether it was detected or not
      obsocc[i,((t-1)*n_vis)+j] <- rbinom(1,1,tocc[i,t]*p)
    }
  }
}

### Output all the randomly generated pieces separately
out_list <- list(true_occ = tocc, obs_occ = obsocc, phi_survival = year.effect.phi, gamma_colonization = r_year, cat_var = cat_var)



