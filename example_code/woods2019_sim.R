
## Fixed study characteristics
nyear = 10
n_sites = sum(emu_veg$hex_num)
n_vis = 2

### Parameters ###
psi1 = list("low" = 0.03, "high" = 0.43) #initial occupancy, from Wood 2019
p = c(0.4, 0.8) #detection probability, from Wood 2019
epsilon1 = .2 #extinction probability, Wood 2019, based on biology?
phi = 1 - epsilon1 # survival
sd_phi = 0.04 # also 0.01 in Woods 2019

perc_red = 0.25 # 25% decline in occupancy

# derive from psi1 = gamma/(gamma + epsilon); occupancy at equilibrium
gamma = (unlist(psi1) * epsilon1)/(1-unlist(psi1))
sd_gamma = 0.01 # woods 2019

# defining categorical variable impacting psi
# named list of number of sites in each category
# must sum to the number of sites
psi_cat_var <- list("low" = n_sites/2, "high" = n_sites/2)


########################################
############ Simulation ################
########################################

woods_sim_dataset <- function(psi1_low, psi1_high, phi, sd_phi, sd_gamma, p, low_n, high_n, perc_red, nyear, n_vis){

  ### psi1 = initial occupancy
  ### phi = local survival in year 1
  ### sd_phi = sd of random normal variable introducing a Season effect
  ### on local survival
  ### sd_gamma = sd of random normal variable introducing a Season effect
  ### on colonization
  ### p = probability of detection
  ### perc_red = total reduction across simulated time in survival
  ### nyear = number of years or study seasons
  ### n_vis = number of visits per site per year (or season)

  # get parameters that have multiple values in correct format
  cat_n <- setNames(c(low_n, high_n), c("low", "high"))
  psi1 <- setNames(c(psi1_low, psi1_high), c("low", "high"))
  # derive from psi1 = gamma/(gamma + epsilon); occupancy at equilibrium
  gamma <- (psi1 * (1-phi))/(1-psi1)
  nsites <- sum(cat_n)

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
  cat_var <-purrr::imap(cat_n, ~rep(.y, .x)) %>% unlist() %>% unname()

  ### simulate seasonal effects on local colonization(using an
  ### additive term on the logit scale drawn from N(0,sd_gamma)

  ### year effects on colonization without trend
  logit.gamma.mean <- qlogis(gamma)
  year.effect.gamma <- rnorm((nyear-1), 0, sd_gamma)

  gamma_year = list()
  gamma_year[[1]] = gamma

  for(j in 2:(nyear-1)){
    gamma_year[[j]] <- plogis(logit.gamma.mean + year.effect.gamma[j-1])
  }

  #r_year <- plogis(logit.r.mean + year.effect.r)

  ### year effects on survival without trend
  # logit.phi.mean <- qlogis(phi)
  # year.effect.phi <- rnorm((nyear-1), 0, sd_phi)
  # phi_year <- plogis(logit.phi.mean + year.effect.phi)

  ### year effects on survival with simulate yearly decreases
  year_phi = rep(0,nyear-1)
  year_phi[1] = phi

  # make survival in a given year a function of previous year and the reduction
  #year_perc_red <- 1 - log(exp(1 -perc_red)/(nyear - 1))

  #if we want occupancy to decline by perc_red, how much should survival decline in that time period
  #survival_perc_red <- (((1-perc_red)*unlist(psi1["high"])*(1-phi))/perc_red)*unlist(gamma["high"])
  # and how much should survival decline yearly
  year_perc_red <- 1 - (exp(log(1 - perc_red)/(nyear)))


  for (j in 2:(nyear-1)) {
    #year_phi[j] = 1 - ((1-year_phi[j-1]) + (1-year_phi[j-1])*perc_red)
    year_phi[j] = year_phi[j-1] * (1 - year_perc_red)
  }

  # add the noise
  year.effect.phi = rep(0,nyear-1)
  for (i in 1:(nyear-1)) {
    year.effect.phi[i] <- rnorm(1, year_phi[i], sd_phi)
  }


  ### Simulate tocc and obsocc from t = 1 to t = nyear
  ### First use the Control parameter values for all sites
  for(i in 1:n_sites) {
    # get initial occupancy and detection histories for year 1
    # generate initial occupancy for each site
    tocc[i,1] = ifelse(cat_var[i] == "low", rbinom(1, 1, psi1["low"]), rbinom(1, 1, psi1["high"]))
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
        tocc[i,t] <- ifelse(cat_var[i] == "low", rbinom(1,1,gamma_year[[t-1]]["low"]), rbinom(1,1,gamma_year[[t-1]]["high"]))
      }
      # for each visit within a year/season
      for(j in 1:n_vis) {
        # figure out whether it was detected or not
        obsocc[i,((t-1)*n_vis)+j] <- rbinom(1,1,tocc[i,t]*p)
      }
    }
  }

  ### Output all the randomly generated pieces separately
  out_list <- list(true_occ = tocc, obs_occ = obsocc, phi_survival = year.effect.phi, gamma_colonization = gamma_year, cat_var = cat_var)
}

test_sim <- sim_dataset(perc_red = perc_red, psi1 = psi1, nyear = nyear, n_sites = n_sites,
                        n_vis = 2, phi = phi, sd_phi = sd_phi, gamma = gamma, sd_gamma = sd_gamma, p = p[1],
                        cat_var_list = psi_cat_var)

measure_decline <- purrr::map(1:10, ~sim_dataset(perc_red = perc_red, psi1 = psi1, nyear = nyear, n_sites = n_sites,
                        n_vis = 2, phi = phi, sd_phi = sd_phi, gamma = gamma, sd_gamma = sd_gamma, p = p[1],
                        cat_var_list = psi_cat_var)$true_occ[15603:31204, 10])
