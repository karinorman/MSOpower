sim_dataset <- function(psi, phi, sd_phi, sd_gamma, p, n_sites, perc_red, nyear, n_vis){

  ### psi = initial occupancy
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
  #cat_n <- setNames(c(low_n, high_n), c("low", "high"))
  #psi1 <- setNames(c(psi1_low, psi1_high), c("low", "high"))
  # derive from psi1 = gamma/(gamma + epsilon); occupancy at equilibrium
  gamma <- (psi * (1-phi))/(1-psi)
  #nsites <- sum(cat_n)

  ### simulate seasonal effects on local colonization(using an
  ### additive term on the logit scale drawn from N(0,sd_gamma)

  ### year effects on colonization without trend
  logit.gamma.mean <- qlogis(gamma)

  # subsequent yearly gammas with noise
  gamma_year <- plogis(logit.gamma.mean + rnorm((nyear-1), 0, sd_gamma))

  #if we want occupancy to decline by perc_red over t-1, this is the annual decrease
  annual_perc_red = exp(log(1-perc_red)/(nyear-1))

  # and this is the annual decrease in extinction to get to that annual survival decrease
  #annual_ext_red = (gamma - (annual_perc_red*psi*gamma))/(annual_perc_red*psi*(1-phi))


  # what should the true occupancy be if we have that annual reduction, start with initial occupancy
  true_psi = c(psi, map(1:9, ~ psi * (annual_perc_red ^ .x))) %>% unlist()

  # initialize vector for survival reduction
  annual_phi_c = rep(0,nyear-1)
  ### year effects on survival with simulate yearly decreases
  phi_year = rep(0,nyear-1)
  phi_year[1] = phi

  for (j in 2:(nyear-1)) {
    annual_phi_c[j] = ((annual_perc_red*true_psi[j-1]) - (gamma*(1-true_psi[j-1])))/(true_psi[j-1]*phi_year[j-1])

    phi_year[j] = phi_year[j-1]*annual_phi_c[j]
  }



  # # create phi for each time step with a reduction from the previous year's phi
  # for (j in 2:(nyear-1)) {
  #   phi_year[j] = phi_year[j-1]*annual_phi_c[j]
  # }

  # get noise around phi at each time step drawn from N(phi, sd_phi)
  phi_year<- purrr::map(phi_year, ~rnorm(1, .x, sd_phi)) %>% unlist()

  ### Simulate tocc and obsocc from t = 1 to t = nyear
  ### First use initial values to generate year 1
  t1 <- rbinom(n_sites, 1,  psi)

  # add site identifier for categorical variable
  tocc <- data.frame(t1) %>%
    mutate(site_id = row_number()) %>%
    select(site_id, t1)

  # function that gets occurrence for next time step
  occ_tplus1 <- function(occ_t, t){
    if (occ_t == 1){
      occ_tplus1 <- rbinom(1,1, phi_year[t-1])
    } else {
      occ_tplus1 <- rbinom(1,1, gamma_year[t-1])
    }
    return(occ_tplus1)
  }

  #create new time step based on previous timestep
  for(t in 2:nyear){
    prev_year <- paste0("t", t-1)
    new_year <- paste0("t", t)

    tocc <- tocc %>%
      rowwise() %>%
      mutate({{new_year}} := occ_tplus1(occ_t = get(!!prev_year), t = t)) %>%
      ungroup()
  }

  obsocc <- purrr::map_dfr(1:n_vis, ~tocc %>%
    rowwise() %>%
    mutate(across(-c(site_id), ~rbinom(1, 1, .x*p))) %>%
      mutate(visit = .x)) %>%
    ungroup()

  ### Output all the randomly generated pieces separately
  out_list <- list(true_occ = tocc, obs_occ = obsocc,
                   phi_survival = phi_year,
                   gamma_colonization = gamma_year)

}

