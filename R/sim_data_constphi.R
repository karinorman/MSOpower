sim_dataset_constphi <- function(psi, phi, sd_phi, sd_gamma, p, gamma, n_sites, nyear, n_vis){

  ### psi = initial occupancy
  ### phi = annual survival
  ### sd_phi = sd of random normal variable introducing a Season effect
  ### on local survival
  ### sd_gamma = sd of random normal variable introducing a Season effect
  ### on colonization
  ### p = probability of detection
  ### nyear = number of years or study seasons
  ### n_vis = number of visits per site per year (or season)

  ### year effects on colonization without trend
  logit.gamma.mean <- qlogis(gamma)

  # subsequent yearly gammas with noise
  gamma_year <- plogis(logit.gamma.mean + rnorm((nyear-1), 0, sd_gamma))
  # enforce no probabilities greater than 1
  gamma_year[gamma_year > 1] <- 1

  # get noise around phi at each time step drawn from N(phi, sd_phi)
  phi_year <- phi + rnorm((nyear-1), 0, sd_phi)
  # enforce no probabilities greater than 1
  phi_year[phi_year > 1] <- 1

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
