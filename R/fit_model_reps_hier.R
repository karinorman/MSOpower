fit_model_reps_hier_working <- function(chunk, reps, n_year, n_visit, data, path, save_ending = ""){

  # get dataframe of scenarios
  map_df <- data %>%
    dplyr::filter(chunk_num == chunk)

  sample_size <- unique(map_df$total_n)

  # create save out directory
  dir.create(path)

  #### Set up nimble model for that sample size ####

  dynoccmod_code <- nimble::nimbleCode({

    # The whole likelihood for the dynamic occupancy model is contained inside
    # dDynOcc_sss. The suffix _sss indicates that persistence, colonization, and
    # detection are provided as scalars (one value for the whole site's
    # detection history). Other variants exist with suffixes like _svm (which
    # would mean that persistence is (s)calar, colonization is a (v)ector
    # varying with season, and detection is a (m)atrix varying with season and
    # with replicate)

    for (i in 1:nsite) {
      y[i, 1:nseason, 1:nrep] ~ dDynOcc_vvs(probPersist = persist[i, 1:(nseason-1)],
                                            probColonize = colonize[1:(nseason-1)],
                                            init = init_occ,
                                            p = detect,
                                            start = start_indexes[1:nseason], # Start and end arguments allow you to provide ragged mtx data
                                            end = end_indexes[1:nseason])


    }

    # Define Priors
    for (i in 1:(nseason-1)){

      persist_int[i] ~ dunif(0,1)
      # random intercept for colonization
      colonize[i] ~ dunif(0,1)

      for (j in 1:nsite){
        logit(persist[j, i]) <- logit(persist_int[i]) + logit(beta) * i + ranef[EMU[j]]
        #logit(persist[j, i]) <- logit(persist_int[i]) + ranef[EMU[j]]
      }
    }

    for (r in 1:num_EMU) {
      # do sd = so life isn't ruined (might think its precision)
      ranef[r] ~ dnorm(0, sd = sigma_ranef)
    }

    beta ~ dunif(0,1)
    init_occ ~ dunif(0,1)
    detect ~ dunif(0,1)
    sigma_ranef ~ dunif(0, 10)

    # Derive posterior for year
    psi[1] <-  init_occ
    for (i in 2:nseason){
      # gives the estimate for year based on mean persistence (not a level of random effect)
      logit(derived_persist[i-1]) <- logit(persist_int[i-1]) + logit(beta) * (i-1)
      psi[i] <- psi[i-1]*(derived_persist[i-1]) + (1-psi[i-1])*colonize[i-1]
    }

    perc_change <- (psi[10] - psi[1])/psi[1]
  })

  # create template model that can be updated with data
  compile_model <- init_model(n = sample_size, year = n_year, visit = n_visit, emu_vec = 1:n_distinct(hex_count$emu), model_obj = dynoccmod_code)

  ## Map across scenarios for EMU and reps (multiple scenarios with the same sample size for each EMU)
  purrr::pmap(map_df, function(chunk_num, psi, phi, sd_phi, sd_gamma, p, perc_red,  total_n,
                               high_n, low_n, high_hex_count, low_hex_count, low_psi, seed,
                               year, visit){
    #browser()
    set.seed(seed)

    low_emu_vec <- unlist(pmap(hex_count, function(emu, high_hex_count, low_hex_count) {rep(emu, time = low_hex_count)}))
    high_emu_vec <- unlist(pmap(hex_count, function(emu, high_hex_count, low_hex_count) {rep(emu, time = high_hex_count)}))

    for (i in 1:reps){

      low_emu_vec <- pmap(hex_count, function(emu, high_hex_count, low_hex_count) {rep(emu, time = low_hex_count)}) %>%
        unlist() %>%
        sample()

      high_emu_vec <- pmap(hex_count, function(emu, high_hex_count, low_hex_count) {rep(emu, time = high_hex_count)}) %>%
        unlist() %>%
        sample()

      #simulate high occupancy
      high_data <- sim_dataset(psi = psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                               n_sites = high_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
        append(c(occ_type = "high", chunk_num = chunk_num, "rep" = i))

      # assign sites to EMU
      high_data$true_occ <- high_data$true_occ %>%
        bind_cols(data.frame(emu = high_emu_vec))

      high_data$obs_occ <- high_data$obs_occ %>%
        arrange(visit, site_id) %>%
        bind_cols(data.frame(emu = rep(high_emu_vec, 2))) %>%
        mutate(landtype = "high")

      #simulate low occupancy
      low_data <- sim_dataset(psi = low_psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                              n_sites = low_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
        append(c(occ_type = "low", chunk_num = chunk_num, "rep" = i))

      # assign sites to EMU
      low_data$true_occ <- low_data$true_occ %>%
        bind_cols(data.frame(emu = low_emu_vec))

      low_data$obs_occ <- low_data$obs_occ %>%
        arrange(visit, site_id) %>%
        bind_cols(data.frame(emu = rep(low_emu_vec, 2))) %>%
        mutate(landtype = "low")

      # save data out
      dir.create(paste0(path, "/emu_simulated_data", save_ending))

      saveRDS(high_data, paste0(path, "/emu_simulated_data", save_ending, "/", "high","_", chunk_num, "_", i, "_simdata.rds"))
      saveRDS(low_data, paste0(path, "/emu_simulated_data", save_ending, "/", "low", "_", chunk_num, "_", i,  "_simdata.rds"))

      # sample observed occupancy
      occ_data <- bind_rows(sample_data_hier(high_data$obs_occ, high_n),
                            sample_data_hier(low_data$obs_occ, low_n)) %>%
        # enforce that visits are ordered as that emu vec is same across visits
        dplyr::arrange(visit, emu)

      occ_array <- occ_data %>%
        dplyr::select(-c(site_id, emu, landtype)) %>%
        split(occ_data$visit) %>%
        purrr::map( ~ .x |> dplyr::select(-visit) |> as.matrix()) %>%
        simplify2array()

      emu_var <-  occ_data %>%
        filter(visit == 1) %>%
        mutate(emu_num = as.numeric(as.factor(emu))) %>%
        pull(emu_num)

      # list of new initialized variables
      new_inits <-  list(
        beta = 0.5,
        colonize = rep(0.5, (n_year-1)),
        init_occ = 0.5,
        detect = 0.5,
        persist_int = rep(0.5, (n_year-1)),
        ranef = rep(0, n_distinct(emu_var)),
        sigma_ranef = 1,
        EMU = emu_var)

      # update model with data
      compile_model$mod$y <- occ_array
      #compile_model$setData("y")
      fit <- nimble::runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                             samplesAsCodaMCMC = TRUE,
                             inits = new_inits)

      summary <- MCMCvis::MCMCsummary(fit, probs = c(0.025, 0.5, 0.95, 0.975)) |>
        dplyr::mutate(chunk_num = chunk_num, rep = i) |>
        tibble::rownames_to_column(var = "parameter")

      # create save out directories
      dir.create(paste0(path, "/emu_summary", save_ending))
      dir.create(paste0(path, "/emu_posterior", save_ending))

      readr::write_csv(summary, paste0(path, "/emu_summary", save_ending, "/", chunk_num, "_", i, "_summary.csv"))
      saveRDS(fit, paste0(path, "/emu_posterior", save_ending, "/", chunk_num, "_", i,"_posterior.rds"))
    }
  }, year = n_year, visit = n_visit)
}

sample_data_hier <- function(data, sample_size){

  # get even sample across emus
  n_emu <- n_distinct(data$emu)
  emu_samp_size <- ceiling(sample_size/n_emu)

  # get site_id's for those sampled
  samp_sites <- data %>%
    filter(visit == 1) %>%
    group_by(emu) %>%
    slice_sample(n = emu_samp_size, replace = FALSE) %>%
    pull(site_id)

  # check that the rounded sample size gets us the right number
  realized_samp_n <- length(samp_sites)

  if (realized_samp_n != sample_size){
    remove_n <- realized_samp_n - sample_size

    remove_index <- sample(1:realized_samp_n, remove_n)
    samp_sites <- samp_sites[-remove_index]
  }

  return(data %>% filter(site_id %in% samp_sites))
}

# function to initialize a model object for a given sample size
# function to initialize a model object for a given sample size
init_model <- function(n, year, visit, emu_vec, model_obj){

  nsite <- n
  nseason <- year
  nrep <- visit

  # make observed occurrence an array with site x year(season) x visit(rep)
  obs_occ_init <- array(sample(c(0,1), (nsite * nseason * nrep), replace = TRUE), c(nsite, nseason, nrep))

  # Build the model
  mod <- nimble::nimbleModel(
    code = model_obj,
    constants = list(nsite = n, nrep = nrep, nseason = nseason,
                     start_indexes = rep(1, nseason),
                     end_indexes = rep(nrep, nseason),
                     num_EMU = n_distinct(emu_vec)
    ),
    data = list(y = obs_occ_init),
    inits = list(
      beta = 0.5,
      colonize = rep(0.5, (nseason-1)),
      init_occ = 0.5,
      detect = 0.5,
      persist_int = rep(0.5, (nseason-1)),
      ranef = rep(0, n_distinct(emu_vec)),
      sigma_ranef = 1,
      EMU = sample(emu_vec, nsite, replace = TRUE)
    )
  )

  # shouldn't NA, infinite, or positive (that's a dist issue)
  #Non-NA means we're fully initialized
  if (!is.finite(mod$calculate())){
    stop("Model did not initialize properly.")
  }

  # Build an MCMC
  conf <- nimble::configureMCMC(mod)
  conf$addMonitors(c("psi", "perc_change", "ranef"))
  mcmc <- nimble::buildMCMC(conf)

  # Compile
  complist <- nimble::compileNimble(mod, mcmc)

  return(complist)
}
