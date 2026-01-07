fit_model_reps <- function(chunk, reps, n_year, n_visit, data, method = c("recursive", "equilibrium", "const_phi"), path, save_ending = ""){
  
  # get dataframe of scenarios
  map_df <- data %>%
    dplyr::filter(chunk_num == chunk) %>%
    select(-chunk_num)
  
  sample_size <- unique(map_df$total_n)
  
  # create save out directory
  dir.create(path)
  
  #### Set up nimble model for that sample size ####
  
  if (method %in% c("recursive", "equilibrium")){
    # Model code for single EMU year estimate
    dynoccmod_code <- nimble::nimbleCode({
      
      # The whole likelihood for the dynamic occupancy model is contained inside
      # dDynOcc_sss. The suffix _sss indicates that persistence, colonization, and
      # detection are provided as scalars (one value for the whole site's
      # detection history). Other variants exist with suffixes like _svm (which
      # would mean that persistence is (s)calar, colonization is a (v)ector
      # varying with season, and detection is a (m)atrix varying with season and
      # with replicate)
      
      for (i in 1:nsite) {
        y[i, 1:nseason, 1:nrep] ~ dDynOcc_vvs(probPersist = persist[1:(nseason-1)],
                                              probColonize = colonize[1:(nseason-1)],
                                              init = init_occ,
                                              p = detect,
                                              start = start_indexes[1:nseason], # Start and end arguments allow you to provide ragged mtx data
                                              end = end_indexes[1:nseason])
        
        
      }
      
      # Define Priors
      for (i in 1:(nseason-1)){
        persist_intercept[i] ~ dunif(0,1)
        logit(persist[i]) <- logit(persist_intercept[i]) + logit(beta) * i
        
        # random intercept for colonization
        colonize[i] ~ dunif(0,1)
      }
      
      # priors
      beta ~ dunif(0,1)
      #colonize ~ dunif(0,1)
      init_occ ~ dunif(0,1)
      detect ~ dunif(0,1)
      
      # Derive posterior for year
      psi[1] <-  init_occ
      for (i in 2:nseason){
        psi[i] <- psi[i-1]*(persist[i-1]) + (1-psi[i-1])*colonize[i-1]
      }
      perc_change <- (psi[10] - psi[1])/psi[1]
    })
  } else{
    # Model code for single EMU year estimate
    dynoccmod_code <- nimble::nimbleCode({
      
      # The whole likelihood for the dynamic occupancy model is contained inside
      # dDynOcc_sss. The suffix _sss indicates that persistence, colonization, and
      # detection are provided as scalars (one value for the whole site's
      # detection history). Other variants exist with suffixes like _svm (which
      # would mean that persistence is (s)calar, colonization is a (v)ector
      # varying with season, and detection is a (m)atrix varying with season and
      # with replicate)
      
      for (i in 1:nsite) {
        y[i, 1:nseason, 1:nrep] ~ dDynOcc_vvs(probPersist = persist[1:(nseason-1)],
                                              probColonize = colonize[1:(nseason-1)],
                                              init = init_occ,
                                              p = detect,
                                              start = start_indexes[1:nseason], # Start and end arguments allow you to provide ragged mtx data
                                              end = end_indexes[1:nseason])
        
        
      }
      
      # Define Priors
      for (i in 1:(nseason-1)){
        persist[i] ~ dunif(0,1)
        
        # random intercept for colonization
        colonize[i] ~ dunif(0,1)
      }
      
      #colonize ~ dunif(0,1)
      init_occ ~ dunif(0,1)
      detect ~ dunif(0,1)
      
      # Derive posterior for year
      psi[1] <-  init_occ
      for (i in 2:nseason){
        psi[i] <- psi[i-1]*(persist[i-1]) + (1-psi[i-1])*colonize[i-1]
      }
      perc_change <- (psi[10] - psi[1])/psi[1]
    })
  }
  
  # create template model that can be updated with data
  compile_model <- init_model(n = sample_size, year = n_year, visit = n_visit, model_obj = dynoccmod_code, method = method)
  
  ## Map across scenarios for EMU and reps (multiple scenarios with the same sample size for each EMU)
  purrr::pmap(map_df, function(high_name, psi, phi, sd_phi, sd_gamma, p, perc_red,  total_n,
                               high_n, low_n, high_hex_count, low_hex_count, low_name, low_psi, seed, 
                               scenario_id, low_gamma, high_gamma, year, visit){
    
    set.seed(seed)
    
    for (i in 1:reps){
      
      if(method == "recursive"){
        
        #simulate high occupancy
        high_data <- sim_dataset(psi = psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                                 n_sites = high_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
          append(c("sim_id" = high_name, "rep" = i))
        
        #simulate low occupancy
        low_data <- sim_dataset(psi = low_psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                                n_sites = low_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
          append(c("sim_id" = low_name, "rep" = i))
        
      } else if (method == "equilibrium"){
        
        #simulate high occupancy
        high_data <- sim_dataset_equil(psi = psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                                       n_sites = high_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
          append(c("sim_id" = high_name, "rep" = i))
        
        #simulate low occupancy
        low_data <- sim_dataset_equil(psi = low_psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p,
                                      n_sites = low_hex_count, perc_red = perc_red, nyear = year, n_vis = visit) %>%
          append(c("sim_id" = low_name, "rep" = i))
        
      } else {
        
        #simulate high occupancy
        high_data <- sim_dataset_constphi(psi = psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p, gamma = high_gamma,
                                          n_sites = high_hex_count, nyear = year, n_vis = visit) %>%
          append(c("sim_id" = high_name, "rep" = i))
        
        #simulate low occupancy
        low_data <- sim_dataset_constphi(psi = low_psi, phi = phi, sd_phi = sd_phi, sd_gamma = sd_gamma, p = p, gamma = low_gamma,
                                         n_sites = low_hex_count, nyear = year, n_vis = visit) %>%
          append(c("sim_id" = low_name, "rep" = i))
        
      }
      
      # save data out
      dir.create(paste0(path, "/emu_simulated_data", save_ending))
      
      saveRDS(high_data, paste0(path, "/emu_simulated_data", save_ending, "/", high_name, "_", i, "_", scenario_id, "_simdata.rds"))
      saveRDS(low_data, paste0(path, "/emu_simulated_data", save_ending, "/", low_name, "_", i, "_", scenario_id, "_simdata.rds"))
      
      #sample observed occupancy as model input
      sample_occ <- dplyr::bind_rows(sample_data(high_data$obs_occ, high_n),
                                     sample_data(low_data$obs_occ, low_n)) %>%
        dplyr::arrange(visit)
      
      occ_array <- sample_occ %>%
        dplyr::select(-site_id) %>%
        split(sample_occ$visit) %>%
        purrr::map( ~ .x |> dplyr::select(-visit) |> as.matrix()) %>%
        simplify2array()
      
      if (method %in% c("recursive", "equilibrium")){
        # list of new initialized variables
        new_inits <-  list(
          colonize = rep(0.5, (year-1)),
          init_occ = 0.5,
          detect = 0.5,
          persist_intercept = rep(0.5, (year-1)),
          beta = 0.5)
      } else{
        new_inits <-  list(
          colonize = rep(0.5, (year-1)),
          init_occ = 0.5,
          detect = 0.5,
          persist = rep(0.5, (year-1)))
      }
      
      # update model with data
      compile_model$mod$y <- occ_array
      #compile_model$setData("y")
      fit <- nimble::runMCMC(compile_model$mcmc, niter = 1000, nchains = 2, nburnin = 500,
                             samplesAsCodaMCMC = TRUE,
                             inits = new_inits)
      
      summary <- MCMCvis::MCMCsummary(fit, probs = c(0.025, 0.5, 0.95, 0.975)) |>
        dplyr::mutate(high_name = high_name, rep = i) |>
        tibble::rownames_to_column(var = "parameter")
      
      # create save out directories
      dir.create(paste0(path, "/emu_summary", save_ending))
      dir.create(paste0(path, "/emu_posterior", save_ending))
      
      readr::write_csv(summary, paste0(path, "/emu_summary", save_ending, "/", high_name, "_", i, "_", scenario_id, "_summary.csv"))
      saveRDS(fit, paste0(path, "/emu_posterior", save_ending, "/", high_name, "_", i, "_", scenario_id, "_posterior.rds"))
    }
  }, year = n_year, visit = n_visit)
}

sample_data <- function(data, sample_size){
  data_ids <- unique(data$site_id)
  samp_data <- data %>%
    filter(site_id %in% sample(data_ids, sample_size, replace = FALSE))
}

# function to initialize a model object for a given sample size
init_model <- function(n, year, visit, model_obj, method = c("recursive", "equilibrium", "const_phi")){
  
  nsite <- n
  nseason <- year
  nrep <- visit
  
  # make observed occurrence an array with site x year(season) x visit(rep)
  obs_occ_init <- array(sample(c(0,1), (nsite * nseason * nrep), replace = TRUE), c(nsite, nseason, nrep))
  
  if (method %in% c("recursive", "equilibrium")){
    inits_list <- list(
      colonize =rep(0.5, (nseason-1)),
      init_occ = 0.5,
      detect = 0.5,
      persist_intercept = rep(0.5, (nseason-1)),
      beta = 0.5
    )
  } else {
    inits_list <- list(
      colonize = rep(0.5, (nseason-1)),
      init_occ = 0.5,
      detect = 0.5,
      persist = rep(0.5, (nseason-1))
    )
  }
  
  # Build the model
  mod <- nimble::nimbleModel(
    code = model_obj,
    constants = list(nsite = n, nrep = nrep, nseason = nseason,
                     start_indexes = rep(1, nseason),
                     end_indexes = rep(nrep, nseason)
    ),
    data = list(y = obs_occ_init),
    inits = inits_list
  )
  
  # shouldn't NA, infinite, or positive (that's a dist issue)
  #Non-NA means we're fully initialized
  if (!is.finite(mod$calculate())){
    stop("Model did not initialize properly.")
  }
  
  # Build an MCMC
  conf <- nimble::configureMCMC(mod)
  conf$addMonitors(c("psi", "perc_change"))
  mcmc <- nimble::buildMCMC(conf)
  
  # Compile
  complist <- nimble::compileNimble(mod, mcmc)
  
  return(complist)
}
