############################################
########### Hierarchical Power #############
############################################

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(spOccupancy)

# real world vegtypes for each emu
emu_veg <- read.csv(here::here("data/EMU_veg_types.csv")) %>%
  # let's say which we think has high or low occupancy
  mutate(occupancy = case_when(
    veg_type_landfire == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "high",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "high",
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
    veg_type_landfire == "Madrean Pinyon-Juniper Woodland" ~ "low",
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

# get dataframe of all possible scenarios
sim_scenarios <- data.frame(
  # these are the parameters that change, taken directly from Woods 2019
  psi = c(0.03, 0.43, 0.6), phi = c(.6, .8, .8), p = c(0.4, 0.8, .8)) %>%
  # get all possible combinations
  tidyr::expand(psi, phi,p) %>%
  # and give each unique combination an ID
  mutate(sim_num = row_number(),
         occupancy = ifelse(psi == 0.03, "low", "high")) %>%
  # these are the same for all scenarios right now, sd's from Wood 2019, psi1_low kinda made up
  mutate(sd_phi = 0.04, sd_gamma = 0.01, perc_red = 0.25)

high_hex_count <- emu_ratio %>% filter(occupancy == "high") %>% pull(hex_count) %>% sum()
#sample_sizes <- seq(log(100), log(3000), by = 0.3) %>% exp() %>% round()
sample_sizes <- seq(log(100), log(high_hex_count/2), by = 0.5) %>% exp() %>% round()

# get scenarios, one for each emu
sim_scenarios_emu <- bind_rows(sim_scenarios %>% mutate(emu = "BRE"),
                               sim_scenarios %>% mutate(emu = "BRW"),
                               sim_scenarios %>% mutate(emu = "CP"),
                               sim_scenarios %>% mutate(emu = "SRM"),
                               sim_scenarios %>% mutate(emu = "UGM")) %>%
  # get sample sizes for low and high occupancy for each emu
  left_join(emu_ratio %>%
              select(emu, occupancy, hex_count)) %>%
  # tidyr::pivot_wider(names_from = occupancy, values_from = hex_count) %>%
  # rename(low_n = low, high_n = high)) %>%
  # get the columns in the right order
  select(sim_num, emu, psi, phi, sd_phi, sd_gamma, p, n = hex_count, perc_red) %>%
  filter(!is.na(n)) %>%
  group_by(emu, sim_num) %>%
  slice(rep(row_number(), length(sample_sizes))) %>% mutate(total_samp = sample_sizes) %>%
  group_by(emu) %>%
  mutate(sim_num = row_number()) %>%
  unite("sim_id", emu, sim_num)

# let's figure out relative area of each emu and total area
emu_ratio %>%
  filter(occupancy == "high") %>%
  select(emu, hex_count) %>%
  mutate(total_hex = sum(hex_count)) %>%
  mutate(proportion = hex_count/total_hex)


###########################################
########### Generate data sets ############
###########################################

# Let's not worry about how the samples are distributed, just get a toy set of simulations to play with

source(here::here("R/sim_dataset.R"))
simn <- 100

plan(multisession, workers = 60)
hier_sim_list <- map(1:simn, ~furrr::future_pmap(sim_scenarios_emu %>%
                                                   select(-sim_id, -total_samp), sim_dataset, nyear = nyear, n_vis = 2) %>%
                       set_names(sim_scenarios_emu$sim_id),
                     .options=furrr_options(seed = TRUE)) %>%
  set_names(paste0("rep", 1:simn))

usethis::use_data(hier_sim_list)

# reorder so top level of nested list is a sim scenario
sim_list_emu <- map(sim_scenarios_emu$sim_id, function(emu) {
  map(1:simn, ~pluck(hier_sim_list, .x, emu)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios_emu$sim_id)

### Let's get the dataset aggregated across emu's
# map high occupancy sims to their low occupancy counterpart
sim_map <- sim_scenarios_emu %>%
  select(sim_id, psi, phi, p, total_samp) %>%
  filter(psi != 0.03) %>%
  separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
  rename(high_name = sim_id) %>%
  left_join(sim_scenarios_emu %>%
              select(sim_id, psi, phi, p, total_samp) %>%
              filter(psi == 0.03) %>%
              separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
              select(-psi) %>%
              rename(low_name = sim_id), by = c("emu", "total_samp", "phi", "p")) %>%
  group_by(phi, p, psi, total_samp) %>%
  mutate(scenario_id = cur_group_id())


# get list of sim id's for each scenario
sim_reps_df <- map(unique(sim_map$scenario_id), function(id) {
  df <- sim_map %>% filter(scenario_id == id)
  
  return(c(df$high_name, df$low_name))
}) %>%
  set_names(unique(sim_map$scenario_id)) %>%
  tibble::enframe() %>%
  unnest(value) %>%
  rename(emu_sim_name = value, scenario_id = name)

# need the true occ for each sim scenario and rep
true_occ <- map_dfr(unique(sim_reps_df$scenario_id), function(scenario){
  emu_sims <- sim_reps_df %>%
    filter(scenario_id == scenario) %>%
    pull(emu_sim_name)
  
  #browser()
  
  map_dfr(1:simn, function(simn, scenario_id) {
    #browser()
    map_dfr(emu_sims, ~pluck(sim_list_emu, .x, simn, "true_occ") %>%
              mutate(emu_sim_name = .x, rep = simn )) %>%
      mutate(scenario_id = scenario_id)
  }, scenario_id = scenario)
}) %>%
  #filter(emu_sim_name %in% sim_map$high_name) %>%
  select(-emu_sim_name, -site_id) %>%
  group_by(rep, scenario_id) %>%
  summarize(across(everything(), mean) )%>%
  mutate(perc_change = (t10-t1)/t1)



###########################################
############### Fit Models ################
###########################################

# get observed detection histories pooled for all EMU's
obs_occ <- map_dfr(unique(sim_reps_df$scenario_id), function(scenario){
  emu_sims <- sim_reps_df %>%
    filter(scenario_id == scenario) %>%
    pull(emu_sim_name)
  
  #browser()
  
  map_dfr(1:simn, function(simn, scenario_id) {
    #$browser()
    map_dfr(emu_sims, ~pluck(sim_list_emu, .x, simn, "obs_occ") %>%
              mutate(emu_sim_name = .x, rep = simn )) %>%
      mutate(scenario_id = scenario_id)
  }, scenario_id = scenario)
}) %>%
  left_join(sim_scenarios_emu %>%
              mutate(landtype = ifelse(psi == 0.03, "low", "high")) %>%
              select(emu_sim_name = sim_id, total_samp, landtype)) %>%
  mutate()

# nested dataframe where detection histories are ID'd by scenario_id, sample size, and the rep number
map_model_df <- obs_occ %>%
  group_by(scenario_id, rep, total_samp) %>%
  nest() %>%
  select(total_samp, scenario_id, repn = rep, model_data = data) %>%
  left_join(true_occ %>% 
              select(scenario_id, repn = rep, true_trend = perc_change)) %>%
  ungroup()

# gotta deal with the giant memory issues
rm(hier_sim_list, sim_list_emu, obs_occ)
gc()


hierarch_model_check <- function(total_samp, scenario_id, repn, model_data, true_trend){
  # get example sample where half the sites are sampled
  #total_samp <- 19552
  
  model_data <- model_data %>%
    group_by(site_id, emu_sim_name) %>%
    mutate(site_id = cur_group_id()) %>%
    ungroup()
  
  samp_sites <- c(sample(model_data %>% filter(landtype == "high") %>% pull(site_id) %>% unique(), round(total_samp*0.75), replace = FALSE),
                  sample(model_data %>% filter(landtype == "low") %>% pull(site_id) %>% unique(), round(total_samp*0.25), replace = FALSE)
  )
  
  model_data_samp <- model_data %>% filter(site_id %in% samp_sites)
  
  obs_occ_array <- model_data_samp %>%
    arrange(site_id) %>%
    select(-c(site_id, emu_sim_name, landtype)) %>%
    split(model_data_samp$visit) %>%
    map(., ~ .x %>% select(-visit) %>% as.matrix()) %>%
    simplify2array()
  
  covars <-  model_data_samp %>%
    arrange(site_id) %>%
    filter(visit == 1) %>%
    separate_wider_delim(emu_sim_name, delim = "_", names = c("emu", "old_sim")) %>%
    select(-old_sim) %>%
    select(emu, landtype) %>%
    group_by(emu) %>%
    mutate(emu_num = cur_group_id())
  
  year_cov <- matrix(1:nyear, nrow = 1)
  year_cov <- year_cov %x% rep(1, dim(obs_occ_array)[1])
  
  covar_list <- list(emu = covars$emu_num, landtype = covars$landtype, year = year_cov)
  
  #fit_model
  n.chains <- 3
  n.thin <- 1
  n.burn <- 2000
  n.batch <- 60
  batch.length <- 50
  
  z.init <- apply(obs_occ_array, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
  inits.list <- list(beta = 0,
                     alpha = 0,
                     z = z.init)
  
  prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                     alpha.normal = list(mean = 0, var = 2.72))
  
  
  test_fit <- tPGOcc(occ.formula = ~ year + (1 | emu),
                     det.formula = ~ 1,
                     data = list(y = obs_occ_array, occ.covs = covar_list),
                     inits = inits.list,
                     priors = prior.list,
                     n.omp.threads = 1,
                     verbose = TRUE,
                     n.report = 750,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains,
                     n.batch = n.batch,
                     batch.length = batch.length)
  
  post <- as.data.frame(test_fit$beta.samples) %>%
    rename(intercept = `(Intercept)`) %>%
    #mutate(across(everything(), plogis)) %>%
    mutate(t10 = plogis(year*10 + intercept),
           t1 = plogis(year + intercept),
           perc_change = (t10-t1)/t1) %>%
    mutate(scenario_id = scenario_id, rep = repn)
  # 
  #   true_trend <- true_occ %>%
  #     filter(scenario_id == !!scenario_id, rep == !!repn) %>%
  #     pull(perc_change)
  
  #check_dist <- between(true_trend, min(post$perc_change),max(post$perc_change)) & !between(0, min(post$perc_change),max(post$perc_change))
  
  readr::write_csv(data.frame("scenario_id" = scenario_id, "rep" = repn,
                              "true_perc_change" = true_trend, "est_perc_change" = mean(post$perc_change),
                              #"success" = check_dist, 
                              "samps_under" = sum(true_trend > post$perc_change),
                              "samps_over" = sum(true_trend < post$perc_change)),
                   here::here("data/hierarchical_output/power_check", paste(scenario_id, repn, "total_samp", "power_check.csv", sep = "_")))
  
  readr::write_csv(post, here::here("data/hierarchical_output/posterior", paste(scenario_id, repn, "total_samp", "posterior.csv", sep = "_")))
  
  # return(list("power_check" = data.frame("scenario_id" = scenario_id, "rep" = repn,
  #                                        "true_perc_change" = true_trend, "est_perc_change" = mean(post$perc_change),
  #                                        #"success" = check_dist, 
  #                                        "samps_under" = sum(true_trend > post$perc_change),
  #                                        "samps_over" = sum(true_trend < post$perc_change)),
  #             "posterior" = post
  #))
}

# perform the power check
dir.create(here::here("data/hierarchical_output"))
dir.create(here::here("data/hierarchical_output/posterior"))
dir.create(here::here("data/hierarchical_output/power_check"))

plan(multisession, workers = 60)
hierarch_check_list <- furrr::future_pmap(map_model_df, hierarch_model_check, .options=furrr_options(seed = TRUE))

# if it doesn't execute them all on first pass
exec_files <- data.frame(file_name = list.files(here::here("data/hierarchical_output/posterior"))) %>%
  separate(file_name, c("scenario_id", "repn"), sep = "_") %>%
  mutate(repn = as.integer(repn))

missing_scenario <- map_model_df %>% select(scenario_id, repn) %>%
  left_join(exec_files %>% mutate(check = "executed")) %>% 
  filter(is.na(check)) %>%
  select(-check) %>% 
  left_join(map_model_df) %>%
  filter(total_samp < 5000)

plan(multisession, workers = 30)
hierarch_check_list <- furrr::future_pmap(missing_scenario, hierarch_model_check, .options=furrr_options(seed = TRUE))

# get the pieces as two seperate dataframes
hier_power_check_df <- map_dfr(list.files(here::here("data/hierarchical_output/power_check/"), full.names = TRUE), read.csv)
hier_power_check_post_df <- map_dfr(list.files(here::here("data/hierarchical_output/posterior/"), full.names = TRUE), read.csv)

# hier_power_check_df <- map_dfr(1:length(hierarch_check_list), ~pluck(hierarch_check_list, .x, "power_check"))
# hier_power_check_post_df <- map_dfr(1:length(hierarch_check_list), ~pluck(hierarch_check_list, .x, "posterior"))

# save out so we don't have to re run
readr::write_csv(hier_power_check_df, here::here("data/hier_power_check_df.csv"))
readr::write_csv(hier_power_check_post_df, here::here("data/hier_power_check_post_df.csv"))

# build checks directly from sim outputs
# let's get the sim's with the sample sizes we care about
sim_keep <- sim_map %>% 
  filter(total_samp < 5000) %>%
  pull(scenario_id)

ci_alpha_list <- c(0.5, 0.90, 0.95)

ci_df <- hier_power_check_post_df %>%
  filter(scenario_id %in% sim_keep) %>%
  select(scenario_id, rep, perc_change) %>%
  group_by(scenario_id, rep) %>%
  reframe(CI_low = map(ci_alpha_list, 
                       ~as.data.frame(bayestestR::ci(perc_change, ci = .x, method = "ETI"))$CI_low) %>% unlist(),
          CI_high = map(ci_alpha_list, 
                        ~as.data.frame(bayestestR::ci(perc_change, ci = .x, method = "ETI"))$CI_high)  %>% unlist(),
          CI_type = ci_alpha_list,
          min_post = min(perc_change),
          max_post = max(perc_change)) %>%
  ungroup() %>%
  mutate(width = CI_high - CI_low)

power_plot_df <- hier_power_check_df %>%
  select(scenario_id, rep, true_perc_change, est_perc_change) %>%
  mutate(bias_comp = true_perc_change < est_perc_change,
         bias = est_perc_change - true_perc_change,
         relative_bias = (est_perc_change - true_perc_change)/true_perc_change) %>%
  left_join(ci_df) %>%
  rowwise() %>%
  mutate(ci_check = between(true_perc_change, CI_low, CI_high) & !between(0,  CI_low, CI_high),
         ci_check_left_tail = true_perc_change > CI_low & true_perc_change < max_post & max_post < 0,
         ci_check_any_decline = CI_high < 0,
         post_check_any_decline = max_post < 0,
         ci_check_interval_twotail = between(true_perc_change, CI_low, CI_high),
         ci_check_interval_righttail = between(true_perc_change, min_post, CI_high),
         # if the true trend isn't in the CI, what direction was the bias?
         bias_out_ci = ifelse(ci_check_interval_twotail == FALSE, bias_comp, NA)) %>%
  select(scenario_id, rep, bias, bias_comp, relative_bias, width, ci_check, ci_check_interval_twotail, ci_check_interval_righttail,
         ci_check_left_tail, ci_check_any_decline, post_check_any_decline,
         bias_out_ci, CI_type) %>%
  group_by(scenario_id, CI_type) %>%
  summarize(bias = mean(bias),
            relative_bias = mean(relative_bias),
            width = mean(width),
            ci_check = sum(ci_check)/simn,
            ci_check_left_tail = sum(ci_check_left_tail)/simn,
            ci_check_any_decline = sum(ci_check_any_decline)/simn,
            post_check_any_decline = sum(post_check_any_decline)/simn,
            ci_check_interval_twotail = sum(ci_check_interval_twotail)/simn,
            ci_check_interval_righttail = sum(ci_check_interval_righttail)/simn,
            percent_trend_lower = sum(bias) /simn,
            percent_exclude_trend_lower = sum(bias_out_ci, na.rm = TRUE) / sum(!is.na(bias_out_ci))) %>%
  left_join(sim_map %>% select(-c(high_name, low_name, emu)) %>% distinct()) #%>%
  # group_by(psi, p, phi) %>%
  # mutate(line_id = cur_group_id()) 

readr::write_csv(power_plot_df, here::here("data/hier_plot_df.csv"))
