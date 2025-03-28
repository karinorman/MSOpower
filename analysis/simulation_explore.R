####################################################################################################################
## Power Analysis Goal: Owl occupancy rates must show a stable or increasing trend after 10 years of monitoring.
## The study design to verify this criterion must have a power of 90% (Type II error rate β = 0.10) to detect a
## 25% decline in occupancy rate over the 10-year period with a Type I error rate (α) of 0.10.
####################################################################################################################

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
  ungroup() %>%
  #let's make up a single cell of low quality in BRW so everything doesn't break
  bind_rows(data.frame(UNIT = "Basin & Range - West", occupancy = "low", hex_count = 1, emu = "BRW"))


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
  select(sim_num, emu, psi, phi, sd_phi, sd_gamma, p, n = hex_count, perc_red)

# Let's get sample size of high quality hexes
# If an emu has enough area, we want the max sample size to be 2500, otherwise max sample is entire high quality area
sample_size_df <- emu_ratio %>%
  filter(occupancy == "high") %>%
  select(emu, hex_count) %>%
  mutate(log_max_samp = ifelse(hex_count > 2500, ceiling(log(2500)), ceiling(log(hex_count)))) %>%
  rowwise() %>%
  mutate(log_samp = list(seq(2.3, log_max_samp, by = 0.5))) %>%
  unnest(log_samp) %>%
  mutate(samp_size = round(exp(log_samp)))

emu_sample_sizes <- emu_ratio %>%
  filter(occupancy == "high") %>%
  select(emu, hex_count) %>%
  left_join(sample_size_df %>% select(emu, n_samp = samp_size)) %>%
  group_by(emu, n_samp) %>%
  # expand again to get a row for each simulation scenario
  slice(rep(row_number() , n_distinct(sim_scenarios$sim_num))) %>%
  mutate(sim_num = 1:n_distinct(sim_scenarios$sim_num))

sim_scenarios_emu <- sim_scenarios_emu %>%
  left_join(emu_sample_sizes %>% select(sim_num, emu, n_samp)) %>%
  group_by(emu) %>%
  mutate(sim_num = row_number()) %>%
  unite("sim_id", emu, sim_num)


###########################################
########### Generate data sets ############
###########################################

simn <- 100

#single_rep <- purrr::pmap(sim_scenarios_emu %>% select(-sim_id), sim_dataset, nyear = nyear, n_vis = 2) %>% set_names(sim_scenarios_emu$sim_id)

set.seed(42)
plan(multisession, workers = 15)
sim_list <- furrr::future_map(1:simn, ~purrr::pmap(sim_scenarios_emu %>%
                                                     select(-sim_id, -n_samp), sim_dataset, nyear = nyear, n_vis = 2) %>%
                                set_names(sim_scenarios_emu$sim_id),
                              .options=furrr_options(seed = TRUE)) %>%
  set_names(paste0("rep", 1:simn))

# reorder so top level of nested list is a sim scenario
sim_list_emu <- map(sim_scenarios_emu$sim_id, function(emu) {
  map(1:simn, ~pluck(sim_list, .x, emu)) %>%
    set_names(paste0("rep", 1:simn))}) %>%
  set_names(sim_scenarios_emu$sim_id)

#get true occurrence for each rep and sim
true_occ <- map_dfr(sim_scenarios_emu$sim_id, function(emu){
  map_dfr(1:simn, ~pluck(sim_list_emu, emu, .x, "true_occ") %>%
            select(-site_id) %>%
            ungroup() %>%
            summarize(across(everything(), mean)) %>%
            mutate(rep = .x)) %>%
    mutate(sim_id = emu)
})

true_occ_stats <- true_occ %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  select(-rep) %>%
  group_by(sim_id, time) %>%
  summarize(mean = mean(occ),
            lower = mean(occ) - qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n()),
            upper = mean(occ) + qt(1- 0.05/2, (n() - 1))*sd(occ)/sqrt(n())) %>%
  ungroup() %>%
  separate(sim_id, c("emu", "sim_num"), sep = "_", remove = FALSE) %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t"))) %>%
  left_join(sim_scenarios_emu)

# let's look at the annual reduction in survival for different scenarios
phi_red <- map_dfr(sim_scenarios_emu$sim_id, function(emu){
  map_dfr(1:simn, ~pluck(sim_list_emu, emu, .x, "phi_reduction") %>%
            data.frame(year = 2:9, phi_reduction = .) %>%
            mutate(rep = .x)) %>%
    mutate(sim_id = emu)
}) %>%
  group_by(sim_id, year) %>%
  summarize(phi_reduction = mean(phi_reduction)) %>%
  mutate(phi_multiplier = phi_reduction, phi_reduction = 1- phi_multiplier)


#This returns giant dataframe, hasn't been processed into encounter histories yet
# obs_occ <- map_dfr(sim_scenarios_emu$sim_id, function(emu){
#   map_dfr(1:3, ~pluck(sim_list_emu, emu, .x, "obs_occ") %>% mutate(rep = .x)) %>%
#     mutate(sim_id = emu)
# })

###########################################
########## Check Realized Trend ###########
###########################################
library(lme4)
library(broom.mixed)

true_occ_model_df <- true_occ %>%
  pivot_longer(starts_with("t"), names_to = "time", values_to = "occ") %>%
  mutate(time = as.numeric(stringr::str_remove(time, "t")))


model_fit_df <- true_occ_model_df %>%
  group_by(sim_id) %>%
  nest() %>%
  # fit model for each sim_id and cat variable
  mutate(model = map(data, ~lmer(occ ~ time + (1|rep), data = .x) %>% broom.mixed::tidy())) %>%
  select(-data) %>%
  unnest(model) %>%
  # get slope and intercept for mean effect
  filter(term %in% c("time", "(Intercept)")) %>%
  select(-std.error, -statistic, -group, -effect) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`) %>%
  mutate(t10 = intercept + (time * 10),
         t1 =  intercept + time,
         percent_change = ((t10 - t1)/abs(t1))) %>%
  separate(sim_id, c("emu", "sim_num"), sep = "_", remove = FALSE) %>%
  left_join(sim_scenarios_emu)



###########################################
########### Visualize True Occ ############
###########################################
library(ggplot2)

true_occ_plot_df <- true_occ_stats %>%
  group_by(phi, psi, p) %>%
  mutate(line_id = cur_group_id()) %>%
  left_join(emu_sample_sizes %>% select(emu, n_samp))

true_occ_plot_df %>%
  group_by(emu, n_samp) %>%
  slice(1) %>%
  #filter(percent_samp == 0.5) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(line_id)), alpha = 0.3) +
  geom_line(aes(color = as.factor(line_id))) +
  theme_classic() +
  facet_wrap(~emu, scales = "free") +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted") +
  geom_hline(yintercept = 0.2, linetype = "dotted") +
  geom_hline(yintercept = 0.15, linetype = "dotted") +
  geom_hline(yintercept = 0.43, linetype = "dotted") +
  geom_hline(yintercept = 0.3225, linetype = "dotted")

true_occ_plot_df %>%
  group_by(emu, n_samp) %>%
  slice(1) %>%
  filter(psi == 0.03) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(line_id)), alpha = 0.3) +
  geom_line(aes(color = as.factor(line_id))) +
  theme_classic() +
  facet_wrap(~emu, scales = "free") +
  scale_color_discrete(name = "Sim Scenario") +
  scale_fill_discrete(name = "Sim Scenario") +
  geom_hline(yintercept = 0.03, linetype = "dotted") +
  geom_hline(yintercept = 0.0225, linetype = "dotted")

# visualize annual reduction in survival to get the desired trend
phi_red %>%
  left_join(sim_scenarios_emu) %>%
  separate(sim_id, c("emu", "sim_num"), sep = "_", remove = FALSE) %>%
  filter(p == 0.8) %>%
  ggplot(aes(x = year, y = phi_reduction)) +
  geom_line(aes(color = sim_num)) +
  facet_wrap(~emu, scales = "free")

###########################################
########### Sampling Protocol ############
###########################################

# map high occupancy sims to their low occupancy counterpart
sim_map <- sim_scenarios_emu %>%
  select(sim_id, psi, phi, p, n_samp) %>%
  filter(psi != 0.03) %>%
  separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
  select(-psi) %>%
  rename(high_name = sim_id) %>%
  left_join(sim_scenarios_emu %>%
              select(sim_id, psi, phi, p, n_samp) %>%
              filter(psi == 0.03) %>%
              separate(sim_id, c("emu"), sep = "_", remove = FALSE) %>%
              select(-psi) %>%
              rename(low_name = sim_id)) %>%
  select(high_name, low_name, sample_size = n_samp)

model_check <- function(high_name, low_name, sample_size, repn) {

  print(c(high_name, repn))

  low_obs <- pluck(sim_list_emu, low_name, repn, "obs_occ") %>% mutate(landtype = "low")
  high_obs <- pluck(sim_list_emu, high_name, repn, "obs_occ") %>% mutate(landtype = "high")

  high_n = sample_size
  # if we have enough low quality samples, this is the ratio
  high_n = sample_size
  low_n = round(high_n*(1/3))

  if(low_n > n_distinct(low_obs$site_id)){
    low_n = n_distinct(low_obs$site_id)
  }

  obs_occ_df <- bind_rows(low_obs %>% filter(site_id %in% sample(unique(low_obs$site_id), low_n, replace = FALSE)),
                          high_obs %>% filter(site_id %in% sample(unique(high_obs$site_id), high_n, replace = FALSE))) %>%
    arrange(visit)

  obs_occ_array <- obs_occ_df %>%
    select(-site_id, landtype) %>%
    split(obs_occ_df$visit) %>%
    map(., ~ .x %>% select(-visit, -landtype) %>% as.matrix()) %>%
    simplify2array()

  if (dim(obs_occ_array)[1] != sample_size){
    stop("incorrect realized sample size")
  }

  landtype_cov <- obs_occ_df %>%
    filter(visit == 1) %>%
    select(landtype)

  year_cov <- matrix(1:nyear, nrow = 1)
  year_cov <- year_cov %x% rep(1, dim(obs_occ_array)[1])

  occ.covs <- list(landtype = landtype_cov, year = year_cov)

  #fit_model
  n.chains <- 3
  n.thin <- 1
  n.burn <- 500
  n.batch <- 30
  batch.length <- 25


  z.init <- apply(obs_occ_array, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
  inits.list <- list(beta = 0,
                     alpha = 0,
                     z = z.init)

  prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                     alpha.normal = list(mean = 0, var = 2.72))

  test_fit <- tPGOcc(occ.formula = ~ year + landtype,
                     det.formula = ~ 1,
                     data = list(y = obs_occ_array, occ.covs = occ.covs),
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

  # null_fit <- tPGOcc(occ.formula = ~ year,
  #                    det.formula = ~ 1,
  #                    data = list(y = obs_occ_array, occ.covs = occ.covs),
  #                    inits = inits.list,
  #                    priors = prior.list,
  #                    n.omp.threads = 1,
  #                    verbose = TRUE,
  #                    n.report = 750,
  #                    n.burn = n.burn,
  #                    n.thin = n.thin,
  #                    n.chains = n.chains,
  #                    n.batch = n.batch,
  #                    batch.length = batch.length)

  # get the posterior for the estimates
  post <- as.data.frame(test_fit$beta.samples) %>%
    rename(intercept = `(Intercept)`) %>%
    #mutate(across(everything(), plogis)) %>%
    mutate(t10 = plogis(year*10 + intercept),
           t1 = plogis(year + intercept),
           perc_change = (t10-t1)/t1) %>%
    mutate(sim_id = high_name, rep = repn)

  # post_null <- as.data.frame(null_fit$beta.samples) %>%
  #   rename(intercept = `(Intercept)`) %>%
  #   #mutate(across(everything(), plogis)) %>%
  #   mutate(t10 = plogis(year*10 + intercept),
  #          t1 = plogis(year + intercept),
  #          perc_change = (t10-t1)/t1) %>%
  #   mutate(sim_id = high_name, rep = repn)

  true_trend <- true_occ %>%
    filter(sim_id == high_name, rep == repn) %>%
    mutate(perc_change = (t10-t1)/t1) %>%
    pull(perc_change)

  check_dist <- between(true_trend, min(post$perc_change),max(post$perc_change)) & !between(0, min(post$perc_change),max(post$perc_change))
  #check_dist_null <- between(true_trend, min(post_null$perc_change),max(post_null$perc_change)) & !between(0, min(post_null$perc_change),max(post_null$perc_change))

  return(list("power_check" = data.frame("sim_id" = high_name, "rep" = repn, "low_n" = low_n, "high_n" = high_n,
             "true_perc_change" = true_trend, "est_perc_change" = mean(post$perc_change), #"est_perc_change_null" = mean(post_null$perc_change),
             "success" = check_dist, "samps_under" = sum(true_trend > post$perc_change),
             "samps_over" = sum(true_trend < post$perc_change)),
             "posterior" = post#,
             #"posterior_null" = post_null
             ))
}

pwr_check <- model_check(high_name = "BRE_21", low_name = "BRE_1", sample_size = 40, repn = 1)

set.seed(42)
plan(multisession, workers = 15)
power_check_list <- furrr::future_map(1:simn, function(x){ pmap(sim_map, model_check, repn = x)}, .options=furrr_options(seed = TRUE))

power_check_df <- map_dfr(1:100, function(y){map_dfr(1:200, ~pluck(power_check_list, y, .x, "power_check"))})

power_plot_df <- power_check_df %>%
  left_join(null_posterior %>% rename(null_success = check_dist)) %>%
  #select(sim_id, simn, low_n, high_n, success, null_success) %>%
  group_by(sim_id) %>%
  summarize(across(c(success, null_success), ~sum(.x)/simn)) %>%
  left_join(sim_scenarios_emu) %>%
  group_by(psi, p, phi) %>%
  mutate(line_id = cur_group_id()) %>%
  separate(sim_id, c("emu", "sim_num"), sep = "_", remove = FALSE)

power_plot_df %>%
  ggplot(aes(x = n_samp, y = success)) +
  geom_line(aes(color = as.factor(line_id))) +
  facet_wrap(~emu, scales = "free") +
  theme_classic() +
  scale_color_discrete(name = "Scenario")

power_plot_df %>%
  ggplot(aes(x = n_samp, y = null_success)) +
  geom_line(aes(color = as.factor(line_id))) +
  facet_wrap(~emu, scales = "free") +
  theme_classic() +
  scale_color_discrete(name = "Scenario")
