library(dplyr)
library(tidyr)
library(lubridate)
library(spOccupancy)

##########################################
### Processing calls validated by Dana ###
##########################################

val_occ <- read.csv(here::here("data/top_conf_all_REID.csv")) %>%
  mutate(date = lubridate::mdy(date)) %>%
  rename(spow = "SPOW.") %>%
  filter(spow == "Y") %>%
  mutate(unit = stringr::str_replace(unit, "_[^_]+$", ""))

val_occ_day <- val_occ %>%
  select(forest, unit, date) %>%
  distinct()# %>%
  #mutate(occ = 1)

# the gila has detection at 7 units, the kaibab at 6
val_occ_day %>% count(forest, unit)

#####################
####### Kaibab ######
#####################

kaibab_meta <- read.csv(here::here("data/bioacoustics.2022.Kaibab.metadata.csv")) %>%
  pivot_longer(cols = c(deploy.swiftID, visit1.swiftID), names_to = "visit", values_to = "unit") %>%
  group_by(unit) %>%
  slice_head() %>%
  ungroup() %>%
  filter(!is.na(unit))

# figure out detection history periods, get start and end of sampling period for each unit
kaibab_meta_dates <- kaibab_meta %>%
  mutate(deploy_date = ymd(paste(cyear, deploy.month, deploy.day, sep = "-")),
         collect_date = ifelse(!is.na(visit2.month), paste(cyear, visit2.month, visit2.day, sep = "-"),
                                                         paste(cyear, visit1.month, visit1.day, sep = "-")),
         collect_date = ymd(collect_date),
         days_deployed = collect_date - deploy_date) %>%
  select(-c(cyear, deploy.month, deploy.day,  visit2.month, visit2.day))


date_ranges <- kaibab_meta_dates %>%
  select(pointID, unit, deploy_date, collect_date) %>%
  distinct() %>%
  group_by(unit) %>%
  nest() %>%
  mutate(week_range_start = purrr::map(data, ~seq(.$deploy_date, .$collect_date, by = "2 week")),
         week_range_end = purrr::map(data, ~c(tail(seq(.$deploy_date, .$collect_date, by = "2 week"), -1), .$collect_date))) %>%
  unnest(cols = c("data", "week_range_start", "week_range_end")) %>%
  filter(week_range_start != week_range_end) %>%
  select(-deploy_date, -collect_date) %>%
  left_join(val_occ_day %>%
              filter(forest == "kaibab") %>%
              select(-forest) %>%
              rename(obs_date = date)) %>%
  mutate(occ = ifelse(obs_date >= week_range_start & obs_date <= week_range_end, 1, 0),
         occ = ifelse(is.na(obs_date), 0, occ)) %>%
  group_by(pointID, week_range_start, week_range_end) %>%
  summarize(occ = sum(occ, na.rm = TRUE)) %>%
  ungroup()

kaibab_detect_hist <- date_ranges %>%
  select(-week_range_start, - week_range_end) %>%
  group_by(pointID) %>%
  mutate(time = paste0("t", row_number())) %>%
  pivot_wider(names_from = "time", values_from = "occ") %>%
  ungroup()

#### Let's try to fit a model ####

n.chains <- 3
n.samples <- 3000
n.thin <- 1
n.burn <- 2000

kaibab_occ_array <- kaibab_detect_hist %>%
  select(-pointID) %>%
  as.matrix()

z.init <- apply(kaibab_occ_array, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))

inits.list <- list(beta = 0,
                   alpha = 0)#,
                   #z = z.init)

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))


kaibab_fit <- PGOcc(occ.formula = ~ 1,
                   det.formula = ~ 1,
                   data = list(y = kaibab_occ_array),
                   inits = inits.list,
                   priors = prior.list,
                   n.omp.threads = 1,
                   verbose = TRUE,
                   n.report = 750,
                   n.burn = n.burn,
                   n.thin = n.thin,
                   n.chains = n.chains,
                   n.samples = n.samples)

#####################
####### Gila ########
#####################

gila_pickup <- readxl::read_excel("/Users/karinorman/Library/CloudStorage/Box-Box/Project R3 Bioacoustics Monitoring/Swift Units Database/Originals/Gila data.xlsx") %>%
  select(uniqueID = UniqueID, retrieved = Retrieved)

gila_meta <- read.csv(here::here("data/bioacoustics.2022.Gila.metadata.csv")) %>%
  left_join(gila_pickup)

# check burned units
burned_units <- gila_meta %>%
  filter(is.na(retrieved)) %>%
  pull(swiftID)

# no occurrences observed from burned units, so just exclude
val_occ %>%
  filter(unit %in% burned_units, forest == "gila") %>%
  View()

# need to fix so time periods are the same for all points in a hex
gila_detect <- gila_meta %>%
  # exclude burned units
  filter(!is.na(retrieved)) %>%
  mutate(deploy_date = ymd(paste(cyear, deploy.month, deploy.day, sep = "-")),
         collect_date = ymd(retrieved)) %>%
  group_by(hexID) %>%
  mutate(hex_deploy_date = min(deploy_date),
         hex_collect_date = max(collect_date)) %>%
  select(uniqueID, swiftID, hexID, hex_deploy_date, hex_collect_date) %>%
  group_by(uniqueID, swiftID, hexID) %>%
  nest() %>%
  mutate(week_range_start = purrr::map(data, ~seq(.$hex_deploy_date, .$hex_collect_date, by = "1 week")),
         week_range_end = purrr::map(data, ~c(tail(seq(.$hex_deploy_date, .$hex_collect_date, by = "1 week"), -1), .$hex_collect_date))) %>%
  unnest(cols = c("data", "week_range_start", "week_range_end")) %>%
  filter(week_range_start != week_range_end) %>%
  select(-hex_deploy_date, -hex_collect_date) %>%
  mutate(days_in_range = week_range_end - week_range_start) %>%
  filter(days_in_range > 2) %>%
  left_join(val_occ_day %>%
              filter(forest == "gila") %>%
              select(-forest) %>%
              rename(obs_date = date), by = c("swiftID" = "unit")) %>%
  mutate(occ = ifelse(obs_date >= week_range_start & obs_date <= week_range_end, 1, 0),
         occ = ifelse(is.na(obs_date), 0, occ)) %>%
  group_by(hexID, week_range_start, week_range_end) %>%
  summarize(occ = ifelse(sum(occ, na.rm = TRUE) > 0, 1, 0)) %>%
  ungroup()

gila_detect_hist <- gila_detect %>%
  select(-week_range_start, - week_range_end) %>%
  group_by(hexID) %>%
  mutate(time = paste0("t", row_number())) %>%
  pivot_wider(names_from = "time", values_from = "occ") %>%
  ungroup()

### Fit model ###
gila_occ_array <- gila_detect_hist %>%
  select(-hexID) %>%
  as.matrix()

z.init <- apply(gila_occ_array, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))

inits.list <- list(beta = 0,
                   alpha = 0)#,
#z = z.init)

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))


gila_fit <- PGOcc(occ.formula = ~ 1,
                    det.formula = ~ 1,
                    data = list(y = gila_occ_array),
                    inits = inits.list,
                    priors = prior.list,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 750,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains,
                    n.samples = n.samples)
