library(dplyr)
library(tidyr)
library(lubridate)

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
  distinct() %>%
  mutate(occ = 1)

# the gila has detection at 7 units, the kaibab at 6
val_occ_day %>% count(forest, unit)

####
# metadata for deployments
kaibab_meta <- read.csv(here::here("data/bioacoustics.2022.Kaibab.metadata.csv")) %>%
  pivot_longer(cols = c(deploy.swiftID, visit1.swiftID), names_to = "visit", values_to = "unit") %>%
  group_by(unit) %>%
  slice_head()

# figure out detection history periods, get start and end of sampling period for each unit
kaibab_meta %>%
  mutate(deploy_date = ymd(paste(cyear, deploy.month, deploy.day, sep = "-")),
         collect_date = ifelse(!is.na(visit2.month), paste(cyear, visit2.month, visit2.day, sep = "-"),
                                                         paste(cyear, visit1.month, visit1.day, sep = "-")),
         collect_date = ymd(collect_date),
         days_deployed = collect_date - deploy_date) %>%
  select(-c(cyear, deploy.month, deploy.day,  visit2.month, visit2.day))



kaibab_detect <- val_occ_day %>%
  filter(forest == "kaibab") %>%
  select(-forest) %>%
  left_join(kaibab_meta)






gila_meta <- read.csv(here::here("data/bioacoustics.2022.Gila.metadata.csv"))
