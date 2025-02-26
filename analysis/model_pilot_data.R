library(dplyr)

##########################################
### Processing calls validated by Dana ###
##########################################

val_occ <- read.csv(here::here("data/top_conf_all_REID.csv")) %>%
  mutate(date = lubridate::mdy(date)) %>%
  rename(spow = "SPOW.") %>%
  filter(spow == "Y")

val_occ_day <- val_occ %>%
  select(forest, unit, date) %>%
  distinct()

# the gila has detection at 7 units, the kaibab at 6
val_occ_day %>% count(forest, unit)
