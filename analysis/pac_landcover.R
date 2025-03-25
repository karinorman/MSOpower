### Explore landcover types in existing PAC's ###

library(dplyr)
library(terra)
library(tidyterra)

pacs <- vect(here::here("data/MSO_PACs/MSO_PACs.shp")) %>%
  project("epsg:5070")

landcover <- rast(here::here("data/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif"))

pacs_landcovers <- extract(landcover, pacs, fun = table)

pac_lc_count <- pacs_landcovers %>%
  select(-`Fill-NoData`) %>%
  pivot_longer(-ID, names_to = "landcover", values_to = "count") %>%
  select(-ID) %>%
  group_by(landcover) %>%
  summarize(count = sum(count)) %>%
  filter(count != 0)


pac_percent <- pacs_landcovers %>%
  mutate(ID = as.character(ID)) %>%
  select(c(ID, unique(pac_lc_count$landcover))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(-c(ID, sum), ~ .x / sum)) %>%
  select(-sum) %>%
  pivot_longer(-ID, names_to = "landcover", values_to = "percent") %>%
  filter(percent != 0)


## average percentage across pacs of landcover types

mean_lc_percent <- pac_percent %>%
  select(-ID) %>%
  group_by(landcover) %>%
  summarize(mean = mean(percent))
