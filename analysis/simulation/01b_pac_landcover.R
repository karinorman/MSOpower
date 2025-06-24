### Explore landcover types in existing PAC's ###

library(dplyr)
library(terra)
library(tidyterra)
library(sf)

pacs <- vect(here::here("data/MSO_PACs/MSO_PACs.shp")) %>%
  project("epsg:5070")

landcover <- rast(here::here("data/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif"))

pacs_landcovers <- extract(landcover, pacs, fun = table)# %>%
  #bind_cols(pacs$PAC_ID)

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


## let's look at the pac area that isn't included in the sample frame
grid_sample_frame <- vect(here::here("data/grid_sample_frame.shp")) %>% project("epsg:5070")
names(grid_sample_frame) <- c("ID", "UNIT", "veg_type_landfire", "habitat_2000", "habitat_2022", "mso_habitat_type", "mso_percent_habitat", "include_patch")

#pacs <- pacs %>% project("epsg:4326")

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

veg_grid <- grid_sample_frame %>% left_join(emu_veg)

sample_grid_footprint <- veg_grid %>%
  #filter(occupancy == "high") %>%
  mutate(frame = 1) %>%
  #select(occupancy) %>%
  select(frame) %>%
  aggregate()

sf_use_s2(FALSE)
pacs_sf <- st_as_sf(pacs)
grid_footprint_sf <- st_as_sf(sample_grid_footprint)

grid_pacs_int <- pacs_sf %>%
  mutate(intersects_grid = st_intersects(pacs_sf, grid_footprint_sf) %>% as.matrix() %>% as.vector())

# what percent of PAC's are included in the grid
dim(grid_pacs_int %>% filter(intersects_grid == FALSE))[1]/dim(grid_pacs_int)[1]


# Let's look at the PACs not contained in the sampling frame, and what their landcover types are
excluded_pacs <- grid_pacs_int %>%
  dplyr::filter(intersects_grid == FALSE) %>%
  vect()


exclude_pacs_lc <- extract(landcover, excluded_pacs, fun = table)

exclude_pac_percent <- exclude_pacs_lc %>%
  mutate(ID = as.character(ID)) %>%
  select(c(ID, unique(pac_lc_count$landcover))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(-c(ID, sum), ~ .x / sum)) %>%
  select(-sum) %>%
  pivot_longer(-ID, names_to = "landcover", values_to = "percent") %>%
  filter(percent != 0)
