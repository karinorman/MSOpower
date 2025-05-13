library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(terra)
library(tigris)
library(tidyterra)
library(sf)

# read in sample frame shapefile
grid_sample_frame <- vect(here::here("data/grid_sample_frame.shp"))
names(grid_sample_frame) <- c("ID", "UNIT", "veg_type_landfire", "habitat_2000", "habitat_2022", "mso_habitat_type", "mso_percent_habitat", "include_patch")

cents <- centroids(grid_sample_frame)

# convert to dataframe, with centroid coordinates
sample_frame <- grid_sample_frame %>%
  filter(include_patch == "yes") %>%
  # get centroids of each hex
  centroids() %>%
  as.data.frame(geom = "XY") %>%
  # label landcovers as high or low occupancy
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
  )) %>%
  mutate(emu = case_when(
    UNIT == "Basin & Range - East" ~ "BRE",
    UNIT == "Basin & Range - West" ~ "BRW",
    UNIT == "Colorado Plateau" ~ "CP",
    UNIT == "Southern Rocky Mountains" ~ "SRM",
    UNIT == "Upper Gila Mountains" ~ "UGM"
  )) %>%
  select(ID, emu, unit = UNIT, veg_type_landfire, occupancy)

sf_use_s2(FALSE)

# emu footprint with 20 mile buffer
emu <- vect(here::here("data/MSO_EMUs/MSO_EMUs.shp")) %>%
  aggregate() %>%
  buffer(width = 32186.9)

# grid footprint with 30 mile buffer
grid_footprint <- grid_sample_frame %>%
  aggregate() %>%
  buffer(width = 48280.3)

# download roads from tigris
state_list <- c("08", "49", "35", "04", "48", "16", "56", "32", "06")

county_shp <- counties(state = state_list) %>%
  st_transform(4326)

grid_counties <- county_shp %>%
  st_intersects(st_as_sf(grid_footprint))

county_shp_grid <- county_shp %>%
  bind_cols(data.frame(intersects = apply(grid_counties, 1, any))) %>%
  filter(intersects == TRUE)

ggplot() + geom_spatvector(data = vect(county_shp_grid)) + geom_spatvector(data = grid_footprint, color = "red", fill = "transparent")

dir.create(here::here('data/road_shp'), showWarnings = FALSE)

# map across states and pull out a list of county codes from the county shapefiles
road_sf_list <- purrr::pmap(county_shp_grid %>% select(STATEFP, COUNTYFP) %>% st_drop_geometry(), function(STATEFP, COUNTYFP) {
  road_shp <- roads(state = STATEFP, county = COUNTYFP) %>% st_transform(4326)

  st_write(road_shp, dsn = here::here('data/road_shp', paste(STATEFP, COUNTYFP, "road.shp", sep = "_")))
  })


grid_sf <- st_as_sf(grid_sample_frame)
ordered_hex_ids <- grid_sample_frame$ID
file_list <- list.files(here::here('data/road_shp'), pattern = ".shp")

# min_road_dist <- map_dfr(file_list, ~ st_read(here::here('data/road_shp', .x)) %>%
#       st_distance(grid_sf, .) %>%
#       apply(MARGIN = 1, FUN = min) %>%
#       data.frame(hex_id = ordered_hex_ids, road_dist = ., county = .x)
#       )


master_roads <- purrr::map_dfr(file_list, ~st_read(here::here("data/road_shp", .x)))
closest_road <- st_nearest_feature(grid_sf, master_roads)

road_df <- data.frame(hex_id = ordered_hex_ids, closest_road = closest_road) %>%
  left_join(master_roads %>%
              select(LINEARID) %>%
              mutate(closest_road = row_number()) %>%
              st_drop_geometry()) %>%
  rowwise() %>%
  mutate(distance = st_distance(grid_sf %>% filter(ID == hex_id), master_roads %>% filter(LINEARID == LINEARID)))

road_distance <- purrr::pmap(road_df %>% select(hex_id, LINEARID), function(hex_id, LINEARID){
  st_distance(grid_sf %>% filter(ID == hex_id), master_roads %>% filter(LINEARID == LINEARID)) %>%
    as.vector() %>%
    data.frame(hex_id = hex_id, distance = .)
})



test_shp <- st_read(here::here('data/road_shp', file_list[1]))

test_dist <- st_distance(st_as_sf(grid_sample_frame), test_shp)

# tigris_road_sf <- bind_rows(road_sf_list) %>%
#   st_transform(4326)

# get min distance
# tigris_road_dist <- purrr::map_dfr(unique(tigris_road_sf$RTTYP), ~st_distance(st_as_sf(grid_sample_frame), tigris_road_sf %>%
#                                                                             filter(RTTYP == .x)))



county_chunk_min_dist <- map_dfr(road_sf_list, ~st_distance(st_as_sf(grid_sample_frame), .x) %>%
              apply(MARGIN = 1, FUN = min) %>%
              data.frame(hex_id = ordered_hex_ids, road_dist = .))

## Let's also get distance to FS roads, in case that data product is more up to date
# get boundary for states
boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Arizona","New Mexico", "Utah", "Colorado")) %>%
  aggregate()

# get min distance to road from hex boundary
nfs_roads <- vect(here::here("data/National_Forest_System_Roads_(Feature_Layer)/National_Forest_System_Roads_(Feature_Layer).shp")) %>%
  project("epsg:4326") %>%
  crop(boundary_states)

nfs_road_dist <- st_distance(st_as_sf(grid_sample_frame), st_as_sf(nfs_roads))

min_nfs_road <- min_road_dist %>%
  apply(MARGIN = 1, FUN = min) %>%
  data.frame(hex_id = grid_sample_frame$ID, road_dist = .)

far_roads <- min_road %>%
  filter(road_dist > 200000) %>%
  pull(hex_id)

ggplot() +
  geom_spatvector(data = roads) +
  geom_spatvector(data = grid_sample_frame %>% filter(ID %in% far_roads) %>% centroids(), color = "red", fill = "red")



# buffer_hex <- grid_sample_frame %>%
#   filter(ID == 8786) %>%
#   buffer(width = 2000)
#
# tigris_buffer <- tigris_road_sf %>% vect() %>% crop(buffer_hex)
