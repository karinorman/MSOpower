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
# emu <- vect(here::here("data/MSO_EMUs/MSO_EMUs.shp")) %>%
#   aggregate() %>%
#   buffer(width = 32186.9)

##########################################
#### Get tigris road data from census ####
##########################################

# grid footprint with 30 mile buffer
grid_footprint <- grid_sample_frame %>%
  aggregate() %>%
  buffer(width = 48280.3)

# download roads from tigris
state_list <- c("08", "49", "35", "04", "48", "16", "56", "32", "06")

# get list of counties in our states of interest
county_shp <- counties(state = state_list) %>%
  st_transform(4326)

# get only the counties that intersect with our sample grid
grid_counties <- county_shp %>%
  st_intersects(st_as_sf(grid_footprint))

county_shp_grid <- county_shp %>%
  bind_cols(data.frame(intersects = apply(grid_counties, 1, any))) %>%
  filter(intersects == TRUE)

#ggplot() + geom_spatvector(data = vect(county_shp_grid)) + geom_spatvector(data = grid_footprint, color = "red", fill = "transparent")

# create folder to save out county level shapefiles
dir.create(here::here('data/road_shp'), showWarnings = FALSE)

# map across states and counties, project and save shapefile
road_sf_list <- purrr::pmap(county_shp_grid %>% select(STATEFP, COUNTYFP) %>% st_drop_geometry(), function(STATEFP, COUNTYFP) {
  road_shp <- roads(state = STATEFP, county = COUNTYFP) %>% st_transform(4326)

  st_write(road_shp, dsn = here::here('data/road_shp', paste(STATEFP, COUNTYFP, "road.shp", sep = "_")))
  })


# read in and combine all shapefiles
grid_sf <- st_as_sf(grid_sample_frame)
ordered_hex_ids <- grid_sample_frame$ID
file_list <- list.files(here::here('data/road_shp'), pattern = ".shp")

# master shapefile, then get name of closest road for each hex
master_roads <- purrr::map_dfr(file_list, ~st_read(here::here("data/road_shp", .x)))
closest_road <- st_nearest_feature(grid_sf, master_roads)

# map the road ID to it's string name from the original shapefile
road_df <- data.frame(hex_id = ordered_hex_ids, closest_road = closest_road) %>%
  left_join(master_roads %>%
              select(LINEARID) %>%
              mutate(closest_road = row_number()) %>%
              st_drop_geometry())

# get distance between each hex and its closest road
road_distance <- pmap_dfr(road_df %>% select(hex_id, road_id = LINEARID), function(hex_id, road_id){
  st_distance(grid_sf %>% filter(ID == hex_id), master_roads %>% filter(LINEARID == road_id)) %>%
    as.vector() %>%
    data.frame(hex_id = hex_id, distance = .)
})

##############################
#### Get dist to FS roads ####
##############################

# get min distance to road from hex boundary
nfs_roads <- vect(here::here("data/National_Forest_System_Roads_(Feature_Layer)/National_Forest_System_Roads_(Feature_Layer).shp")) %>%
  project("epsg:4326") %>%
  crop(grid_footprint) %>%
  st_as_sf()

nfs_closest_road <- st_nearest_feature(grid_sf, nfs_roads)

# map the road ID to it's string name from the original shapefile
nfs_road_df <- data.frame(hex_id = ordered_hex_ids, closest_road = nfs_closest_road) %>%
  left_join(nfs_roads %>%
              select(OBJECTID) %>%
              mutate(closest_road = row_number()) %>%
              st_drop_geometry())

# get distance between each hex and its closest road
nfs_road_distance <- pmap_dfr(nfs_road_df %>% select(hex_id, road_id = OBJECTID), function(hex_id, road_id){
  st_distance(grid_sf %>% filter(ID == hex_id), nfs_roads %>% filter(OBJECTID == road_id)) %>%
    as.vector() %>%
    data.frame(hex_id = hex_id, distance = .)
})

usethis::use_data(nfs_road_distance)

# get min distance between the two road data sources
min_distance <- road_distance %>%
  rename(tigris_dist = distance) %>%
  left_join(nfs_road_distance %>%
              rename(nfs_dist = distance)) %>% 
  rowwise() %>%
  mutate(road_distance = min(tigris_dist, nfs_dist)) %>%
  select(hex_id, road_distance)

hex_metadata <- sample_frame %>%
  rename(hex_id = ID) %>%
  left_join(min_distance)

readr::write_csv(hex_metadata, here::here("data/hex_metadata.csv"))

