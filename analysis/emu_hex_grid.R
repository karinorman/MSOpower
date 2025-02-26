#####################################
######## FIT GRID TO EMU'S ##########
#### Get number of cells per EMU ####
#####################################

library(dplyr)
library(terra)
library(tidyterra)
library(ggplot2)
library(sf)
library(dggridR)

emus <- st_read(here::here("data/MSO_EMUs/MSO_EMUs.shp")) %>%
  select(UNIT) %>%
  # equal area projection
  st_transform(crs = 9822)


# get 4 km2 grid, make hex ID
st_grid <- st_make_grid(emus, cellsize =  2150, square = FALSE) %>%
  st_as_sf() %>%
  mutate(ID = row_number())

# this  doesn't give full grids cells on the edges, but it does assign emu lables to grid cells
grid_labels <-  st_grid %>% st_intersection(., emus)
# this one does!
emu_grid <- st_grid[emus,] %>%
  left_join(st_drop_geometry(grid_labels), by = "ID")

# check that they're all the same size
hex_sizes <- as.vector(st_area(grid_labels)) / 1000000

# get count of number
emu_grid %>%
  as.data.frame() %>%
  count(UNIT)


# zoom in on a little piece to see that the hexagons aren't distorted or cropped at the boundary
# example emu
#gila_emu <- emus %>% filter(UNIT == "Upper Gila Mountains")
# cropping bbox
#bb <- st_bbox(c(xmin = -112, xmax = -111, ymax = 35, ymin = 34.5))
bb <- st_bbox(c(xmin = -112, xmax = -111, ymax = 33, ymin = 31))

emu_grid %>%
  st_transform(crs = 4326) %>%
  st_crop(bb) %>%
  ggplot() +
  geom_sf(color = "grey", fill = "transparent") +
  geom_sf(data = emus %>%
            st_transform(crs = 4326) %>% st_crop(bb), fill = "transparent") #+
#st_transform(crs = 4326), fill = "transparent") #+
  #coord_sf(datum = 9822)

#plot grid on top of emu's
grid_labels %>%
  st_transform(crs = 4326) %>%
  ggplot() +
  geom_sf(color = "grey", fill = "transparent") +
  geom_sf(data = emus %>%
            st_transform(crs = 4326), fill = "transparent")



###############################################################
### Extract grid cell attributes from external data sources ###
###############################################################

# project the emus to the GFC raster projection
emus_proj <- emus %>% st_transform(crs = 4326)
emu_grid_proj <- emu_grid %>% st_transform(crs = 4326)
grid_bbox <- st_bbox(emu_grid_proj)


### Get the tree cover from year 2000 ###
if(!file.exists(here::here('data/gfc_treecover2000_study_area.tif'))){

  tile_files <- list.files(here::here("data/hanson_global_forest_change/"), pattern = ".tif", full.names = TRUE)

  tile_list <- purrr::map(tile_files, ~ rast(.x) %>% crop(., grid_bbox))

  treecover <- mosaic(sprc(tile_list))

  writeRaster(treecover, here::here('data/gfc_treecover2000_study_area.tif'))

} else{

  treecover <- rast(here::here('data/gfc_treecover2000_study_area.tif'))

}

# # plot to make sure bbox isn't too stringent
# bb <- bb <- st_bbox(c(xmin = -114.76548, xmax = -114, ymax = 37, ymin = 35))
# ggplot() +
#   geom_spatraster(data = treecover %>% crop(bb)) +
#   geom_spatvector(data = vect(emus_proj) %>% crop(bb), color = "grey", fill = "transparent") +
#   geom_spatvector(data = vect(emu_grid_proj) %>% crop(bb), color = "white", fill = "transparent")


# get grid level summaries of both percent tree cover
emu_grid_vect <- vect(emu_grid_proj)
grid_treecover_zonal <- zonal(treecover, emu_grid_vect, fun = "mean", na.rm = TRUE, as.polygons = TRUE)

### Get the landcover data ###
landcover_type <- rast(here::here("data/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif")) %>%
  crop(., emu_grid %>% st_transform(crs = 5070) %>% st_bbox) %>%
  project("epsg:4326")

# get metadata for landcover types
landcover_meta <- read.csv(here::here("data/LF2023_EVT_240_CONUS/CSV_Data/LF23_EVT_240.csv"))

# pull out landcover types and add back in a bunch of meta data, exclude some current landcover groups
grid_attr <- zonal(landcover_type, grid_treecover_zonal, fun = "modal", na.rm = TRUE, as.polygons = TRUE) %>%
  left_join(grid_attr %>% rename(LFRDB = EVT_NAME), landcover_meta)

#writeVector(grid_attr, here::here("data/grid_attr.shp"))


# Let's look at the living map stuff
lm_2022 <- rast(here::here("data/MSO_SDM_SDM_2022.tif")) %>%
  project("epsg:4326") %>%
  crop(vect(emus_proj))

lm_2000 <- rast(here::here("data/MSO_SDM_SDM_2000.tif")) %>%
  project("epsg:4326") %>%
  crop(vect(emus_proj))

habitat_rast <- c(lm_2022 %>% rename(habitat_2022 = SDM), lm_2000 %>% rename(habitat_2000 = SDM))

grid_attr_habitat <- zonal(habitat_rast, grid_attr, fun = "mean", na.rm = TRUE, as.polygons = TRUE)

# let's get values within grid so we can be a little more precise with what we want than simply mean
mso_poly_values <- extract(round(lm_2022), grid_attr, fun = table)
poly_ids <- grid_attr %>% pull(ID)

write.csv(mso_poly_values, here::here("data/mso_poly_values.csv"))

total_cells <- mso_poly_values %>%
  select(-ID) %>%
  rowSums(na.rm = TRUE) %>%
  data.frame(poly_id = poly_ids, extract_ID = mso_poly_values$ID, total_cell = .)

# get percent of each polygon above .5 cut off
threshold_count <- mso_poly_values %>%
  pivot_longer(-ID, names_to = "mso_value") %>%
  mutate(mso_value = as.numeric(mso_value)) %>%
  filter(mso_value >= 5000) %>%
  group_by(ID) %>%
  summarize(threshold_cells = sum(value)) %>%
  left_join(total_cells, by = c("ID" = "extract_ID")) %>%
  rename(extract_ID = ID) %>%
  mutate(percent = threshold_cells/total_cell) %>%
  filter(percent > 0.1)

grid_attr_habitat <- grid_attr_habitat %>%
  mutate(mso_percent_habitat = ifelse(ID %in% threshold_count$poly_id, 1, NA))
#####################################################
### Which hexes are included in our sample frame ####
#####################################################

grid_attr_habitat <- grid_attr_habitat %>%
  # initial pass, create two variables that give veg type and veg group if the grid is in the sample frame
   mutate(sample_frame_veg = ifelse(EVT_LF == "Tree" | Hansen_GFC > 10, EVT_NAME, NA),
          sample_frame_type = ifelse(EVT_LF == "Tree" | Hansen_GFC > 10, EVT_LF, NA)) %>%
   mutate(across(starts_with("sample_frame"), ~replace(., sample_frame_type %in%  c("Agriculture", "Barren", "Developed", "Sparse", "Water", "Snow-Ice"), NA))) %>%
  # Identify tree veg types, and whether or not they're included in the sample frame
  mutate(mso_habitat_type = case_when(
    sample_frame_veg == "Madrean Pinyon-Juniper Woodland" ~ "no",
    sample_frame_veg == "Madrean Encinal" ~ "no",
    sample_frame_veg == "Interior West Ruderal Riparian Forest" ~ "no",
    sample_frame_veg == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "yes",
    sample_frame_veg == "Western Warm Temperate Orchard" ~ "no",
    sample_frame_veg == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "yes",
    sample_frame_veg == "North American Warm Desert Riparian Woodland" ~ "no",
    sample_frame_veg == "North American Warm Desert Riparian Mesquite Bosque Woodland" ~ "no",
    sample_frame_veg == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "yes",
    sample_frame_veg == "Madrean Upper Montane Conifer-Oak Forest and Woodland" ~ "yes",
    sample_frame_veg == "North American Warm Desert Lower Montane Riparian Woodland" ~ "no",
    sample_frame_veg == "Madrean Juniper Savanna" ~ "no",
    sample_frame_veg == "Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland" ~ "yes",
    sample_frame_veg == "Colorado Plateau Pinyon-Juniper Woodland" ~ "no",
    sample_frame_veg == "Rocky Mountain Aspen Forest and Woodland" ~ "yes",
    sample_frame_veg == "Inter-Mountain Basins Aspen-Mixed Conifer Forest and Woodland" ~ "yes",
    sample_frame_veg == "Southern Rocky Mountain Juniper Woodland and Savanna" ~ "no",
    sample_frame_veg == "Southern Rocky Mountain Ponderosa Pine Savanna" ~ "yes",
    sample_frame_veg == "Inter-Mountain Basins Subalpine Limber-Bristlecone Pine Woodland" ~ "yes",
    sample_frame_veg == "Inter-Mountain Basins Curl-leaf Mountain Mahogany Woodland" ~ "no",
    sample_frame_veg == "Rocky Mountain Bigtooth Maple Ravine Woodland" ~ "yes",
    sample_frame_veg == "Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland" ~ "yes",
    sample_frame_veg == "Rocky Mountain Subalpine-Montane Riparian Woodland" ~ "yes",
    sample_frame_veg == "Western Cool Temperate Urban Evergreen Forest" ~ "no",
    sample_frame_veg == "Rocky Mountain Lodgepole Pine Forest" ~ "yes",
    sample_frame_veg == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "yes",
    sample_frame_veg == "Western Cool Temperate Orchard" ~ "no",
    sample_frame_veg == "Rocky Mountain Lodgepole Pine Forest" ~ "yes",
    sample_frame_veg == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "yes",
    sample_frame_veg == "Western Cool Temperate Orchard" ~ "no",
    sample_frame_veg == "Rocky Mountain Foothill Limber Pine-Juniper Woodland" ~ "no",
    sample_frame_veg == "Rocky Mountain Lower Montane-Foothill Riparian Woodland" ~ "no",
    sample_frame_veg == "Western Great Plains Riparian Woodland" ~ "no",
    sample_frame_veg == "Southern Rocky Mountain Mesic Montane Mixed Conifer Forest and Woodland" ~ "yes",
    sample_frame_veg == "Great Basin Pinyon-Juniper Woodland" ~ "no",
    sample_frame_veg == "Southern Rocky Mountain Pinyon-Juniper Woodland" ~ "no",
    .default = NA
  ))

writeVector(grid_attr_habitat, here::here("data/grid_attr_habitat.shp"))

#### Plot predicted MSO habitat against forested areas ####

habitat_2000_map <- ggplot() +
  geom_spatvector(data = emus_proj, color = "black", fill = "transparent") +
  geom_spatvector(data = grid_attr_habitat %>%
                    filter(!is.na(sample_frame_type)),
                  aes(color = sample_frame_type, fill = sample_frame_type)) +
  geom_spatvector(data = grid_attr_habitat %>% filter(habitat_2000 > 4000),
                  color = "red", fill = "red", alpha = 0.5) +
  scale_fill_discrete(na.value = "transparent") +
  theme_void()# +
#theme(legend.position = "none")
ggsave("figures/habitat_2000.jpeg", habitat_2000_map)


habitat_2022_map <-  ggplot() +
  geom_spatvector(data = emus_proj, color = "black", fill = "transparent") +
  geom_spatvector(data = grid_attr_habitat %>%
                    filter(!is.na(sample_frame_type)),
                  aes(color = sample_frame_type, fill = sample_frame_type)) +
  geom_spatvector(data = grid_attr_habitat %>% filter(habitat_2022 > 4000),
                  color = "red", fill = "red", alpha = 0.5) +
  scale_fill_discrete(na.value = "transparent") +
  theme_void()# +
#theme(legend.position = "none")
ggsave("figures/habitat_2022.jpeg", habitat_2022_map)

### limit the study frame to the types of forest MSO could concievably be in
ggplot() +
  geom_spatvector(data = emus_proj, color = "black", fill = "transparent") +
  geom_spatvector(data = grid_attr_habitat %>%
                    filter(mso_habitat_type == "yes"), color = 'grey') +
  geom_spatvector(data = grid_attr_habitat %>% filter(habitat_2000 > 4000),
                  color = "red", fill = "red", alpha = 0.5) +
  scale_fill_discrete(na.value = "transparent") +
  theme_void()# +

###############################################
### Let's figure out the area of each patch ###
###############################################

sample_grids <- grid_attr_habitat %>%
  filter(mso_habitat_type == "yes")

sample_polys <- sample_grids %>%
  select(ID) %>%
  aggregate() %>%
  disagg()

id_map <- relate(centroids(sample_grids, inside=TRUE), sample_polys, "intersects", pairs=TRUE) %>%
  as.data.frame() %>%
  bind_cols(grid_id = sample_grids$ID) %>%
  select(-id.x)

sample_poly_area <- sample_polys %>%
  project("epsg:9822") %>%
  expanse(unit = "ha") %>%
  as.data.frame() %>%
  rename(poly_area = ".") %>%
  mutate(agg_poly_id = row_number()) %>%
  left_join(id_map, by = c("agg_poly_id" = "id.y"))

small_patch <- sample_poly_area %>%
  count(agg_poly_id) %>%
  filter(n < 5) %>%
  left_join(sample_poly_area)

# final sample frame,
grid_sample_frame <- sample_grids %>%
  mutate(include_patch = as.factor(ifelse(ID %in% small_patch$grid_id, "no", "yes"))) %>%
  select(ID, UNIT, veg_type_landfire = sample_frame, habitat_2000, habitat_2022, mso_habitat_type, mso_percent_habitat, include_patch)

writeVector(grid_sample_frame, here::here("data/grid_sample_frame.shp"))

### let's look at stuff ###
# map with smaller patches
pal <- list("yes" = "#82A6B1", "no" = "#2F394D")
ggplot() +
  geom_spatvector(data = emus_proj, color = "black", fill = "transparent") +
  geom_spatvector(data = grid_sample_frame, aes(fill = include_patch, color = include_patch)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  geom_spatvector(data = grid_attr_habitat %>% filter(mso_percent_habitat == 1),
                  color = "#BC4749", fill = "#BC4749", alpha = 0.5) +
  #scale_fill_discrete(na.value = "transparent") +
  theme_void()# +

# map with only included patches
gila_emu <- emus_proj %>% filter(UNIT == "Upper Gila Mountains")

ggplot() +
  geom_spatvector(data = gila_emu, color = "black", fill = "transparent") +
  geom_spatvector(data = grid_sample_frame %>% crop(gila_emu), aes(fill = include_patch, color = include_patch)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  #geom_spatvector(data = grid_sample_frame %>% filter(include_patch == "yes") %>% crop(gila_emu), fill = "#82A6B1", color =  "#82A6B1") +
  geom_spatvector(data = grid_attr_habitat %>% filter(mso_percent_habitat == 1) %>% crop(gila_emu),
                  color = "#BC4749", fill = "#BC4749", alpha = 0.5) +
  #scale_fill_discrete(na.value = "transparent") +
  theme_void()


##################################################
########### Count of veg types by EMU ############
##################################################

grid_sample_frame %>%
  as.data.frame() %>%
  group_by(UNIT, veg_type_landfire) %>%
  summarize(hex_num = n_distinct(ID)) %>%
  readr::write_csv(here::here("data/EMU_veg_types.csv"))
