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



##################################################
## Which hexes are included in our sample frame ##
##################################################

# project the emus to the GFC raster projection
emus_proj <- emus %>% st_transform(crs = 4326)
emu_grid_proj <- emu_grid %>% st_transform(crs = 4326)
grid_bbox <- st_bbox(emu_grid_proj)

if(!file.exists(here::here('data/gfc_treecover2000_study_area.tif'))){

  tile_files <- list.files(here::here("data/hanson_global_forest_change/"), pattern = ".tif", full.names = TRUE)

  tile_list <- purrr::map(tile_files, ~ rast(.x) %>% crop(., grid_bbox))

  treecover <- mosaic(sprc(tile_list))

  writeRaster(treecover, here::here('data/gfc_treecover2000_study_area.tif'))

} else{

  treecover <- rast(here::here('data/gfc_treecover2000_study_area.tif'))

}

# plot to make sure bbox isn't too stringent
bb <- bb <- st_bbox(c(xmin = -114.76548, xmax = -114, ymax = 37, ymin = 35))
ggplot() +
  geom_spatraster(data = treecover %>% crop(bb)) +
  geom_spatvector(data = vect(emus_proj) %>% crop(bb), color = "grey", fill = "transparent") +
  geom_spatvector(data = vect(emu_grid_proj) %>% crop(bb), color = "white", fill = "transparent")


# get grid level summaries of both percent tree cover and landcover type
emu_grid_vect <- vect(emu_grid_proj)
grid_treecover_zonal <- zonal(treecover, emu_grid_vect, fun = "mean", na.rm = TRUE, as.polygons = TRUE)

# read in data
landcover_type <- rast(here::here("data/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif")) %>%
  crop(., emu_grid %>% st_transform(crs = 5070) %>% st_bbox) %>%
  project("epsg:4326")

# get metadata for landcover types
landcover_meta <- read.csv(here::here("data/LF2023_EVT_240_CONUS/CSV_Data/LF23_EVT_240.csv"))

# pull out landcover types and add back in a bunch of meta data, exclude some current landcover groups
grid_attr <- zonal(landcover_type, grid_treecover_zonal, fun = "modal", na.rm = TRUE, as.polygons = TRUE) %>%
  left_join(grid_attr %>% rename(LFRDB = EVT_NAME), landcover_meta) %>%
  mutate(sample_frame = ifelse(EVT_LF == "Tree" | Hansen_GFC > 10, EVT_NAME, NA),
         sample_frame_type = ifelse(EVT_LF == "Tree" | Hansen_GFC > 10, EVT_LF, NA)) %>%
  mutate(across(starts_with("sample_frame"), ~replace(., sample_frame_type %in%  c("Agriculture", "Barren", "Developed", "Sparse", "Water", "Snow-Ice"), NA)))

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


##### Find MSO habitat types #####

grid_attr_habitat <- grid_attr_habitat %>%
  #filter(sample_frame_type == "Tree") %>%
  mutate(MSO_habitat = case_when(
    sample_frame == "Madrean Pinyon-Juniper Woodland" ~ "no",
    sample_frame == "Madrean Encinal" ~ "no",
    sample_frame == "Interior West Ruderal Riparian Forest" ~ "no",
    sample_frame == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "yes",
    sample_frame == "Western Warm Temperate Orchard" ~ "no",
    sample_frame == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "yes",
    sample_frame == "North American Warm Desert Riparian Woodland" ~ "no",
    sample_frame == "North American Warm Desert Riparian Mesquite Bosque Woodland" ~ "no",
    sample_frame == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "yes",
    sample_frame == "Madrean Upper Montane Conifer-Oak Forest and Woodland" ~ "yes",
    sample_frame == "North American Warm Desert Lower Montane Riparian Woodland" ~ "no",
    sample_frame == "Madrean Juniper Savanna" ~ "no",
    sample_frame == "Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland" ~ "yes",
    sample_frame == "Colorado Plateau Pinyon-Juniper Woodland" ~ "no",
    sample_frame == "Rocky Mountain Aspen Forest and Woodland" ~ "yes",
    sample_frame == "Inter-Mountain Basins Aspen-Mixed Conifer Forest and Woodland" ~ "yes",
    sample_frame == "Southern Rocky Mountain Juniper Woodland and Savanna" ~ "no",
    sample_frame == "Southern Rocky Mountain Ponderosa Pine Savanna" ~ "yes",
    sample_frame == "Inter-Mountain Basins Subalpine Limber-Bristlecone Pine Woodland" ~ "yes",
    sample_frame == "Inter-Mountain Basins Curl-leaf Mountain Mahogany Woodland" ~ "no",
    sample_frame == "Rocky Mountain Bigtooth Maple Ravine Woodland" ~ "yes",
    sample_frame == "Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland" ~ "yes",
    sample_frame == "Rocky Mountain Subalpine-Montane Riparian Woodland" ~ "yes",
    sample_frame == "Western Cool Temperate Urban Evergreen Forest" ~ "no",
    sample_frame == "Rocky Mountain Lodgepole Pine Forest" ~ "yes",
    sample_frame == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "yes",
    sample_frame == "Western Cool Temperate Orchard" ~ "no",
    sample_frame == "Rocky Mountain Lodgepole Pine Forest" ~ "yes",
    sample_frame == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "yes",
    sample_frame == "Western Cool Temperate Orchard" ~ "no",
    sample_frame == "Rocky Mountain Foothill Limber Pine-Juniper Woodland" ~ "no",
    sample_frame == "Rocky Mountain Lower Montane-Foothill Riparian Woodland" ~ "no",
    sample_frame == "Western Great Plains Riparian Woodland" ~ "no",
    sample_frame == "Southern Rocky Mountain Mesic Montane Mixed Conifer Forest and Woodland" ~ "yes",
    sample_frame == "Great Basin Pinyon-Juniper Woodland" ~ "no",
    sample_frame == "Southern Rocky Mountain Pinyon-Juniper Woodland" ~ "no",
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
                    filter(MSO_habitat == "yes"), color = 'grey') +
  geom_spatvector(data = grid_attr_habitat %>% filter(habitat_2000 > 4000),
                  color = "red", fill = "red", alpha = 0.5) +
  scale_fill_discrete(na.value = "transparent") +
  theme_void()# +

###############################################
### Let's figure out the area of each patch ###
###############################################

# sample_polys <- grid_attr_habitat %>% filter(!is.na(sample_frame)) %>% mutate(sample_area = 1, sample_area = 2) %>%
#   select(sample_area) %>%
#   aggregate() %>%
#   mutate(sample_area = 1)
#   #expanse()

# Get the grid cells in the sample frame that are adjacent
adj_grids <- grid_attr_habitat %>%
  filter(MSO_habitat == "yes") %>%
  adjacent(type = "rook")

adj_grids_emu <- grid_attr_habitat %>%
  filter(MSO_habitat == "yes", UNIT == "Basin & Range - West") %>%
  adjacent(type = "rook", symmetrical = TRUE)

emu_ids <-  grid_attr_habitat %>%
  as.data.frame() %>%
  filter(!is.na(sample_frame), UNIT == "Basin & Range - West") %>%
  select(ID) %>%
  mutate(adj_id = row_number())

sym_mat <- adj_grids_emu %>%
  as.data.frame() %>%
  mutate(ind_value = 1) %>%
  # fill in the missing index so we don't have suprious adjacencies in the matrix
  tidyr::complete(to = 1:dim(emu_ids)[1]) %>%
  arrange(from) %>%
  pivot_wider(names_from = to, values_from = ind_value, values_fill = 0) %>%
  tidyr::complete(from = 1:dim(emu_ids)[1])  %>%
  filter(!is.na(from)) %>%
  select(-from)

id_names <- as.integer(colnames(sym_mat)) %>% sort() %>% as.character()
sym_mat <- sym_mat[, id_names]

# adjacency matrix to raster, find clumps
sym_rast <- raster::raster(as.matrix(sym_mat))
clumps <- raster::clump(sym_rast, directions = 4)

# back to data frame, add ID's back in, find the clump id for each grid hex
clump_df <- as.data.frame(as.matrix(clumps))
colnames(clump_df) <- id_names
clump_df <- clump_df %>%
  janitor::remove_empty(which = "cols") %>%
  mutate(adj_id = as.integer(id_names)) %>%
  rowwise() %>%
  # there's only one clump ID per grid cell, so we can use max to get that unique value (unique() doesn't drop na's)
  mutate(clump_id = max(c_across(-adj_id), na.rm = TRUE)) %>%
  select(adj_id, clump_id) %>%
  filter(!is.infinite(clump_id)) %>%
  left_join(emu_ids, by = "adj_id")






clump_ids <- unique(as.vector(clump_mat), na.rm = TRUE)

tot <- max(clump_mat, na.rm = TRUE)
res <- vector("list", tot)
for (i in 1:tot){
  res[i] <- list(which(clump_mat == i, arr.ind = TRUE))
}

# Get count of number of cells in each adjacency group, get area of aggregated polygons






