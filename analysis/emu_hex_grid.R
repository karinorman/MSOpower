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

  # treecover_tile1 <- rast(here::here("data/hanson_global_forest_change/Hansen_GFC-2023-v1.11_treecover2000_40N_120W.tif")) %>%
  #   crop(., emus_bbox)
  #
  # treecover_tile2 <- rast(here::here("data/hanson_global_forest_change/Hansen_GFC-2023-v1.11_treecover2000_40N_110W.tif")) %>%
  #   crop(., emus_bbox)
  #
  # treecover <- mosaic(treecover_tile1, treecover_tile2)

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

emu_grid_vect <- vect(emu_grid_proj)

grid_forest_cover_zonal <- zonal(treecover, emu_grid_vect, fun = "mean", na.rm = TRUE, as.polygons = TRUE)
#grid_forest_cover <- extract(treecover, emu_grid_vect)

