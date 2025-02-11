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
  st_transform(crs = 9822)
  # project("epsg:4326") %>%
  # sf::st_as_sf(v)

gila_emu <- emus %>% filter(UNIT == "Upper Gila Mountains")

# try the st approach
# original projection
#st_grid <- st_make_grid(emus, cellsize =  0.021, square = FALSE)
st_grid <- st_make_grid(emus, cellsize =  2150, square = FALSE)# %>%

# check that they're all the same size
hex_sizes <- as.vector(st_area(st_grid)) / 1000000

# zoom in on a little piece to see that the hexagons aren't distorted
bb <- st_bbox(c(xmin = -112, xmax = -111, ymax = 35, ymin = 34.5))
st_grid %>%
  st_intersection(gila_emu) %>%
  st_transform(crs = 4326) %>%
  st_crop(bb) %>%
  ggplot() +
  geom_sf(color = "grey", fill = "transparent") +
  geom_sf(data = gila_emu %>%
            st_transform(crs = 4326) %>% st_crop(bb), fill = "transparent") #+
  #coord_sf(datum = 9822)


emu_grids <- st_grid %>% st_as_sf() %>% mutate(ID = row_number()) %>% st_intersection(., emus)

# get count of number
emu_grids %>%
  as.data.frame() %>%
  count(UNIT)


#### Attempt at dggrid version
# # construct a grid with resolution close to desired area
# grid <- dgconstruct(area = 4, metric = TRUE, topology = "HEXAGON")
# # snap it to the
# emus_grid <- dgshptogrid(grid, gila_emu) #%>% vect()
# #cellcenters <- dgSEQNUM_to_GEO(emus_grid, in_seqnum)
#
# # I want to work in terra please
# #emus <- vect(emus)
#
# grid_by_emus <- emus_grid %>% terra::intersect(emus)
#
# # Visualize grid
# ggplot() +
#   geom_sf(data = emus  %>% filter(UNIT == "Upper Gila Mountains")) +
#   geom_sf(data = emus_grid %>% filter(UNIT == "Upper Gila Mountains"))#, fill= "grey", color= "white", alpha = 0.4)
