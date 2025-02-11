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

bb <- st_bbox(c(xmin = -112, xmax = -111, ymax = 35, ymin = 34.5))
gila_emu <- emus %>% filter(UNIT == "Upper Gila Mountains") %>%
  # let's get a tiny portion of it so we can zoom in
  st_crop(bb)

# construct a grid with resolution close to desired area
grid <- dgconstruct(area = 4, metric = TRUE, topology = "HEXAGON")
# snap it to the
emus_grid <- dgshptogrid(grid, gila_emu) #%>% vect()
#cellcenters <- dgSEQNUM_to_GEO(emus_grid, in_seqnum)

# I want to work in terra please
#emus <- vect(emus)

grid_by_emus <- emus_grid %>% terra::intersect(emus)

# Visualize grid
ggplot() +
  geom_sf(data = emus  %>% filter(UNIT == "Upper Gila Mountains")) +
  geom_sf(data = emus_grid %>% filter(UNIT == "Upper Gila Mountains"))#, fill= "grey", color= "white", alpha = 0.4)



# Let's generate some random unit locations
emus_bbox <- st_bbox(gila_emu) %>%
  st_as_sfc() %>%
  st_sf()

pts <- st_sample(emus_bbox, 500) %>%
  st_sf(as.data.frame(st_coordinates(.)), geometry = .) %>%
  rename(lat = Y, lon = X)
pts$cell <- dgGEO_to_SEQNUM(emus_grid, pts$lon, pts$lat)$seqnum


hexagons <- dgcellstogrid(emus_grid, unique(pts$cell)) %>%
  st_as_sf()

ggplot() + geom_sf(data = hexagons) +
  geom_sf(data = gila_emu) +
  geom_sf(data = pts, size = 0.5) +
  theme_bw()

st_grid %>% st_as_sf() %>% st_crop(bb) %>% ggplot() + geom_sf() + geom_sf(data = gila_emu, fill = "transparent")


# try the st approach
# original projection
#st_grid <- st_make_grid(emus, cellsize =  0.021, square = FALSE)
st_grid <- st_make_grid(emus, cellsize =  2150, square = FALSE)

hex_sizes <- as.vector(st_area(st_grid)) / 1000000

