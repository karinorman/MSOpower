library(dplyr)
library(terra)
library(sf)
library(tidyterra)
library(geosphere)
library(ggplot2)

############################################################################
## Let's make an example figure showing hexes with samples in those hexes ##
############################################################################


# ex_emu <- st_read(here::here("data/MSO_EMUs/MSO_EMUs.shp")) %>%
#   select(UNIT) %>%
#   filter(UNIT == "Colorado Plateau") %>%
#   # equal area projection
#   st_transform(crs = 9822)


# get 4 km2 grid, make hex ID
# st_grid <- st_make_grid(ex_emu, cellsize =  2150, square = FALSE) %>%
#   st_as_sf() %>%
#   mutate(ID = row_number())


# Let's make a square example bounding box
# Create midpoint
midpoint <-  st_point(c(-5697609, 6434438)) %>%
  st_sfc() %>%
  st_set_crs(9822)
# make a circular buffer w/5km radius
buf_round <- st_buffer(midpoint, dist = 10000)
# bounding box (square of the round buffer)
box_buf <- st_bbox(buf_round) %>%
  st_as_sfc() %>%
  st_as_sf()

# Length of one side of the box
box_buf %>%
  st_cast('LINESTRING') %>%
  st_length() / 4

#Make grid for box
st_grid <- st_make_grid(box_buf, cellsize =  2150, square = FALSE) %>%
  st_as_sf() %>%
  mutate(ID = row_number(),
         sampled = ifelse(ID %in% c(39,106,72,32), TRUE, FALSE),
         sample_frame = ifelse(ID %in% c(1:54, 56:60, 64:66, 71:72, 76:78,
                                         82:84, 87:90, 93:96, 99:115),
                               TRUE, FALSE),
         sample_frame = ifelse(ID %in% c(24, 29, 30, 36, 42, 48, 41, 102, 114, 108),
                               FALSE, sample_frame),
         fill_var = case_when(
           sample_frame == TRUE & sampled == TRUE ~ "sampled",
           sample_frame == TRUE & sampled == FALSE ~ "included",
           .default = "excluded"
         )) %>%
  # make it more rectangular
  filter(!ID %in% c(1:18, 115:126))# %>%
  # st_crop(box_buf)

# set.seed(452)
# points <- st_sample(st_grid %>% filter(sampled == TRUE), size = c(3,3,3,3), type = "random")

# Get centroids of sample hexes
points <- st_centroid(st_grid %>% filter(sampled == TRUE))
st_geometry(points) <- "geometry"
#point_buffers <- st_buffer(points, dist = 600)
points_proj <- points %>% st_transform(crs = 4326)

# Get points in the angled cardinal directions
dist_points <- purrr::map_dfr(#c(45, 135, 225, 315),
  c(60, 150, 240, 330),
                                   ~destPoint(st_coordinates(points_proj), .x, 600) %>%
  as.data.frame() %>%
    bind_cols(., ID = points_proj$ID) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 9822)) %>%
  rbind(., points %>% select(geometry, ID)) %>%
  mutate(point_id = row_number())
  #rbind(., points %>% st_geometry() %>% st_as_sf() %>% rename(x = geometry))

set.seed(432)
# get "sampled" points
points_samp <- dist_points %>%
  group_by(ID) %>%
  slice_sample(n = 3) %>%
  mutate(aru = TRUE) #%>%
  #rbind(., dist_points %>% filter(!point_id %in% .$point_id) %>% mutate(aru = FALSE))

points_plot <- points_samp %>%
  rbind(., dist_points %>% filter(!point_id %in% points_samp$point_id) %>% mutate(aru = FALSE))

# plot(st_geometry(st_grid))
# plot(st_geometry(points), pch = 20, add = TRUE)
# plot(st_geometry(dist_points), pch = 20, add = TRUE)

##### Plot hexes with ID labels
hex_cent <- st_centroid(st_grid) %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

sampling_example <- ggplot() +
  geom_spatvector(data = vect(st_grid), aes(fill = fill_var)) +
  geom_spatvector(data = points_plot, aes(shape = aru), size = 3.5) +
  #geom_text(data = hex_cent, aes(label = ID, x = lon, y = lat)) +
  ggthemes::theme_map() +
  scale_fill_manual(values = list("excluded" = "transparent", "included" = "#dbe5d3", "sampled" = "#87A96B")) +
  scale_shape_manual(values = list("TRUE" = 16, "FALSE" = 1)) +
  theme(legend.position = "none")

ggsave(here::here("figures/aru_sampling_example.jpeg"), sampling_example, height = 10.5, width = 10.3)
