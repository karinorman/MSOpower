library(dplyr)
library(ggplot2)

# Create example image of occupancy process and simulation

# Occupancy grid
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
set.seed(312415)
st_grid <- st_make_grid(box_buf, cellsize =  2150, square = FALSE) %>%
  st_as_sf() %>%
  # get the hexes that are occupied
  mutate(ID = row_number(),
         occupied = ifelse(ID %in% sample(ID, 20), TRUE, FALSE)) %>%
  # get the hexes that are sampled
  mutate(sampled = ifelse(ID %in% sample(ID, 34), TRUE, FALSE)) %>%
  # make it more rectangular
  filter(!ID %in% c(1:18, 115:126))#

occ_fig <- ggplot() +
  geom_spatvector(data = vect(st_grid), aes(fill = occupied)) +
  #geom_text(data = hex_cent, aes(label = ID, x = lon, y = lat)) +
  ggthemes::theme_map() +
  scale_fill_manual(values = list("FALSE" = "transparent", "TRUE" = "#d88993")) +
  scale_shape_manual(values = list("TRUE" = 16, "FALSE" = 1)) +
  theme(legend.position = "none")

ggsave(here::here("figures/occ_grid_example.jpg"), occ_fig, height = 10.5, width = 10.3)

# now add the sampling process on top
occ_samp_fig <- ggplot() +
  geom_spatvector(data = vect(st_grid), aes(fill = occupied, linewidth = sampled)) +
  #geom_text(data = hex_cent, aes(label = ID, x = lon, y = lat)) +
  ggthemes::theme_map() +
  scale_fill_manual(values = list("FALSE" = "transparent", "TRUE" = "#d88993")) +
  scale_linewidth_manual(values = c(0.5,2)) +
  theme(legend.position = "none")

ggsave(here::here("figures/occ_grid_samp_example.jpg"), occ_samp_fig, height = 10.5, width = 10.3)

## Trend example
ggplot() +
  # scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  geom_segment(aes(x = 1, xend = 10, y = 1, yend = 0.75), color = "#28587B", size = 1) +
  theme_classic() +
  ylab("Occupancy, \u03A8") +
  xlab("Year") +
  ylim(c(0.7, 1))

ggsave(here::here("figures/ex_25decline.jpg"), height = 2, width = 4)

## Power curve example

ggplot() +
  scale_x_continuous(expand=c(0,0), limits = c(0, 500) ) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  #geom_segment(aes(x = 1, xend = 10, y = 1, yend = 0.75), color = "#28587B", size = 1) +
  theme_classic() +
  ylab("Power") +
  xlab("Sample Size")

ggsave(here::here("figures/power_template.jpg"), height = 2, width = 4)
