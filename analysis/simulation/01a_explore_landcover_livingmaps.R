library(dplyr)
library(terra)
library(tidyterra)
library(ggplot2)
library(sf)

### read in the data

emus <- st_read(here::here("data/MSO_EMUs/MSO_EMUs.shp")) %>%
  select(UNIT) %>%
  # equal area projection
  st_transform(crs = 9822)

emus_proj <- emus %>% st_transform(crs = 4326)

# get hex level information
grid_attr <- vect(here::here("data/grid_attr_habitat.shp"))
landcover_meta <- read.csv(here::here("data/LF2023_EVT_240_CONUS/CSV_Data/LF23_EVT_240.csv"))

# get only the hexes that were id'd by living maps
lm_hexes <- grid_attr %>% filter(mso_perce0 == 1) %>% select(ID, mso_perce0)

# extract landcover values in the hexes
lm_landcover <- extract(landcover_type, lm_hexes, fun = "table")

# get total counts across hexes for landcover types
lm_landcover_counts <- lm_landcover %>%
  bind_cols(poly_id = lm_hexes$ID) %>%
  select(-`Fill-NoData`) %>%
  pivot_longer(cols = -c("ID", "poly_id"), names_to = "landcover", values_to = "count") %>%
  filter(count != 0) %>%
  select(-poly_id) %>%
  group_by(landcover) %>%
  summarize(count = sum(count))

# landcover raster cropped to the living map hexes
habitat_landcover_rast <- crop(landcover_type, lm_hexes, mask = TRUE)

# landcover_plot <- plot(habitat_landcover_rast)
#
# lc_pal <- landcover_plot$leg$fill
# names(lc_pal) <- landcover_plot$leg$legend

rast_vals <- unique(values(habitat_landcover_rast))

# pick an EMU
ug <- emus_proj %>% filter(UNIT == "Upper Gila Mountains")

# plot
rast_plot <- ggplot() +
  geom_spatvector(data = ug, fill = "white") +
  geom_spatraster(data = habitat_landcover_rast %>% crop(ug)) +
  #scale_fill_manual(breaks = rast_vals) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(#legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6))

## this didn't work
# legend <- cowplot::get_legend(rast_plot)
#
# grid.newpage()
# grid.draw(legend)
#
# cowplot::save_plot(here::here("figures/landcover_legend.jpeg"), landcover_legend)

