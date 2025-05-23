library(dplyr)
library(terra)
library(tidyterra)
library(ggplot2)
library(sf)
library(cowplot)


emus <- st_read(here::here("data/MSO_EMUs/MSO_EMUs.shp")) %>%
  select(UNIT) %>%
  # equal area projection
  st_transform(crs = 4326)

# if we have to read in again, names are messed up
grid_attr_habitat <- vect(here::here("data/grid_attr_habitat.shp"))
names(grid_attr_habitat)[23:29] <- c("sample_frame_veg", "sample_frame_type",  "habitat_2000", "habitat_2022", "mso_habitat_type", "mso_percent_habitat")

grid_sample_frame <- vect(here::here("data/grid_sample_frame.shp"))
names(grid_sample_frame) <- c("ID", "UNIT", "veg_type_landfire", "habitat_2000", "habitat_2022", "mso_habitat_type", "mso_percent_habitat", "include_patch")

grid_veg_samp <- grid_sample_frame %>%
  filter(include_patch == "yes") %>%
  select(veg_type_landfire) %>%
  mutate(plot_lc = case_when(
    veg_type_landfire == "Madrean Pinyon-Juniper Woodland" ~ "Pinyon-Juniper",
    veg_type_landfire == "Madrean Lower Montane Pine-Oak Forest and Woodland" ~ "Pine/Oak",
    veg_type_landfire == "Southern Rocky Mountain Dry-Mesic Montane Mixed Conifer Forest and Woodland" ~ "Mixed Conifer",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Woodland" ~ "Pine",
    veg_type_landfire == "Madrean Upper Montane Conifer-Oak Forest and Woodland" ~ "Pine/Oak",
    veg_type_landfire == "Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland" ~ "Subalpine",
    veg_type_landfire == "Rocky Mountain Aspen Forest and Woodland" ~ "Mixed Conifer",
    veg_type_landfire == "Inter-Mountain Basins Aspen-Mixed Conifer Forest and Woodland" ~ "Mixed Conifer",
    veg_type_landfire == "Southern Rocky Mountain Ponderosa Pine Savanna" ~ "Pine",
    veg_type_landfire == "Inter-Mountain Basins Subalpine Limber-Bristlecone Pine Woodland" ~ "Subalpine",
    veg_type_landfire == "Rocky Mountain Bigtooth Maple Ravine Woodland" ~ "Mixed Conifer",
    veg_type_landfire == "Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland" ~ "Subalpine",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Riparian Woodland" ~ "Subalpine",
    veg_type_landfire == "Rocky Mountain Lodgepole Pine Forest" ~ "Subalpine",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "Subalpine",
    veg_type_landfire == "Rocky Mountain Lodgepole Pine Forest" ~ "Pine",
    veg_type_landfire == "Rocky Mountain Subalpine-Montane Limber-Bristlecone Pine Woodland" ~ "Subalpine",
    veg_type_landfire == "Southern Rocky Mountain Mesic Montane Mixed Conifer Forest and Woodland" ~ "Mixed Conifer",
    .default = NA
  ))

ext_box <- ext(emus)
ext_box[4] <- 41

# landcover <- rast(here::here("data/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif"))

# lc_plot <- plot(landcover)
# pal <- lc_plot$leg$fill
# names(pal) <- lc_plot$leg$legend

# n_lc_types <- n_distinct(grid_sample_frame$veg_type_landfire)
# pal <- colorspace::terrain_hcl(n_lc_types + 1)

pal <- c("#0B6623", "#87A96B", "#D8E4BC", "#CABA9C", "#8A6240")
names(pal) <- c("Mixed Conifer", "Pine", "Pine/Oak", "Pinyon-Juniper", "Subalpine")

boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Arizona", "Colorado", "New Mexico", "Utah", "Texas")) %>%
  crop(ext_box)

crop_box <- st_bbox(c(xmin = -105.7, xmax = -104.5, ymax = 39.3, ymin = 38.5), crs = st_crs(4326))

base_map <- ggplot() +
  geom_spatvector(data = boundary_states, fill = "transparent") +
  geom_spatvector(data = emus, color = "white", fill = "lightgrey", alpha = 0.5, linewidth = .5) +
  geom_spatvector(data = boundary_states %>% aggregate(), fill = "transparent" ) +
  geom_spatvector(data = grid_veg_samp, aes(fill = plot_lc, color = plot_lc)) +
  geom_spatvector(data = vect(st_as_sfc(crop_box)), color = "black", fill = "transparent", linewidth = 0.6) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  theme_void() +
  theme(legend.position = "inside",
        legend.position.inside = c(1.15, .35),
        legend.title = element_blank(),
        plot.margin = margin(2, 8, .5, 0.5, "cm")) +
  annotate("label", x = -112.5, y = 42, label = "Colorado \nPlateau") +
  annotate("label", x = -106, y = 42, label = "Southern Rocky \nMountains") +
  annotate("label", x = -112.5, y = 30.5, label = "Basin and \nRange West") +
  annotate("label", x = -108, y = 30.5, label = "Upper Gila \nMountains") +
  annotate("label", x = -108, y = 30.5, label = "Upper Gila \nMountains") +
  annotate("label", x = -103, y = 35, label = "Basin and \nRange East")

ggsave(here::here("figures/emu_basemap.jpg"), base_map)

inset_map <- ggplot() +
  geom_spatvector(data = emus %>% st_crop(crop_box), color = "white", fill = "lightgrey", alpha = 0.5, linewidth = .8) +
  geom_spatvector(data = grid_veg_samp %>% crop(crop_box), aes(fill = plot_lc), color = "white") +
  scale_fill_manual(values = pal) +
  theme_void()+
  theme(legend.position = "none") +
  geom_spatvector(data = vect(st_as_sfc(crop_box)), color = "black", fill = "transparent", linewidth = 0.5)

study_map <- ggdraw(base_map) +
  draw_plot(
    {inset_map},
    x = .59,
    y = .55,
    width = 0.4,
    height = 0.4
  )

ggsave(here::here("figures/sample_frame_map.jpeg"), study_map, width = 210, height = 180, units = "mm")

