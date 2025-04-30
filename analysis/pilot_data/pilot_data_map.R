library(dplyr)
library(ggplot2)
library(cowplot)
library(terra)
library(tidyterra)
library(sf)
library(ggspatial)
library(purrr)

#### Plot of all ARU locations ####

kaibab_meta <- read.csv(here::here("data/bioacoustics.2022.Kaibab.metadata.csv"))
gila_meta <- read.csv(here::here("data/bioacoustics.2022.Gila.metadata.csv"))
cibola_locs <- read.csv("data/allspp_obs.csv") %>%
  filter(forest == "cibola")


all_locs <- cibola_locs %>%
  select(location) %>%
  distinct() %>%
  mutate(location = ifelse(location == "N,3832062,E369715", "N3832062,E369715", location)) %>%
  tidyr::separate_wider_delim(location, delim = ",", names = c("utmN", "utmE")) %>%
  mutate(across(c(utmN, utmE), ~as.integer(stringr::str_remove(.x, "[a-zA-Z]+")))) %>%
  mutate(forest = "cibola", utmZone = 13, datum = "WGS84") %>%
  bind_rows(kaibab_meta %>%
              select(utmN, utmE, datum, utmZone) %>%
              distinct() %>%
              mutate(forest = "kaibab"),
            gila_meta %>%
              select(utmN, utmE, datum = Datum, utmZone) %>%
              distinct() %>%
              mutate(forest = "gila")) %>%
  mutate(epsg = case_when(
    datum == "NAD83" ~ 26912,
    datum == "WGS84" & utmZone == 13 ~ 32613,
    datum == "WGS84" & utmZone == 12 ~ 32612,
    .default = NA
  ))

all_locs_vect <- all_locs %>%
  group_by(epsg) %>%
  group_map(~as_spatvector(x = .x, geom = c("utmE", "utmN"), crs = paste0("EPSG:", .y)) %>% project("epsg:4326")) %>%
  unlist() %>%
  vect()

writeVector(all_locs_vect, filename = here::here("data/pilot_aru_locs.shp"))

bb <- st_bbox(all_locs_vect) %>% as.list()
crop_bb <- st_bbox(c(xmin = -113, xmax = bb$xmax, ymax = bb$ymax, ymin = bb$ymin), crs = st_crs(4326))

boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Arizona","New Mexico"))

forest_boundaries <- vect(here::here("data/S_USA.AdministrativeForest/S_USA.AdministrativeForest.shp")) %>%
  filter(FORESTNAME %in% c("Cibola National Forest", "Gila National Forest", "Kaibab National Forest")) %>%
  project("epsg:4326") %>%
  crop(boundary_states)

landcover <- rast(here::here("data/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif")) %>%
  crop(boundary_states %>% project("epsg:5070")) %>%
  project("epsg:4326")

pal <- c(
  "Cibola National Forest"   = "#3F6634",
  "Kaibab National Forest"    = "#586994",
  "Gila National Forest" = "#C8963E"
)

# let's get the boxes that represent inset zooms
bb_list <- pmap(data.frame(forest_name = c("kaibab", "cibola", "gila"),
                           buffer_dist = c(10, 18, 12)),
                           function(forest_name, buffer_dist) {all_locs_vect %>%
                 crop(crop_bb) %>%
                 filter(forest == forest_name) %>%
                 st_as_sf() %>%
                 terrainr::add_bbox_buffer(distance = buffer_dist, distance_unit = "km") %>%
                 st_bbox()}
                 ) %>% set_names(c("kaibab", "cibola", "gila"))

# map of the entire study area (New Mexico and Arizona)
study_map <- ggplot() +
  geom_spatvector(data = boundary_states, color = "black", fill = "transparent") +
  geom_spatvector(data = forest_boundaries, aes(fill = FORESTNAME), color = "transparent", alpha = 0.60) + #, color = "grey") +
  geom_spatvector(data = all_locs_vect %>% crop(crop_bb), color = "black", size = 0.005) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  ggthemes::theme_map() +
  theme(legend.title = element_blank()) +
  theme(plot.margin = margin(7,7,7,1, "cm"),
            legend.position.inside = c(.14, .17),
            legend.text = element_text(size = 5),
            theme(legend.key.size = unit(1,"line")),
            legend.background = element_rect(fill = "transparent")) +
  geom_rect(aes(xmin = bb_list$kaibab$xmin, xmax = bb_list$kaibab$xmax, ymin = bb_list$kaibab$ymin, ymax = bb_list$kaibab$ymax),
            fill = "transparent", color = "black",
            linewidth = 0.4) +
  geom_rect(aes(xmin = bb_list$gila$xmin, xmax = bb_list$gila$xmax, ymin = bb_list$gila$ymin, ymax = bb_list$gila$ymax),
            fill = "transparent", color = "black",
            linewidth = 0.4) +
  geom_rect(aes(xmin = bb_list$cibola$xmin, xmax = bb_list$cibola$xmax, ymin = bb_list$cibola$ymin, ymax = bb_list$cibola$ymax),
            fill = "transparent", color = "black",
            linewidth = 0.4)

# let's zoom in on kaibab points

get_inset_map <- function(forest_name, buffer_dist, pt_size){
  bb <- all_locs_vect %>%
    crop(crop_bb) %>%
    filter(forest == forest_name) %>%
    st_as_sf() %>%
    terrainr::add_bbox_buffer(distance = buffer_dist, distance_unit = "km") %>%
    vect()

  shape_color <- ifelse(forest_name == "kaibab", "transparent", "black")
  filter_forest <- paste(tools::toTitleCase(forest_name), "National Forest", sep = " ")

  ggplot() +
    geom_spatraster(data = landcover %>% crop(bb), alpha = 0.5) +
    geom_spatvector(data = forest_boundaries %>% crop(bb), linewidth = 0.4,
                    fill = "transparent", color = shape_color) +
    geom_spatvector(data =all_locs_vect %>% crop(bb), size = pt_size, color = "black")+
    geom_spatvector(data = bb, color = pal[filter_forest], fill = "transparent", linewidth = 1) +
    #scale_color_manual(values = pal) +
    #scale_fill_manual(values = pal) +
    ggthemes::theme_map() +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill = "yellow"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          plot.margin = margin(0,0,0,0)
          ) +
    annotation_scale(pad_x = unit(0.5, "cm"),
                    pad_y = unit(0.5, "cm"),
                    height = unit(0.1, "cm"),
                    text_cex = 0.5)
}

get_inset_map("kaibab", 2, 2)

plt_df <- data.frame(forest_name = c("kaibab", "cibola", "gila"), buffer_dist = c(4, 15, 12), pt_size = c(1,.15, .5))

forest_insets <- pmap(plt_df, get_inset_map) %>% set_names(c("kaibab", "cibola", "gila"))

study_map_insets <- ggdraw(study_map) +
  draw_plot(
    {forest_insets$kaibab},
    x = 0.1,
    y = 0.64,
    width = 0.35,
    height = 0.35
  ) +
  draw_plot(
    {forest_insets$cibola},
    x = 0.64,
    y = 0.12,
    width = 0.35,
    height = 0.85
  ) +
  draw_plot(
    {forest_insets$gila},
    x = 0.35,
    y = 0.01,
    width = 0.35,
    height = 0.35
  )

ggsave(here::here("figures/pilot_locations.jpeg"), study_map_insets, width = 180, height = 200, units = "mm")


# Map of just the Gila w/PACs
#
# pacs <- vect(here::here("data/MSO_PACs/MSO_PACs.shp")) %>%
#   project("epsg:4326")
#
# ggplot() +
#   geom_spatvector(data = forest_boundaries %>% filter(FORESTNAME ==  "Gila National Forest")) +
#   geom_spatvector(data =all_locs_vect %>% crop(crop_bb) %>% filter(forest == "gila"), size = 0.005, color = "red") +
#   geom_spatvector(data = pacs %>% crop(forest_boundaries %>% filter(FORESTNAME == "Gila National Forest"))) +
#   ggthemes::theme_map()
