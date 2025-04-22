library(dplyr)
library(ggplot2)

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
  crop(boundary_states)

pal <- c(
  "Cibola National Forest"   = "#3F6634",
  "Kaibab National Forest"    = "#586994",
  "Gila National Forest" = "#C8963E"
)

# map of the entire study area (New Mexico and Arizona)
study_map <- ggplot() +
  geom_spatvector(data = boundary_states, color = "black", fill = "transparent") +
  geom_spatvector(data = forest_boundaries, aes(fill = FORESTNAME), color = "transparent", alpha = 0.60) + #, color = "grey") +
  geom_spatvector(data = all_locs_vect %>% crop(crop_bb), color = "black", size = 0.005) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  ggthemes::theme_map() +
  theme(legend.title = element_blank()) +
  geom_rect(
    xmin = -112.2353,
    ymin = 35.18814,
    xmax = -112.1877,
    ymax = 35.22706,
    fill = NA,
    colour = "black",
    linewidth = 1
  ) + theme(plot.margin = margin(5,7,5,1, "cm"),
            legend.position.inside = c(.2, .17))

# let's zoom in on kaibab points

get_inset_map <- function(forest_name, buffer_dist){
  bb <- all_locs_vect %>%
    crop(crop_bb) %>%
    filter(forest == forest_name) %>%
    st_as_sf() %>%
    terrainr::add_bbox_buffer(distance = buffer_dist, distance_unit = "km") %>%
    vect()

  shape_color <- ifelse(forest_name == "kaibab", "transparent", "black")
  filter_forest <- paste(tools::toTitleCase(forest_name), "National Forest", sep = " ")

  ggplot() +
    geom_spatvector(data = forest_boundaries %>% crop(bb),
                    fill = "transparent", color = shape_color) +
                    #aes(fill = FORESTNAME, color = FORESTNAME), alpha = 0.75)  +
    geom_spatvector(data =all_locs_vect %>% crop(bb), size = 0.5, color = pal[filter_forest]) +
    #scale_color_manual(values = pal) +
    #scale_fill_manual(values = pal) +
    ggthemes::theme_map() +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = pal[filter_forest], size=2))
}

get_inset_map("kaibab", 6)

plt_df <- data.frame(forest_name = c("kaibab", "cibola", "gila"), buffer_dist = c(4, 15, 12))

forest_insets <- pmap(plt_df, get_inset_map) %>% set_names(c("kaibab", "cibola", "gila"))

library(cowplot)

ggdraw(study_map) +
  draw_plot(
    {forest_insets$kaibab},
    x = 0.2,
    y = 0.7,
    width = 0.3,
    height = 0.3
  ) +
  draw_plot(
    {forest_insets$cibola},
    x = 0.75,
    y = 0.3,
    width = 0.2,
    height = 0.5
  ) +
  draw_plot(
    {forest_insets$gila},
    x = 0.4,
    y = 0.01,
    width = 0.3,
    height = 0.3
  )

# ggdraw(study_map) +
#   draw_plot(
#     {
#       study_map +
#         coord_sf(
#           xlim = c(-112.2353, -112.1877),
#           ylim = c(35.18814, 35.22706),
#           expand = FALSE) +
#         theme(legend.position = "none")
#     },
#     x = 0.05,
#     y = 0.5,
#     width = .8,
#     height = .8
#   ) +
#   draw_plot(
#     {
#       study_map +
#         coord_sf(
#           xlim = c(-106.7749, -106.0331),
#           ylim = c(34.59445, 35.20282),
#           expand = FALSE) +
#         theme(legend.position = "none")
#     },
#     x = 0.5,
#     y = 0.5,
#     width = .8,
#     height = .8
#   )

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
