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

bb <- st_bbox(all_locs_vect) %>% as.list()
crop_bb <- st_bbox(c(xmin = -113, xmax = bb$xmax, ymax = bb$ymax, ymin = bb$ymin), crs = st_crs(4326))

boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Arizona","New Mexico"))

forest_boundaries <- vect(here::here("data/S_USA.AdministrativeForest/S_USA.AdministrativeForest.shp")) %>%
  filter(FORESTNAME %in% c("Cibola National Forest", "Gila National Forest", "Kaibab National Forest")) %>%
  crop(boundary_states)


ggplot() +
  geom_spatvector(data = boundary_states, color = "black", fill = "transparent") +
  geom_spatvector(data = forest_boundaries, aes(fill = FORESTNAME), color = "grey") +
  geom_spatvector(data = all_locs_vect %>% crop(crop_bb), color = "black", size = 0.005) +
  ggthemes::theme_map()

ggplot() +
  geom_spatvector(data = forest_boundaries %>% filter(FORESTNAME ==  "Kaibab National Forest")) +
  geom_spatvector(data =all_locs_vect %>% crop(crop_bb) %>% filter(forest == "kaibab"), size = 0.005) +
  ggthemes::theme_map()
