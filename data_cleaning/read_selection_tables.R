##########################################
## Read in and format selection tables ###
##########################################

library(dplyr)
library(Rraven)


# get an example table to work with

file_path <- here::here("data/selection_tables/gila/S01169_000/S01169_2022-05-11/S01169_20220511_180000.BirdNET.selection.table.txt")

process_cibola_seltable <- function(file_path){

  sel_table <- utils::read.delim(file_path, header = TRUE,
                                 sep = "\t", fill = TRUE, stringsAsFactors = FALSE,
                                 check.names = FALSE)

  if(dim(sel_table)[1] > 1){
    sel_table %>%
      mutate(path = stringr::str_replace(file_path, ".*selection_tables/", "")) %>%
    tidyr::separate(path, c("forest", "sd_id", "unit", "date"), "/") %>%
    janitor::clean_names() %>%
    mutate(location = stringr::str_split(sd_id, "_")[[1]][3],
           sd_id = stringr::str_replace(sd_id, "_[^_]+$", ""),
           date = stringr::str_replace(date, ".*_", ""))
  } else (return(data.frame()))

}

process_gila_seltable <- function(file_path){

  sel_table <- utils::read.delim(file_path, header = TRUE,
                                 sep = "\t", fill = TRUE, stringsAsFactors = FALSE,
                                 check.names = FALSE)

  if(dim(sel_table)[1] > 1){
    sel_table %>%
      mutate(path = stringr::str_replace(file_path, ".*selection_tables/", "")) %>%
      tidyr::separate(path, c("forest", "unit", "date"), "/") %>%
      janitor::clean_names() %>%
      mutate(date = stringr::str_replace(date, ".*_", ""))
  } else (return(data.frame()))

}

process_kaibab_seltable <- function(file_path){

  sel_table <- utils::read.delim(file_path, header = TRUE,
                                 sep = "\t", fill = TRUE, stringsAsFactors = FALSE,
                                 check.names = FALSE)

  if(dim(sel_table)[1] > 1){
    sel_table %>%
      mutate(path = stringr::str_replace(file_path, ".*selection_tables/", "")) %>%
      tidyr::separate(path, c("forest", NA, NA, "unit", "date"), "/") %>%
      janitor::clean_names() %>%
      mutate(date = stringr::str_replace(date, ".*_", ""))
  } else (return(data.frame()))

}

cibola_paths <- list.files(here::here("data/selection_tables/cibola"), pattern = ".txt", full.names = TRUE, recursive = TRUE)
cibola_occ <- purrr::map_dfr(cibola_paths, process_cibola_seltable)

kaibab_paths <- list.files(here::here("data/selection_tables/kaibab"), pattern = ".txt", full.names = TRUE, recursive = TRUE)
kaibab_occ <- purrr::map_dfr(kaibab_paths, process_kaibab_seltable)

kaibab_spow_occ <- kaibab_occ %>%
  filter(species_code == "spoowl") %>%
  #90% confidence threshold
  filter(confidence >  0.80)

# let's get some metadata
spow_occ %>%
  select(unit)

##################
###### GILA ######
##################

gila_paths <- list.files(here::here("data/selection_tables/gila"), pattern = ".txt", full.names = TRUE, recursive = TRUE)
gila_occ <- purrr::map_dfr(gila_paths, process_gila_seltable)

gila_spow_occ <- gila_occ %>%
  filter(species_code == "spoowl") %>%
  #90% confidence threshold
  filter(confidence >  0.9895) %>%
  mutate(date = lubridate::ymd(date))

# let's try to figure out the sampling scheme (which dates units were deployed)
deployment_info <- gila_spow_occ %>%
  select(unit, date) %>%
  distinct() %>%
  group_by(unit) %>%
  summarise(deploy_start = min(date), deploy_end = max(date)) %>%
  mutate(date_range = deploy_end - deploy_start)

# get range between observations for a unit
gila_spow_occ %>%
  group_by(forest, unit) %>%
  summarize(date_range = max(date) - min(date))
