#### Explore survival (phi)
#### Are the simulated annual survival values linearly declining

library(dplyr)

files <- data.frame(files = list.files(here::here("data/nimble/emu_simulated_data/")),
                    file_paths = list.files(here::here("data/nimble/emu_simulated_data/"), full.names = TRUE))

check_phi_df <- files %>%
  separate(files, c("EMU", "sim_id", "rep", "scenario_id")) %>%
  unite("high_name", c("EMU", "sim_id")) %>%
  mutate(scenario_id = as.integer(scenario_id)) %>%
  left_join(sim_map_names %>% select(high_name, low_name, scenario_id)) %>%
  filter(!is.na(low_name)) %>%
  separate(high_name, c("emu", "sim_num"), sep = "_", remove = FALSE) %>%
  filter(emu == "SRM") %>%
  left_join(sim_map_names) #%>%
  #filter(psi == 0.43, phi == 0.8)


annual_phi <- pmap_dfr(check_phi_df %>% select(high_name, rep, file_paths), function(high_name, rep, file_paths){
  #browser()
  data.frame(survival = readRDS(here::here(file_paths))$phi_survival,
             year = 1:9) %>%
    mutate(high_name, rep)
})


ggplot(annual_phi) +
  geom_point(aes(x = year, y = survival)) +
  geom_smooth(aes(x = year, y = survival), method = "loess")

ggplot(annual_phi) +
  geom_point(aes(x = year, y = car::logit(survival))) +
  geom_smooth(aes(x = year, y = car::logit(survival)), method = "loess")
