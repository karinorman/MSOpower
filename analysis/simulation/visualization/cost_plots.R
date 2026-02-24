#author: Jamie Sanderlin
#date: 1/26/2026

#purpose: this script takes results from the MSO power analysis
#         and combines the information with costs to provide some
#         optimal design scenarios

###################
# load functions  #
###################
source(here::here("R/logistics.sims.load.packages.R"))

################
# upload data  #
################
thresholds.sample <- read.csv(here::here('data/thresholds.sample.costs.v2.csv')) %>%
  mutate(simulation_type = case_when(
    simulation_type == "BRE" ~ "Basin & Range - East",
    simulation_type == "BRW" ~ "Basin & Range - West",
    simulation_type == "CP" ~ "Colorado Plateau",
    simulation_type == "SRM" ~ "Southern Rocky Mountains",
    simulation_type == "UGM" ~ "Upper Gila Mountains",
    simulation_type == "hierarchical" ~ "Range-wide"
  ),
  simulation_type = factor(simulation_type, levels = c("Basin & Range - East","Basin & Range - West","Colorado Plateau",
                 "Southern Rocky Mountains","Upper Gila Mountains", 'Range-wide')),
  ARUS = factor(paste(ARUS, "ARUs per hexagonal unit")))


#############
# plotting  #
#############

pal <- paletteer::paletteer_d("ggsci::default_uchicago")

my_plot <- thresholds.sample %>%
  rowwise() %>%
  mutate(plotcost = meancost / 1000000) %>%
  mutate(n.obs.type = str_replace_all(n.obs.type,
                                pattern = "(FS)", replacement = "GE")) %>%
  ggplot(aes(x=factor(ndeploy),y=plotcost,#ymax=costmax,ymin=costmin,
                     color=factor(simulation_scenario),shape=factor(n.obs.type)))+
  geom_jitter(size=1.5)+
  # labs(title="Optimal study design with respect to costs differs by number of deployments
  #      and number and types of observers over all simulation scenarios")+
  facet_grid(factor(ARUS)~factor(simulation_type))+
  scale_y_continuous(labels = function(y) format(y, scientific = FALSE))+
  scale_color_manual(name="Simulation Scenario",labels=c( expression("1: " * psi ~ "0.43, " * phi ~ "0.6, p 0.4"),
                                                          expression("2: " * psi ~ "0.43, " * phi ~ "0.6, p 0.8"),
                                                          expression("3: " * psi ~ "0.43, " * phi ~ "0.8, p 0.4"),
                                                          expression("4: " * psi ~ "0.43, " * phi ~ "0.8, p 0.8"),
                                                          expression("5: " * psi ~ "0.60, " * phi ~ "0.6, p 0.4"),
                                                          expression("6: " * psi ~ "0.60, " * phi ~ "0.6, p 0.8"),
                                                          expression("7: " * psi ~ "0.60, " * phi ~ "0.8, p 0.4"),
                                                          expression("8: " * psi ~ "0.60, " * phi ~ "0.8, p 0.8")),values=pal)+
  labs(x="Number of deployments", y="Cost (millions of USD)",
       shape="Number and types of observers")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="transparent"))


ggsave(file=here::here("figures/cost_plot.jpg"), plot = my_plot, dpi=600, width=350, height=250, units='mm')

test_data <-  thresholds.sample %>%
  filter(simulation_scenario == 8) %>%
  rowwise() %>%
  mutate(plotcost = meancost / 1000000) %>%
  ungroup() %>%
  mutate(n.obs.type = str_replace_all(n.obs.type,
                                      pattern = "(FS)", replacement = "GE")) %>%
  group_by(simulation_type, n.obs.type, ARUS) %>%
  mutate(opt_ind = ifelse(meancost == min(meancost), "Optimal Design", "Suboptimal Design"))

cost_scen8_plot <- test_data %>%
  #filter(ARUS == "2 ARUs per hexagonal unit", simulation_type == "Range-wide") %>%
  mutate(n.obs.type = factor(n.obs.type), opt_ind = factor(opt_ind)) %>%
  ggplot(aes(x=ndeploy)) +
  #geom_pointrange(size=0.3)+
  geom_line(aes(y=plotcost, color = n.obs.type, group = paste(n.obs.type, ARUS, simulation_type)))+
  geom_point(aes(y=plotcost, color = n.obs.type, shape = opt_ind, size = opt_ind)) +
  guides(color = guide_legend(override.aes = list(shape = NA))) +
  facet_grid(factor(ARUS)~simulation_type) +
  scale_y_continuous(labels = function(y) format(y, scientific = FALSE))+
  scale_x_continuous(breaks = seq(1,3, by=1))+
  guides(color = guide_legend(override.aes = list(shape = NA))) +
  scale_color_manual("", values=pal)+
  scale_shape_manual("", values = c(8, 16))+
  scale_size_manual("", values = c(2, 1)) +
  labs(x="Number of deployments", y="Cost (millions of USD)",
       color ="Number and types of observers")+
  theme_bw()+
  theme(#legend.position = 'bottom',
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="transparent"))


ggsave(file=here::here("figures/cost_scen8.jpg"), plot = cost_scen8_plot, dpi=600,width=350,height=150,units='mm')
