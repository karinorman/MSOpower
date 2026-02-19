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
  mutate(simulation_type = factor(ifelse(simulation_type == "hierarchical", "Range-wide", simulation_type),
                                  levels = c('BRE','BRW','CP','SRM','UGM', 'Range-wide')),
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


ggsave(file=here::here("figures/cost_plot.jpg"), plot = my_plot, dpi=600, width=250, height=250, units='mm')

cost_scen8_plot <- thresholds.sample %>%
  filter(simulation_scenario == 8) %>%
  rowwise() %>%
  mutate(plotcost = meancost / 1000000) %>%
  ungroup() %>%
  mutate(n.obs.type = str_replace_all(n.obs.type,
                                      pattern = "(FS)", replacement = "GE")) %>%
  ggplot(aes(x=ndeploy,y=plotcost,ymax=costmax/1000000,ymin=costmin/1000000,shape=factor(n.obs.type),col=factor(n.obs.type)))+
  geom_pointrange(size=0.3)+
  geom_line()+
  # labs(title="Optimal study design with respect to costs differs by number of deployments
  #      and number and types of observers with simulation scenario 8")+
  facet_grid(factor(ARUS)~simulation_type)+
  scale_y_continuous(labels = function(y) format(y, scientific = FALSE))+
  scale_x_continuous(breaks = seq(1,3, by=1))+

  #guides(color='none')+
  scale_color_manual("", values=pal)+
  scale_shape_manual("", values = c(16,17, 15, 8, 3, 17))+
  labs(x="Number of deployments", y="Cost (millions of USD)",
       shape="Number and types of observers")+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill="transparent"))


ggsave(file=here::here("figures/cost_scen8.jpg"), plot = cost_scen8_plot, dpi=600,width=210,height=150,units='mm')
