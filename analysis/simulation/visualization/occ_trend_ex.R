library(dplyr)
library(ggplot2)

ggplot() +
  # scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  geom_segment(aes(x = 1, xend = 5, y = 1, yend = 0.25), color = "#28587B", size = 1) +
  theme_classic() +
  ylab("Occupancy, \u03A8") +
  xlab("Year")

ggsave(here::here("figures/ex_occ_trend.jpg"), height = 3, width = 4)
