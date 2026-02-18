library(dplyr)
library(ggplot2)
library(fGarch)

x <- seq(-0.6, 0.1, length.out = 1000)
y <- dsnorm(x, mean = -.35, sd = .07, xi = 1.5)

values <- rsnorm(1000, mean = -.35, sd = .07, xi = 1.5)

ci95 <- bayestestR::ci(values, ci = 0.95, method = "ETI")
ci90 <- bayestestR::ci(values, ci = 0.90, method = "ETI")

ci_plt <- ggplot(data.frame(x = x, y = y), aes(x, y)) +
  geom_line() +
  theme_classic() +
  scale_y_continuous(limits = c(0,6), expand=c(0,0)) +
  xlim(c(-0.63, 0.05)) +
  xlab("Estimated Percent Change") +
  ylab("Density") +
  geom_vline(xintercept = ci95$CI_low, color = "#87A96B", linetype = "dashed") +
  geom_vline(xintercept = ci95$CI_high, color = "#87A96B", linetype = "dashed") +
  geom_vline(xintercept = -0.06, color = "#d0b740", linetype = "dashed") +
  geom_vline(xintercept = ci90$CI_low, color = "#d0b740", linetype = "dashed") +
  geom_segment(aes(x = ci95$CI_low, y = 0.5, xend = ci95$CI_high, yend = 0.5), color = "#87A96B") +
  geom_segment(aes(x = 0, y = 0.2, xend = -.6, yend = 0.2), color = "#28587B", arrow = arrow(length=unit(0.2,"cm"), ends="last", type = "closed")) +
  theme(#axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #axis.line = element_line(colour = 'black', size = 2)
        ) +
  annotate("text", label = "\u03B1 = 0.025", x = ci95$CI_high - 0.05, y = 5, color = "#87A96B") +
  annotate("text", label = "\u03B1 = 0.05", x = -0.11, y = 4.8, color = "#d0b740") +
  geom_segment(aes(x = ci90$CI_low, y = 0.35, xend = -0.06, yend = 0.35), color = "#d0b740")

ggsave(here::here("figures/power_check_conceptual.jpg"), ci_plt, width = 6, height = 4)
