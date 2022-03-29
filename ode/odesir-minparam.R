# Library
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)

# Read data file
r <- read.csv("odesir-minparam.csv")

# Plot R0
beta <- ggplot(r, aes(y = Beta)) + 
  geom_boxplot(coef=3, fill="dark blue", alpha = 0.4) +
  scale_y_continuous(breaks = seq(0.06, 0.12, by = .01),
                     limits=c(0.06, 0.12)) +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=15)) +
  labs(y=expression(paste("Effective Transmission Rate (", beta, ")")))

gamma <- ggplot(r, aes(y = Gamma)) + 
  geom_boxplot(coef=3, fill="dark blue", alpha = 0.4) +
  scale_y_continuous(breaks = seq(0.03, 0.09, by = .01),
                     limits=c(0.03, 0.09)) +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=15)) +
  labs(y=expression(paste("Removal Rate (", gamma, ")")))

rnaught <- ggplot(r, aes(y = R0)) + 
  geom_boxplot(coef=3, fill="dark blue", alpha = 0.4) +
  scale_y_continuous(breaks = seq(1.3, 1.7, by = .05),
                     limits=c(1.3, 1.7)) +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=15)) +
  labs(y=expression(paste("Basic Reproductive Number (", R[0], ")")))
 
title <- ggdraw() + 
  draw_label("Epidemiological Parameter Estimates from the SIR model", x = 0, hjust = 0) +
  theme( plot.margin = margin(0, 0, 0, 7))
all <- plot_grid(beta, gamma, rnaught, ncol=3, align="v",
                 labels = c('A', 'B', 'C'),
                 scale = 0.8) 
allTitle <- plot_grid(title, all, ncol = 1, rel_heights = c(0.1, 1))

#Save Plot
ggsave(plot = allTitle,
       filename = "odesir-minparam.png",
       width = 7, height = 5, units = "in", dpi = 300)