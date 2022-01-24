# Library
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)

# Read data file
r <- read.csv("r-best-10.csv")

# Plot R0
p <- ggplot(r, aes(y = r)) + 
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
  labs(y=expression("R"[0]))
 
#Save Plot
ggsave(plot = p,
       filename = "odesir-r.png",
       width = 5, height = 7, units = "in", dpi = 300)