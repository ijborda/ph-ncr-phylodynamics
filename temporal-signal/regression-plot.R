# Library
library(ggtree)
library(treeio)
library(ggplot2)

# TempEst Output
slope = 6.1877E-4
tmrca = 2019.8348
rsqr	= 0.415
genDist <- read.csv("genetic-distance.csv")

# Plot Limits
xmin = 2019.8
xmax = 2020.7
ymin = 0
ymax = 0.0006

# Regression
serror = 4.85075E-05;

#Plot  regression line
p1 <- ggplot(genDist, aes(x= date, y = distance)) + 
  geom_point(size = 5, color='blue', alpha = 0.5) + 
  stat_function(fun = function(x){slope*x + -slope*tmrca}, color = "red", size = 1) + 
  stat_function(fun = function(x){slope*x + -slope*tmrca + 2*serror}, color = "red", size = 1, linetype = "dashed") +
  stat_function(fun = function(x){slope*x + -slope*tmrca - 2*serror}, color = "red", size = 1, linetype = "dashed") +
  labs(y = "Genetic Distance (Root-to-tip divergence)", 
       x = "Sampling Time (in years)",
       title = "Temporal Signal")  + 
  scale_x_continuous(limits = c(xmin, xmax)) + 
  scale_y_continuous(limits = c(ymin, ymax), labels = scales::number_format(accuracy = 0.0001)) +
  annotate(geom = "text", label = ~~italic(R)^2~" = 0.415;"~~italic(p-val)~"= 0.005", x = 2020, y = 0.00055, size = 4) +
  theme_bw() 

ggsave(plot = p1,
       filename = "temporal-signal.png",
       width = 5, height = 4, units = "in", dpi = 300)

