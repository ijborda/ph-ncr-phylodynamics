# Library
library(beastio)
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)


# Read data file
coal <- read.csv("ncr-coalsky.csv")

# Plot Ne
p <- ggplot(data = coal) +
  geom_line(aes(x=as.Date(date, "%m/%d/%y"), y=median), color = "dark blue", size = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  labs(x="Date", y="Effective Population Size") +
  scale_x_date(date_breaks = "1 month", 
               labels=date_format("%b-%Y")) +
  geom_ribbon(data = coal, aes(x=as.Date(date, "%m/%d/%y"), ymin=lower, ymax=upper), fill="dark blue", alpha=0.4) +
  scale_y_log10(limits=c(0.01,10)) +
  geom_vline(xintercept = as.Date(c("2020/3/17")), 
             linetype = 1, 
             color = "black") +
  geom_vline(xintercept = as.Date(c("2020/5/16")), 
             linetype = 1, 
             color = "black") +
  geom_vline(xintercept = as.Date(c("2020/6/1")), 
             linetype = 1, 
             color = "black") 

#Save Plot
ggsave(plot = p,
       filename = "ncr-coalsky.png",
       width = 7, height = 5, units = "in", dpi = 300)


