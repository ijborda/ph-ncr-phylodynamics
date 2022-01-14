# Library
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)
library(bdskytools)
library(tibble)

# Read file
lf = readLogfile("ncr-bdsky-time-post.log", burnin=0.1)

# Smooth Re
Re_sky    <- getSkylineSubset(lf, "reproductiveNumber")
Re_hpd    <- getMatrixHPD(Re_sky)
plotSkyline(1:6, Re_hpd, type='step', ylab="R")
timegrid       <- seq(0,median(lf$origin_BDSKY_Serial),length.out=101)
Re_gridded     <- gridSkyline(Re_sky, lf$origin, timegrid)
Re_gridded_hpd <- getMatrixHPD(Re_gridded)
times     <- 2020.57377-timegrid
plotSkyline(times, Re_gridded_hpd, type='smooth', xlab="Time", ylab="R")

# Convert to dataframe
data <- t(Re_gridded_hpd)
rownames(data) <- NULL
data <- data.table::as.data.table(data)
data$times <- date_decimal(times)
colnames(data) <- c("lower", "median", "upper", "date")

# Plot Re
p <- ggplot(data = data) +
  geom_line(aes(x=as.Date(date, "%m/%d/%y"), y=median), color = "dark blue", size = 1) +
  theme_classic() +
  labs(x="Date", y="Effective Reproduction Number") +
  scale_x_date(date_breaks = "1 month", 
               labels=date_format("%b-%Y")) +
  geom_ribbon(data = data, aes(x=as.Date(date, "%m/%d/%y"), ymin=lower, ymax=upper), fill="dark blue", alpha=0.4) +
  theme(panel.background = element_rect(fill = 'whitesmoke'), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0,10)) +
  geom_hline(yintercept = 1.0, linetype = 1, color = "red") +
  scale_x_date(date_breaks = "1 month", 
               labels=date_format("%b-%Y"),
               limits = as.Date(c('2020/01/01','2020/07/28'))) +
  geom_vline(xintercept = as.Date(c("2020/3/17")), 
             linetype = 1, 
             color = "black") +
  geom_vline(xintercept = as.Date(c("2020/5/16")), 
             linetype = 1, 
             color = "black") +
  geom_vline(xintercept = as.Date(c("2020/6/1")), 
             linetype = 1, 
             color = "black") 

# Save Plot
ggsave(plot = p,
       filename = "ncr-bdsky-time.png",
       width = 7, height = 5, units = "in", dpi = 300)
