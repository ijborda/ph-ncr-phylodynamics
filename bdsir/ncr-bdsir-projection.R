# Library
library(pomp)
library(ggplot2)
library(reshape2)
library(dplyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)

# Set the model
pomp(
  data=data.frame(
    time=seq(0,737,by=1),
    reports=NA
  ),
  times="time",t0=0,
  rprocess=gillespie_hl(
    birth=list("rate=mu*N;",c(N=1,X=1,Y=0,Z=0,cases=0)),
    deathS=list("rate=mu*X;",c(N=-1,X=-1,Y=0,Z=0,cases=0)),
    deathI=list("rate=mu*Y;",c(N=-1,X=0,Y=-1,Z=0,cases=0)),
    deathR=list("rate=mu*Z;",c(N=-1,X=0,Y=0,Z=-1,cases=0)),
    infection=list("rate=Beta*X*Y/N;",c(N=0,X=-1,Y=1,Z=0,cases=1)),
    recovery=list("rate=gamma*Y;",c(N=0,X=0,Y=-1,Z=1,cases=0)),
    hmax=0.01
  ),
  rmeasure=Csnippet("reports=rbinom(cases,rho);"),
  paramnames=c("rho","mu","Beta","gamma"),
  statenames=c("N","X","Y","Z","cases"),
  accumvars=c("cases"),
  params=c(X.0=6999988.2, Y.0=11.8, Z.0=0, N.0=7000000, cases.0=0,
           mu=0, Beta=0.1062, gamma=0.0284, rho=0.5)
) -> sir

# Produce simulations
simulate(sir, nsim=10, format="data.frame") -> sims
sde <- subset(sims, select = c(time, .id, Y))
sde <- reshape(sde, idvar = "time", timevar = ".id", direction = "wide")

# Calculate median of simulations
sde <- sde %>% 
  rowwise() %>% 
  mutate(median = median(c(Y.1, Y.2, Y.3, Y.4, Y.5, Y.6, Y.7, Y.8, Y.9, Y.10), na.rm = TRUE))

# Load simulation results from other models
models <- read.delim(file = "ncr-bdsir-ode-projection.csv" , sep = ',', header = TRUE)
sde$date <- models$date;

# Load reported data
reported <- read.delim(file = "ncr-reported-projection-data.csv" , sep = ',', header = TRUE)

# Define scale multiplier
mult <- 400;

# Plot models
p <- ggplot() +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.1, colour = "Stochastic SIR",), size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.2), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.3), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.4), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.5), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.6), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.7), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.8), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.9), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=Y.10), color = "gray70", size = 1.0, alpha = 1) +
  geom_line(models, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=genomic, colour = "Deterministic SIR"), size = 1.0) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=median, colour = "Stochastic SIR (median)"), size = 1.0) +
  geom_line(models, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=reported, colour = "Deterministic SIR (using reported data)"), size = 1.0) +
  geom_col(reported, mapping = aes(x=as.Date(date, format="%d-%b-%y") , y=reported*mult, fill = "Reported data")) +
  theme_classic() +
  labs(x="Date", y="Number of Projected Cases") +
  scale_x_date(date_breaks = "3 month", 
               labels=date_format("%b-%Y")) +
  scale_colour_manual(" ",
                      values = c("Deterministic SIR" = "dark blue", 
                                 "Stochastic SIR (median)" = "dark green", 
                                 "Stochastic SIR" = "gray70", 
                                 "Deterministic SIR (using reported data)" = "dark red"))  +
  scale_fill_manual("",values=alpha("orange",0.5)) + 
  theme(legend.position = "top") +
  geom_vline(xintercept = as.Date("9/5/20", format="%m/%d/%y"), linetype = 4, color="black") +
  theme(axis.text.x = element_text(hjust = 1),
        text = element_text(size=15),
        axis.title.y.left = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        legend.key=element_blank()) +
  scale_y_continuous(label = comma,
                     sec.axis=sec_axis(~.*1/mult,name="Number of Reported Cases")) +
  annotate("rect", xmin = as.Date("7/18/20", format="%m/%d/%y"), xmax = as.Date("3/8/22", format="%m/%d/%y"), 
           ymin = 0, ymax = 2800000,
           alpha = .2) +
  annotate("text", x = as.Date("1/8/22", format="%m/%d/%y"), y = 2600000, label = "Projection")

# Save plot
ggsave(plot = p,
       filename = "ncr-bdsir-projection.png",
       width = 13, height = 5, units = "in", dpi = 300)
