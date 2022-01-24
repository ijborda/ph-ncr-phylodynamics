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
    time=seq(0,139,by=1),
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
  params=c(X.0=6999970.8, Y.0=29.2, Z.0=0, N.0=7000000, cases.0=0,
           mu=0, Beta=0.0891, gamma=0.0611, rho=0.5)
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
models <- read.delim(file = "model-sims.csv" , sep = ',', header = TRUE)
sde$date <- models$date;

# Plot models
p <- ggplot() +
  geom_col(models, mapping = aes(x=as.Date(date, format="%m/%d/%y"), y=reported), color = "khaki",  fill = "khaki", width = 1, alpha = 0.5) +
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
  geom_line(models, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=ode, colour = "Deterministic SIR"), size = 1.0) +
  geom_line(models, mapping = aes(x=as.Date(bdsir_date,format="%m/%d/%y"), y=bdsir, colour = "BDSIR"), size = 1.0) +
  geom_line(sde, mapping = aes(x=as.Date(date,format="%m/%d/%y"), y=median, colour = "Stochastic SIR (median)"), size = 1.0) +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1),
        text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  labs(x="Date", y="Number of Infected Individuals") +
  scale_x_date(date_breaks = "1 month", 
               labels=date_format("%b-%Y"),
               limits = as.Date(c("3/1/20", "7/18/2020"), format="%m/%d/%y")) +
  scale_colour_manual("", 
                      breaks = c("Deterministic SIR", "Stochastic SIR (median)", "BDSIR", "Stochastic SIR"),
                      values = c("dark blue", "dark green", "dark red", "gray70")) +
  theme(legend.position = "top")

# Save plot
ggsave(plot = p,
       filename = "models.png",
       width = 7, height = 5, units = "in", dpi = 300)
