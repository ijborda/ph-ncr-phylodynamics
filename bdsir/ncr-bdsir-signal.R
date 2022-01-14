# Llibrary
library(beastio)
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)

# Read log files
df.pos = readLog(filename = "ncr-bdsir-post.log", 
                 burnin = 0.1, 
                 maxsamples = -1, 
                 as.mcmc =FALSE, 
                 burninAsSamples = FALSE)

df.pri = readLog(filename = "ncr-bdsir-prior.log", 
                 burnin = 0.1, 
                 maxsamples = -1, 
                 as.mcmc =FALSE, 
                 burninAsSamples = FALSE)

# Read only needed columns
df.pos.sam <- df.pos[,c("samplingProportionEs")]
df.pri.sam <- df.pri[,c("samplingProportionEs")]
df.sam <- data.frame(df.pos.sam,df.pri.sam)
colnames(df.sam) <- c("posterior", "prior")

df.pos.uni <- df.pos[,c("becomeUninfectiousRateEs")]
df.pri.uni <- df.pri[,c("becomeUninfectiousRateEs")]
df.uni <- data.frame(df.pos.uni,df.pri.uni)
colnames(df.uni) <- c("posterior", "prior")

df.pos.ori <- df.pos[,c("originEs")]
df.pri.ori <- df.pri[,c("originEs")]
df.ori <- data.frame(df.pos.ori,df.pri.ori)
colnames(df.ori) <- c("posterior", "prior")

df.pos.rep <- df.pos[,c("reproductiveNumberEs")]
df.pri.rep <- df.pri[,c("reproductiveNumberEs")]
df.rep <- data.frame(df.pos.rep,df.pri.rep)
colnames(df.rep) <- c("posterior", "prior")

# Plot graphs
sam <- ggplot(data=df.sam) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "a. Sampling Proportion")

uni <- ggplot(data=df.uni) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "b. Become Uninfectious Rate")

ori <- ggplot(data=df.ori) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "c. Origin")

rep <- ggplot(data=df.rep) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "d. Reproductive Number")


all <- plot_grid(sam, uni, ori, rep, ncol=2,align="v")


ggsave(plot = all,
       filename = "ncr-bdsir-signal.png",
       width = 7, height = 5, units = "in", dpi = 300)