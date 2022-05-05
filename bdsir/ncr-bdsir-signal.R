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

# Transform origin datafram to dates
latest.sample = ymd("2020-07-18")
df.ori$posterior = latest.sample - df.ori$posterior*365 
df.ori$prior = latest.sample - df.ori$prior*365

# Plot graphs
ori <- ggplot(data=df.ori) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "Date of Origin") +
  scale_x_date(date_breaks = "2 month", labels=date_format("%b-%Y")) + 
  geom_vline(data=df.ori, aes(xintercept=mean(posterior)),  colour="dark blue") +
  geom_vline(data=df.ori, aes(xintercept=quantile(posterior, 0.05, type = 1)), colour="dark blue", linetype="dashed") +
  geom_vline(data=df.ori, aes(xintercept=quantile(posterior, 0.95, type = 1)), colour="dark blue", linetype="dashed")

rep <- ggplot(data=df.rep) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "Basic Reproductive Number") +
  geom_vline(data=df.rep, aes(xintercept=mean(posterior)),  colour="dark blue") +
  geom_vline(data=df.rep, aes(xintercept=quantile(posterior, 0.05, type = 1)), colour="dark blue", linetype="dashed") +
  geom_vline(data=df.rep, aes(xintercept=quantile(posterior, 0.95, type = 1)), colour="dark blue", linetype="dashed")

uni <- ggplot(data=df.uni) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "Become Uninfectious Rate") +
  geom_vline(data=df.uni, aes(xintercept=mean(posterior)),  colour="dark blue") +
  geom_vline(data=df.uni, aes(xintercept=quantile(posterior, 0.05, type = 1)), colour="dark blue", linetype="dashed") +
  geom_vline(data=df.uni, aes(xintercept=quantile(posterior, 0.95, type = 1)), colour="dark blue", linetype="dashed")

sam <- ggplot(data=df.sam) +
  geom_density(aes(x=posterior,y=..density..), color=FALSE , fill="blue", alpha=.3) +
  geom_density(aes(x=prior,y=..density..), size=1, color='black', linetype="dashed") +
  theme_classic() +
  labs(x = " ", y = "Marginal Density", title = "Sampling Proportion") +
  geom_vline(data=df.sam, aes(xintercept=mean(posterior)),  colour="dark blue") +
  geom_vline(data=df.sam, aes(xintercept=quantile(posterior, 0.05, type = 1)), colour="dark blue", linetype="dashed") +
  geom_vline(data=df.sam, aes(xintercept=quantile(posterior, 0.95, type = 1)), colour="dark blue", linetype="dashed")

all <- plot_grid(ori, rep, uni, sam, ncol=2,align="v",
                 labels = c('A', 'B', 'C', 'D'))

# Kolmogorov-Smirnov test
pval_ori <- ks.test(as.numeric(df.ori$posterior), as.numeric(df.ori$prior)) #p-value<2.2e-16
pval_rep <- ks.test(as.numeric(df.rep$posterior), as.numeric(df.rep$prior)) #p-value<2.2e-16
pval_uni <- ks.test(as.numeric(df.uni$posterior), as.numeric(df.uni$prior)) #p-value<2.2e-16
pval_sam <- ks.test(as.numeric(df.sam$posterior), as.numeric(df.sam$prior)) #p-value<2.2e-16

ggsave(plot = all,
       filename = "ncr-bdsir-signal.png",
       width = 7, height = 5, units = "in", dpi = 300)