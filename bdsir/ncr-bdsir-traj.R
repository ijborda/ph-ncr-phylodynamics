# Library
library(beastio)
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(gridExtra)
library(cowplot)
library(s20x)
library(boa)
library(ggforce)
library(data.table)

# Compute time and SIR
colnameindex = function(M , colname0) 
{ 
colsnames = names(M[1,]); 
theindex = which(colsnames==colname0); 
return(theindex); 
}
loglist = read.table("loglist.txt", as.is=TRUE, header=F) #list of log filename(s)
finalSampleTime = 2020.546448 # time of most recent sample (18 July 2020)
dim =length(loglist[,1])
intervals=101
for(i in 1:dim){
  
  log = read.table(loglist[i,], header=TRUE)
  samples=length(log[,1])
  burnin = round(.1 * samples)
  ds0 = rep(NA,samples-burnin+1)
  dr0 = rep(NA,samples-burnin+1)
  
  R_perInterval=matrix(data = NA, nrow = dim, ncol = intervals)
  R_perInterval_HPD=matrix(data = NA, nrow = 2, ncol = intervals)
  allS = matrix(NA, ncol=intervals, nrow=samples-burnin+1)
  
  attach(log)
  
  S0 = get(names(log)[which(regexpr("S0\\w*", names(log))>0)])
  
  allS[,1]=S0[burnin:samples]
  
  finalS = matrix(NA,nrow=(samples-burnin+1),ncol=5)
  finalI = matrix(NA,nrow=(samples-burnin+1),ncol=5)
  finalR = matrix(NA,nrow=(samples-burnin+1),ncol=5)
  
  N= median(S0[burnin:samples])
  T = get(names(log)[which(regexpr("origin\\w*", names(log))>0)])
  
  # Calculations for plotting trajectories:	
  
  S = I = inc = R = t = rep(NA,intervals)
  S[1] = median(S0[burnin:samples])
  I[1] = 1
  inc[1] = 1
  R[1] = 0
  t[1] = finalSampleTime - median(T[burnin:samples])
  step = median(T[burnin:samples])/intervals
  
  for (j in 2:intervals){
    dS = get(names(log)[which(grepl("dS",names(log)))[j-1]])
    dR = get(names(log)[which(grepl("dR",names(log)))[j-1]])
    S[j] =S[j-1] - median(dS[burnin:samples])
    R[j] = R[j-1] + median(dR[burnin:samples])
    I[j] = I[j-1] + median(dS[burnin:samples]) - median(dR[burnin:samples])
    inc[j] = median(dS[burnin:samples])
    t[j] = t[j-1] + step
  }
}

# Create dataframe
df <- data.frame(t, S, I, R)
df$t <- as.Date((df$t-2020)*366, origin="2020-01-01")

# Define scale
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Define Squish Function
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # Get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # Apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # Get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # Apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # Return the transformation
  return(trans_new("squished", trans, inv))
}

# Transform data to long format
df.long <- df %>% gather(comp, value, -c(t))
df.long$comp[df.long$comp == "S"] <- "Suscepptible"
df.long$comp[df.long$comp == "I"] <- "Infected"
df.long$comp[df.long$comp == "R"] <- "Removed"
df.long$comp_f = factor(df.long$comp, levels=c('Suscepptible','Infected','Removed'))

# Plot
p <- ggplot(data = df.long) +
  geom_line(aes(x=t, y=value), size = 1.25) +
  guides(color=FALSE) +
  facet_grid(comp_f~., scales="free_y") +
  theme_classic() +
  labs(x="Date", y="Number of Individuals") +
  theme(panel.background = element_rect(fill = 'whitesmoke'), 
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  scale_x_date(date_breaks = "1 month", 
               labels=date_format("%b-%Y")) +
  scale_y_continuous(label = comma)

# Save Plot
ggsave(plot = p,
       filename = "ncr-bdsir-traj.png",
       width = 8, height = 5, units = "in", dpi = 300)

# Get infected values
df.infected <- subset(df.long, comp == "Infected")
df.infected$comp_f <- NULL
df.infected$comp <- NULL
colnames(df.infected) <- c("bdsir_date", "bdsir_infected")
df.infected <- subset(df.infected, bdsir_date >= as.Date("2020-03-01"))

# Save infected values 
fwrite(df.infected, "ncr-bdsir-infected.csv", na = "NA") 

