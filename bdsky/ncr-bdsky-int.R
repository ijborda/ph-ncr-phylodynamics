# Library
library(beastio)
library(ggplot2)
library(tidyr)
library(lubridate) 
library(scales)
library(tidyverse)
library(forcats)
library(hrbrthemes)
library(viridis)
library(see)

# Read logfiles
res1 = readLog(filename = "ncr-bdsky-int-post.log", 
              burnin = 0.1, 
              maxsamples = -1, 
              as.mcmc =FALSE, 
              burninAsSamples = FALSE)
res2 = readLog(filename = "ncr-bdsky-int-prior.log", 
              burnin = 0.1, 
              maxsamples = -1, 
              as.mcmc =FALSE, 
              burninAsSamples = FALSE)

# Create dataframe
data1 <- res1[, c("reproductiveNumber_BDSKY_Serial.1", "reproductiveNumber_BDSKY_Serial.2", "reproductiveNumber_BDSKY_Serial.3", "reproductiveNumber_BDSKY_Serial.4")]
colnames(data1) <- c("None.post", "ECQ.post", "MECQ.post", "GCQ.post")
data2 <- res2[, c("reproductiveNumber_BDSKY_Serial.1", "reproductiveNumber_BDSKY_Serial.2", "reproductiveNumber_BDSKY_Serial.3", "reproductiveNumber_BDSKY_Serial.4")]
colnames(data2) <- c("None.prior", "ECQ.prior", "MECQ.prior", "GCQ.prior")
data <-  data.frame(data1, data2)

# Tidy
data_tidy <- gather(data, condition, rep, None.post:GCQ.prior, factor_key=TRUE)

# Separate prior and posterior
data_tidy$distribution <- ifelse(grepl("post",data_tidy$condition),'Posterior','Prior')

# Separate intervals
data_tidy$interval <- ""
data_tidy$interval <- ifelse(grepl("None",data_tidy$condition),'None0000', data_tidy$interval)
data_tidy$interval <- ifelse(grepl("ECQ",data_tidy$condition),'ECQ', data_tidy$interval)
data_tidy$interval <- ifelse(grepl("MECQ",data_tidy$condition),'MECQ', data_tidy$interval)
data_tidy$interval <- ifelse(grepl("GCQ",data_tidy$condition),'GCQ', data_tidy$interval)
data_tidy$interval <- ifelse(grepl("None0000",data_tidy$interval),'None', data_tidy$interval)


# Define function for splitting violion plots, source: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# Plot violin plot
data_tidy$interval <- factor(data_tidy$interval, levels = c('None','ECQ', 'MECQ','GCQ'), ordered = TRUE)
p <- ggplot(data_tidy, aes(x=interval, y=rep, fill=distribution)) + geom_split_violin() +
  theme_classic() +
  ylab("Effective Reproductive Number") +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0,10)) +
  scale_fill_manual(values=c("darkblue", "8c8c8c")) +
  scale_x_discrete(labels = c('None \n Before Mar 17','ECQ \n Mar 17 - May 15', 'MECQ \n May 16 - 31 May','GCQ \n Jun 1 - 18 Jul')) +
  geom_hline(yintercept = 1.0, linetype = 1, color = "red") +
  theme(legend.position="top", legend.box = "horizontal", legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_blank())

# Kolmogorov-Smirnov test
pval_none <- ks.test(data_tidy[data_tidy$condition == "None.prior", ]$rep, data_tidy[data_tidy$condition == "None.post", ]$rep)
# p-value < 2.2e-16
pval_ecq  <- ks.test(data_tidy[data_tidy$condition == "ECQ.prior", ]$rep,  data_tidy[data_tidy$condition == "ECQ.post", ]$rep)
# p-value < 2.2e-16
pval_mecq <- ks.test(data_tidy[data_tidy$condition == "MECQ.prior", ]$rep, data_tidy[data_tidy$condition == "MECQ.post", ]$rep)
#p-value = 0.0215
pval_gcq  <- ks.test(data_tidy[data_tidy$condition == "GCQ.prior", ]$rep,  data_tidy[data_tidy$condition == "GCQ.post", ]$rep)
#p-value < 2.2e-16%

# Put error bars
none.prior.mid = median(data_tidy[data_tidy$condition == "None.prior", ]$rep)
none.prior.bot = quantile(data_tidy[data_tidy$condition == "None.prior", ]$rep, 0.05)
none.prior.top = quantile(data_tidy[data_tidy$condition == "None.prior", ]$rep, 0.95)
none.post.mid = median(data_tidy[data_tidy$condition == "None.post", ]$rep)
none.post.bot = quantile(data_tidy[data_tidy$condition == "None.post", ]$rep, 0.05)
none.post.top = quantile(data_tidy[data_tidy$condition == "None.post", ]$rep, 0.95)

ecq.prior.mid = median(data_tidy[data_tidy$condition == "ECQ.prior", ]$rep)
ecq.prior.bot = quantile(data_tidy[data_tidy$condition == "ECQ.prior", ]$rep, 0.05)
ecq.prior.top = quantile(data_tidy[data_tidy$condition == "ECQ.prior", ]$rep, 0.95)
ecq.post.mid = median(data_tidy[data_tidy$condition == "ECQ.post", ]$rep)
ecq.post.bot = quantile(data_tidy[data_tidy$condition == "ECQ.post", ]$rep, 0.05)
ecq.post.top = quantile(data_tidy[data_tidy$condition == "ECQ.post", ]$rep, 0.95)

mecq.prior.mid = median(data_tidy[data_tidy$condition == "MECQ.prior", ]$rep)
mecq.prior.bot = quantile(data_tidy[data_tidy$condition == "MECQ.prior", ]$rep, 0.05)
mecq.prior.top = quantile(data_tidy[data_tidy$condition == "MECQ.prior", ]$rep, 0.95)
mecq.post.mid = median(data_tidy[data_tidy$condition == "MECQ.post", ]$rep)
mecq.post.bot = quantile(data_tidy[data_tidy$condition == "MECQ.post", ]$rep, 0.05)
mecq.post.top = quantile(data_tidy[data_tidy$condition == "MECQ.post", ]$rep, 0.95)

gcq.prior.mid = median(data_tidy[data_tidy$condition == "GCQ.prior", ]$rep)
gcq.prior.bot = quantile(data_tidy[data_tidy$condition == "GCQ.prior", ]$rep, 0.05)
gcq.prior.top = quantile(data_tidy[data_tidy$condition == "GCQ.prior", ]$rep, 0.95)
gcq.post.mid = median(data_tidy[data_tidy$condition == "GCQ.post", ]$rep)
gcq.post.bot = quantile(data_tidy[data_tidy$condition == "GCQ.post", ]$rep, 0.05)
gcq.post.top = quantile(data_tidy[data_tidy$condition == "GCQ.post", ]$rep, 0.95)

p <- p + annotate("pointrange", x = 1.05, y = none.prior.mid, ymin = none.prior.bot, ymax = none.prior.top,
             colour = "green", size = .5) + 
         annotate("pointrange", x = 0.95, y = none.post.mid, ymin = none.post.bot, ymax = none.post.top,
                colour = "green", size = .5) +
         annotate("pointrange", x = 2.05, y = ecq.prior.mid, ymin = ecq.prior.bot, ymax = ecq.prior.top,
                  colour = "green", size = .5) + 
         annotate("pointrange", x = 1.95, y = ecq.post.mid, ymin = ecq.post.bot, ymax = ecq.post.top,
                colour = "green", size = .5) +
         annotate("pointrange", x = 3.05, y = mecq.prior.mid, ymin = mecq.prior.bot, ymax = mecq.prior.top,
                  colour = "green", size = .5) + 
         annotate("pointrange", x = 2.95, y = mecq.post.mid, ymin = mecq.post.bot, ymax = mecq.post.top,
                colour = "green", size = .5) +
         annotate("pointrange", x = 4.05, y = gcq.prior.mid, ymin = gcq.prior.bot, ymax = gcq.prior.top,
                  colour = "green", size = .5) + 
         annotate("pointrange", x = 3.95, y = gcq.post.mid, ymin = gcq.post.bot, ymax = gcq.post.top,
                colour = "green", size = .5)

# Save plot
ggsave(plot = p,
       filename = "ncr-bdsky-int.png",
       width = 7, height = 5, units = "in", dpi = 300)