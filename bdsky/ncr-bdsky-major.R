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

# Read log files
df.post = readLog(filename = "ncr-bdsky-major-post.log", 
              burnin = 0.1, 
              maxsamples = -1, 
              as.mcmc =FALSE, 
              burninAsSamples = FALSE)
df.prior = readLog(filename = "ncr-bdsky-major-prior.log", 
                 burnin = 0.1, 
                 maxsamples = -1, 
                 as.mcmc =FALSE, 
                 burninAsSamples = FALSE)

# Create dataframe
data1 <- df.post[, c("reproductiveNumber_BDSKY_Serial.1", "reproductiveNumber_BDSKY_Serial.2")]
colnames(data1) <- c("first.post", "second.post")
data2 <- df.prior[, c("reproductiveNumber_BDSKY_Serial.1", "reproductiveNumber_BDSKY_Serial.2")]
colnames(data2) <- c("first.prior", "second.prior")
data <-  data.frame(data1, data2)

# Tidy
data_tidy <- gather(data, condition, rep, first.post:second.prior, factor_key=TRUE)

# Separate prior and posterior
data_tidy$distribution <- ifelse(grepl("post",data_tidy$condition),'Posterior','Prior')

# Separate intervals
data_tidy$interval <- ""
data_tidy$interval <- ifelse(grepl("first",data_tidy$condition),'Before major change in Re', data_tidy$interval)
data_tidy$interval <- ifelse(grepl("second",data_tidy$condition),'After major change in Re', data_tidy$interval)

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
p <- ggplot(data_tidy, aes(x=interval, y=rep, fill=distribution)) + geom_split_violin() +
  theme_classic() +
  xlab("Date") +
  ylab("Effective Reproductive Number") +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0,10)) +
  scale_fill_manual(values=c("darkblue", "8c8c8c")) +
  geom_hline(yintercept = 1.0, linetype = 1, color = "red") +
  theme(legend.position="top", legend.box = "horizontal", legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

# Save Plot
ggsave(plot = p,
       filename = "ncr-bdsky-major.png",
       width = 7, height = 5, units = "in", dpi = 300)