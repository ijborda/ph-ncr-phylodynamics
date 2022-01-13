# Library
library(ggtree)
library(treeio)

# Set seed
set.seed(12345)

#source confounding effect scripts
source("./functions/mantel_counding.r")
source("./functions/rand_regression.r")
source("./functions/rand_xmls.r")

# Read files
tree <- read.tree("ncr.treefile")
df <- read.delim("genetic-distance.csv", stringsAsFactors = FALSE, sep = ",")

# Data frame for tree information
dt <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
dt <- merge(x = dt, y = df, all.x = TRUE, all.y = FALSE, by  = "tip", sort = FALSE)

# Root-to-tip against sampling date regression on a dated phylogeny
r <- pathogen.permutation.test(tree, dt$date, auto.round.dates=F, rounded.dates= NULL,  reroot=T, stat='nve-rms', nreps=1000, use.clusters = F, print.progress=T, output.rooted.tree=F)
