# Libraries
library(ggplot2)
library(ggtree)
library(treeio)
library(phytools)
library(tidyr)
library(lubridate) 
library(scales)
library(tidyverse)
library(forcats)
library(hrbrthemes)
library(viridis)
library(see)

# Read tree
tree <- read.iqtree("ncr-global.treefile") 
metadata <- read.delim("ncr-global.tsv", sep="\t")                         

# Plot ML tree
phylo_tree <- ggtree(tree, root.position = 1) %<+% metadata +      
  geom_tippoint(aes(color=area),                            
                alpha=.9,                                       
                size=3) +                                       
  scale_color_brewer("island",                                
                     palette="Paired",                         
                     direction = -1) +                          
  guides(color = guide_legend(override.aes = list(size = 8))) + 
  geom_nodelab(size = 2)
phylo_tree

# Save ML tree
ggsave(plot = phylo_tree,
       filename = "ncr-global.png",
       width = 5.5, height = 20, units = "in", dpi = 500)
