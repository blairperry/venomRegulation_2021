library(ggplot2)
library(RColorBrewer)
library(reshape)
library(gridExtra)
library(cowplot)
#library(ape)
#library(geiger)
#library(nlme)
#library(phytools)

val <- read.csv('pi.csv')
val$Source <- factor(val$Source, levels = val$Source)


#blank theme
blank_theme <- theme_bw() + theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size = 6)) +
  theme(axis.title = element_text(size = 10))
#+
#  theme(legend.position = "none")

## Boxplot

p <- ggplot(val, aes(x=Source, y=pi))
p_II <- p + geom_boxplot() + blank_theme + coord_flip()
p_II + geom_jitter(position=position_jitter(0.15)) 
