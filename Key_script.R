#Set working directory - LAB
setwd("/Users/tm582/Desktop/part2_project/R Introduction")

#Load Libraries
library(RColorBrewer)
library(gplots)
library(DESeq2)

#Colours loaded for heatmaps
hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
