setwd("/Macintosh HD/Users/Tara/Desktop/part2_project/Spermatogenesis")
library(RColorBrewer)
library(gplots)
library(DESeq2)

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

#load count data
names=read.table(file.choose())

#load pdata
pdata=read.table(file.choose(),row.names=1,header=TRUE)
colnames(pdata)
rownames(pdata)

#set conditions
conds=as.vector(pdata$condition)
#Set colours
cond_colours = brewer.pal(length(unique(conds)),"Accent")

