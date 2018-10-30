#Set working directory - LAB
setwd("/Users/tm582/Desktop/part2_project/R Introduction")

#Load Libraries
library(RColorBrewer)
library(gplots)
library(DESeq2)

#Colours loaded for heatmaps
hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

# DESeq from HTSeq
ddsHTSeq=DESeqDataSetFromHTSeqCount(sampleTable = pdata, directory = '.', design= ~ condition)
colData(ddsHTSeq)$condition=factor(colData(ddsHTSeq)$condition, levels=levels(pdata$condition))

#Data normalisation
dds=estimaddteSizeFactors(ddsHTSeq)
dds=estimateDispersions(dds)

normcounts <- counts(dds, normalized=TRUE)
rawcounts=counts(dds,normalized=FALSE)
log2counts=log2(normcounts+1)

#VST
vsd <- varianceStabilizingTransformation(dds)
#from matrix
vstMat <- assay(vsd)
vstcounts <- vstMat[order(apply(vstMat,1,sum),decreasing =TRUE),]

#PCA plots
plot(pca$loadings, main='PCA Variance Stabilised', pch=21, col='black', bg=cond_colours,cex=1)
text(pca$loadings, conds, pos=1, cex=0.8)


#Volcano Plots
plot(RESULTSDATA$log2FoldChange,-log(RESULTSDATA$padj,10), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= 'Volcano Plot: Control vs CONDITION', pch=19, cex=0.2)
text(RESULTSDATA[HITSDATA,]$log2FoldChange,-log(RESULTSDATA[HITSDATA,]$padj,10),labels=names[rownames(RESULTSDATA[HITSDATA,]),'COLUMN WITHIN DATA'],pos=3,cex=0.4)
points(RESULTSDATA[HITSDATA,'log2FoldChange'],-log(RESULTSDATA[HITSDATA,'padj'],10),pch=21,cex=0.4,col='turquoise1')
#Add horizotnal and vertical boudnary lines
abline(v=-1, lty=3)
abline(v=-1, lty=-3)
abline(v=1, lty=3)