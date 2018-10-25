setwd("/Macintosh HD/Users/Tara/Desktop/part2_project/Spermatogenesis")
library(RColorBrewer)
library(gplots)
library(DESeq2)

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

#load count data
names=read.table('gene_names.txt')

#load pdata
pdata=read.table('pdata.txt',header=TRUE)
colnames(pdata)
rownames(pdata)

#set conditions
conds=as.factor(pdata$condition)
#Set colours
cond_colours = brewer.pal(length(unique(conds)),"Accent")

#PROBLEM
> ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = pdata, directory = '"~/Desktop/part2_project/Spermatogenesis"', design= ~ condition)
Error in file(file, "rt") : cannot open the connection
In addition: Warning message:
  In file(file, "rt") :
  cannot open file '"~/Desktop/part2_project/Spermatogenesis"/tophat.C34M8ACXX_j3LT4_15s008602-1-1_Vasiliauskaite_lane115s008602.2238.counts': No such file or directory

colData(ddsHTSeq)$condition = factor(colData(ddsHTSeq)$condition,levels=levels(pdata$condition))

#further code
#Pre Normalisation
#barplot....Work out code for this
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])

dds = estimateSizeFactors(ddsHTSeq)
dds = estimateDispersions(dds)

dds = nbinomWaldTest(dds)

counts_table=counts(dds,normalized=TRUE)
plotDispEsts(dds,main='Dispersion Plot')

#Data normalisation
quartz()
par(mfrow=c(2,))
barplot(colSums(counts(dds, normalized=F)), col=cond_colours, las=2,cex.names=0.5,main='Pre Normalised Counts')

#TO EXPLAIN
plot(1, type='n', axes=F, xlab='', ylab='')
legend('center',levels(as.factor(conds)),fill=cond_colours)
barplot(colSums(counts(dds, normalized=T)), col=condcols[as.factor(conds)], las=2,cex.names=0.4,main='Post Normalised Counts')

#QC - VST
vsd = varianceStabilizingTransformation(dds)
vstMat = assay(vsd)
counts_table=data.frame(counts_table,vstMat)

#Plot Sample to Sample correlation on Heatmap
heatmap.2(cor(assay(vsd)),trace='none',main='Sample Correlation Variance Stabilised',col=hmcol,cexRow=0.5,cexCol=0.5)

#PCA Plot - to EXPLAIN
pca = princomp(assay(vsd))
plot(pca$loadings, main='PCA Variance Stabilised', col='black', bg=condcols[pdata$condition],  pch=21, cex=1)
text(pca$loadings, conds, pos=3, cex=0.5)