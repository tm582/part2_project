setwd("/Macintosh HD/Users/Tara/Desktop/part2_project/Spermatogenesis")
library(RColorBrewer)
library(gplots)
library(DESeq2)

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

#load count data
names=read.table('gene_names.txt')

#load pdata
pdata=read.table('pdata2.txt',header=TRUE)
colnames(pdata2)
rownames(pdata2)

#set conditions
conds=as.factor(pdata2$condition)
#Set colours
cond_colours = brewer.pal(length(unique(conds)),"Accent")

#Insert counts from HTSeq 
ddsHTSeq=DESeqDataSetFromHTSeqCount(sampleTable = pdata2, directory = '/Users/tm582/Desktop/part2_project/Spermatogenesis', design= ~ condition)

#further code
colData(ddsHTSeq)$condition=factor(colData(ddsHTSeq)$condition, levels=levels(pdata$condition))

#pre normalisation
coldata=as.data.frame(t(t(conds)))
colnames(coldata)='condition'

#Data normalisation
dds=estimaddteSizeFactors(ddsHTSeq)
dds=estimateDispersions(dds)

quartz()
plotDispEsts(dds, main="Dispersion Plot")

normcounts <- counts(dds, normalized=TRUE)
rawcounts=counts(dds,normalized=FALSE)
log2counts=log2(normcounts+1)

#work out code to barplot
#Raw counts Bar Plot
barplot(apply(rawcounts,2,sum), las=2, col=cond_colours, main='Raw Counts', cex.names = 0.5)
legend("topleft",levels((conds)),cex=0.5,fill=cond_colours)
#Normalised Bar plot
barplot(apply(normcounts,2,sum), las=2, col=cond_colours, main='Normalised Counts')
legend('topleft',levels((conds)),cex=0.5,fill=cond_colours)

dds=nbinomWaldTest(dds)
counts_table=counts(dds, normalized=TRUE)
plotDispEsts(dds,main='Dispersion Plot')

#QC
quartz()
barplot(colSums(counts(dds, normalized=FALSE)), col=cond_colours, las=2,cex.names=0.5,main='Pre Normalised Counts')
barplot(colSums(counts(dds, normalized=TRUE)), col=cond_colours, las=2,cex.names=0.5,main='Post Normalised Counts')

#VST
vsd <- varianceStabilizingTransformation(dds)
vstMat <- assay(vsd)
vstcounts <- vstMat[order(apply(vstMat,1,sum),decreasing =TRUE),]

#Sample to Sample Correlation using VST
quartz()
heatmap.2(cor(assay(vsd)),trace='none',main='Sample Correlation Variance Stabilised',col=hmcol,cexRow=0.5,cexCol=0.5)

#Sample to Sample PCA
plot(pca$loadings, main='PCA Variance Stabilised', col='black', bg=cond_colours, pch=21, cex=1, text(pca$loadings, conds, pos=3, cex=0.8))
#QUESTIONS: What is 'pch', 'pos'
#PROBLEM = NEED TO LABEL POINTS WITH CONDITION NAME

#Statistical Analysis
##CONFUSED FROM HERE
res.dnmt3l=results(dds, contrast=c('condition', 'control', 'dnmt3l'))
res.mili=results(dds, contrast=c('condition', 'control', 'mili'))
res.miwi2=results(dds, contrast=c('condition', 'control', 'miwi2'))

#List of Statistical hits
#QUESTION: Do the hits have to be separate as pairwise comparisons or are they compiled together - how does this follow through when plotting...
hits.dnmt3l=(rownames(res.dnmt3l[((res.dnmt3l$padj<=0.05) & (abs(res.dnmt3l$log2FoldChange)>=1) & (!is.na(res.dnmt3l$padj))),]))
hits.mili=(rownames(res.mili[((res.mili$padj<=0.05) & (abs(res.mili$log2FoldChange)>=1) & (!is.na(res.mili$padj))),]))
hits.miwi2=(rownames(res.miwi2[((res.miwi2$padj<=0.05) & (abs(res.miwi2$log2FoldChange)>=1) & (!is.na(res.miwi2$padj))),]))

#Analysis of Sample Median Data
control_median=apply(vstMat[,7:9],1,median)
dnmt3l_median=apply(vstMat[,1:3],1,median)
mili_median=apply(vstMat[,4:6],1,median)
miwi2_median=apply(vstMat[,10:12],1,median)

#Plotting - Scatterplot
#PROBLEM - need to work out how to highlight outliers from HITS ^ see above
#Control vs. dnmt3l
plot(control_median,dnmt3l_median, main="Control vs. Dnmt3l", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)

#Control vs. mili
plot(control_median,mili_median, main="Control vs. mili", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)

#Control vs. miwi2
plot(control_median,miwi2_median, main="Control vs. miwi2", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)