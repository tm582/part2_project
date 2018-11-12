#Set working directory - LAB
setwd("/Users/tm582/Desktop/part2_project/R Introduction")

#Load Libraries
library(RColorBrewer)
library(gplots)
library(DESeq2)

#Colours loaded for heatmaps
hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

#Load data
NAME=read.table('',header=TRUE,row.names = 1)

#Set conditions & colours
conds=as.factor(pdata$condition)
cond_colours=brewer.pal(length(conds),'COLOURPALETTE')
names(cond_colours)=unique(conds)

# DESeq from HTSeq
ddsHTSeq=DESeqDataSetFromHTSeqCount(sampleTable = pdata, directory = '.', design= ~ condition)
colData(ddsHTSeq)$condition=factor(colData(ddsHTSeq)$condition, levels=levels(pdata$condition))

#DESeq from Matrix
coldata=as.data.frame(t(t(conds)))
colnames(coldata)='condition'
dds<-DESeqDataSetFromMatrix(countdata=tablename,coldata=coldata,design=~treatment)

#Data normalisation
dds=estimateDispersions(dds)

normcounts <- counts(dds, normalized=TRUE)
rawcounts=counts(dds,normalized=FALSE)
log2counts=log2(normcounts+1)

#QC
quartz()
par(mfrow=c(2,1))
barplot(colSums(counts(dds, normalized=FALSE)), col=cond_colours[pdata$Sample_Type], las=2,cex.names=0.5,main='Pre Normalised Counts')
barplot(colSums(counts(dds, normalized=TRUE)), col=cond_colours[pdata$Sample_Type], las=2,cex.names=0.5,main='Post Normalised Counts')
#Add legends if required


#Dispersion Plot
quartz()
dds=nbinomWaldTest(dds)
counts_table=counts(dds, normalized=TRUE)
plotDispEsts(dds,main='Dispersion Plot')

#VST
vsd <- varianceStabilizingTransformation(dds)
vstMat <- assay(vsd)
vstcounts <- vstMat[order(apply(vstMat,1,sum),decreasing =TRUE),]

#PCA plots
pca=princomp(assay(vsd))
plot(pca$loadings, main='PCA Variance Stabilised', pch=21, col='black', bg=cond_colours,cex=1)
text(pca$loadings, conds, pos=1, cex=0.8)

#Statistical Analysis
p_threshold=0.05
lfc_threshold=0.7

#Results & hits calcs
res=results(dds, contrast=c('condition', 'control', 'comparecondition'))
res=res[order(res$padj),]

hits=(rownames(res[((res$padj<=p_threshold) & (abs(res$log2FoldChange)>=lfc_threshold) & (!is.na(res$padj))),]))

#Median
median=apply(vstMat[rows,columns],1,median)

#Volcano Plots
plot(RESULTSDATA$log2FoldChange,-log(RESULTSDATA$padj,10), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= 'Volcano Plot: Control vs CONDITION', pch=19, cex=0.2)
text(RESULTSDATA[HITSDATA,]$log2FoldChange,-log(RESULTSDATA[HITSDATA,]$padj,10),labels=names[rownames(RESULTSDATA[HITSDATA,]),'COLUMN WITHIN DATA'],pos=3,cex=0.4)
points(RESULTSDATA[HITSDATA,'log2FoldChange'],-log(RESULTSDATA[HITSDATA,'padj'],10),pch=21,cex=0.4,col='turquoise1')

#Add horizotnal and vertical boudnary lines
abline(v=-1, lty=3)
abline(v=-1, lty=-3)
abline(v=1, lty=3)

#Scatterplots
plot(controlmedian,conditionmedian, main="Control vs. Condition", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)
points(controlmedian[hits],conditionmedian[hits], pch=21, col='turquoise1')
text(controlmedian[hits],conditionmedian[hits],cex=4,pos=4,labels = table[hits,column])

#Heatmaps
heatmap.2(vstMat[hits,],trace='none',col=hmcol,labRow=names[hits, column],cexRow=0.4,cexCol=0.6,las=2,Colv=FALSE,dendrogram='row',main=substitute(paste("Heatmap Sig. Hits for condition")))
#FOR LOOPS - results, hits, volcano, heatmap

condition_list=levels(conds)
for (i in 2:length(condition_list)){
  
  control=condition_list[1]
  compare=condition_list[i]
  print(paste("Comparing:",control," vs ",compare))
  
  rm(res)
  res=results(dds, contrast=c('condition', control, compare))
  res=res[order(res$padj),]
  
  rm(hits)
  hits=(rownames(res[((res$padj<=p_threshold) & (abs(res$log2FoldChange)>=lfc_threshold) & (!is.na(res$padj))),]))
  
  print(paste("Found: ",length(hits)," Statistical hits"))
  
  quartz()
  plot(res$log2FoldChange,-log10(res$padj), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= paste('Volcano Plot:', control, 'vs', compare), pch=19, cex=0.2)
  text(res[hits,]$log2FoldChange,-log(res[hits,]$padj,10),labels=names[hits,]$V2,pos=3,cex=0.4)
  points(res[hits,'log2FoldChange'],-log(res[hits,'padj'],10),pch=21,cex=0.4,col='turquoise1')
  abline(h=-log10(0.05), lty=3)
  abline(v=-1, lty=3)
  abline(v=1, lty=3)
  
  quartz()
  heatmap.2(vstMat[hits,],trace='none',col=hmcol,labRow=names[hits,'V2'],cexRow=0.4,cexCol=0.6,las=2,Colv=FALSE,dendrogram='row',main=substitute(paste("Heatmap Sig. Hits for", compare)))
}