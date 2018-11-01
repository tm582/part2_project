#Laptop
setwd("/Macintosh HD/Users/Tara/Desktop/part2_project/Spermatogenesis")
#Lab
setwd("~/Desktop/part2_project/Spermatogenesis")

library(RColorBrewer)
library(gplots)
library(DESeq2)

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

#load count data
names=read.table('gene_names.txt')
rownames(names)=names$V1

#load pdata
pdata=read.table('pdata2.txt',header=TRUE)
colnames(pdata)
rownames(pdata)

#set conditions
conds=as.factor(pdata$condition)
conds2=as.vector(conds)

#Set colours
cond_colours = brewer.pal(length(unique(conds)),"Accent")
names(cond_colours)=unique(conds)

#Insert counts from HTSeq 
ddsHTSeq=DESeqDataSetFromHTSeqCount(sampleTable = pdata, directory = '/Users/tm582/Desktop/part2_project/Spermatogenesis', design= ~ condition)

colData(ddsHTSeq)$condition=factor(colData(ddsHTSeq)$condition, levels=levels(pdata$condition))

#pre normalisation
coldata=as.data.frame(t(t(conds)))
colnames(coldata)='condition'

#Data normalisation
dds=estimateSizeFactors(ddsHTSeq)
dds=estimateDispersions(dds)

normcounts <- counts(dds, normalized=TRUE)
rawcounts=counts(dds,normalized=FALSE)
log2counts=log2(normcounts+1)

#Raw counts Bar Plot
quartz()
par(mfrow=c(2,1))
barplot(apply(rawcounts,2,sum), las=2, col=cond_colours, main='Raw Counts', cex.names = 0.5, cex.axis = 0.8)
legend("topleft",levels((conds)),cex=0.5,fill=cond_colours)

#Normalised Bar plot
barplot(apply(normcounts,2,sum), las=2, col=cond_colours, main='Normalised Counts', cex.names=0.5, cex.axis = 0.8)
legend('topleft',levels((conds)),cex=0.5,fill=cond_colours)

#Dispersion Plot
quartz()
dds=nbinomWaldTest(dds)
counts_table=counts(dds, normalized=TRUE)
plotDispEsts(dds,main='Dispersion Plot')

#QC
#QUESTION: Is this the same as bar plots above?
quartz()
par(mfrow=c(2,1))
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
quartz()
pca=princomp(assay(vsd))
plot(pca$loadings, main='PCA Variance Stabilised', pch=21, col='black', bg=cond_colours,cex=1)
text(pca$loadings, conds2, pos=1, cex=0.8)

#Statistical Analysis

p_threshold=0.05
lfc_threshold=0.7

#Long - results calculation
res.dnmt3l=results(dds, contrast=c('condition', 'control', 'dnmt3l'))
res.dnmt3l=res.dnmt3l[order(res.dnmt3l$padj),]

res.mili=results(dds, contrast=c('condition', 'control', 'mili'))
res.mili=res.mili[order(res.mili$padj),]

res.miwi2=results(dds, contrast=c('condition', 'control', 'miwi2'))
res.miwi2=res.miwi2[order(res.miwi2$padj),]

#Long - List of Statistical hits
hits.dnmt3l=(rownames(res.dnmt3l[((res.dnmt3l$padj<=p_threshold) & (abs(res.dnmt3l$log2FoldChange)>=lfc_threshold) & (!is.na(res.dnmt3l$padj))),]))

hits.mili=(rownames(res.mili[((res.mili$padj<=p_threshold) & (abs(res.mili$log2FoldChange)>=lfc_threshold) & (!is.na(res.mili$padj))),]))

hits.miwi2=(rownames(res.miwi2[((res.miwi2$padj<=p_threshold) & (abs(res.miwi2$log2FoldChange)>=lfc_threshold) & (!is.na(res.miwi2$padj))),]))

#For loop - results, hits, volcano plots, heatmaps
#PROBLEM: labelling graphs with condition

condition_list=levels(conds)
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
  plot(res$log2FoldChange,-log10(res$padj), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= substitute(paste('Volcano Plot:', control, 'vs', compare)), pch=19, cex=0.2)
  text(res[hits,]$log2FoldChange,-log(res[hits,]$padj,10),labels=names[hits,]$V2,pos=3,cex=0.4)
  points(res[hits,'log2FoldChange'],-log(res[hits,'padj'],10),pch=21,cex=0.4,col='turquoise1')
  abline(h=-log10(0.05), lty=3)
  abline(v=-1, lty=3)
  abline(v=1, lty=3)
  
  quartz()
  heatmap.2(vstMat[hits,],trace='none',col=hmcol,labRow=names[hits,'V2'],cexRow=0.4,cexCol=0.6,las=2,Colv=FALSE,dendrogram='row',main=substitute(paste("Heatmap Sig. Hits for", compare)))

  fullres=merge(names,counts_table,by.x=1,by.y=0)
  fullres=merge(fullres, as.matrix(res), by.x=1,by.y=0)
  fullres=fullres[order(fullres$log2FoldChange,decreasing = TRUE),]
  fullres=fullres[!is.na(fullres$log2FoldChange),]
  rownames(fullres)=fullres$Row.names
  write.table(fullres,'FullRes.txt',sep='\t', quote=F)
  
  }

#Analysis of Sample Median Data
control_median=apply(vstMat[,7:9],1,median)

dnmt3l_median=apply(vstMat[,1:3],1,median)

mili_median=apply(vstMat[,4:6],1,median)

miwi2_median=apply(vstMat[,10:12],1,median)

#Plotting - Scatterplot
#PROBLEM - need to work out how to highlight outliers from HITS ^ see above
#PROBLEM - text function is not labelling the points with correct names want to label with gene names from names$V2 that correspond to hits 
#Control vs. dnmt3l
quartz()
par(mfrow=c(3,1))
plot(control_median,dnmt3l_median, main="Control vs. Dnmt3l", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)
points(control_median[hits.dnmt3l],dnmt3l_median[hits.dnmt3l], pch=21, col='turquoise1')
text(control_median[hits.dnmt3l],dnmt3l_median[hits.dnmt3l],cex=4,pos=4,labels = names[hits.dnmt3l,'V2'])

#Control vs. mili
plot(control_median,mili_median, main="Control vs. mili", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)
points(control_median[hits.mili],mili_median[hits.mili], pch=21, col='turquoise1')
text(control_median[hits.mili],mili_median[hits.mili],cex=4,pos=4,labels = names[hits.mili,'V2'])

#Control vs. miwi2
plot(control_median,miwi2_median, main="Control vs. miwi2", pch=19, col='darkblue', cex=0.4)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)
points(control_median[hits.miwi2],miwi2_median[hits.miwi2], pch=21, col='turquoise1')
text(control_median[hits.miwi2],miwi2_median[hits.miwi2],cex=4,pos=4,labels = names[hits.miwi2,'V2'])

#Volcano Plots
quartz()
par(mfrow=c(3,1))

plot(res.dnmt3l$log2FoldChange,-log(res.dnmt3l$padj,10), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= 'Volcano Plot: Control vs dnmt3l', pch=19, cex=0.2)
text(res.dnmt3l[hits.dnmt3l,]$log2FoldChange,-log(res.dnmt3l[hits.dnmt3l,]$padj,10),labels=names[hits.dnmt3l,]$V2,pos=3,cex=0.4)
points(res.dnmt3l[hits.dnmt3l,'log2FoldChange'],-log(res.dnmt3l[hits.dnmt3l,'padj'],10),pch=21,cex=0.4,col='turquoise1')
abline(h=-log10(0.05), lty=3)
abline(v=-1, lty=3)
abline(v=1, lty=3)

plot(res.mili$log2FoldChange,-log(res.mili$padj,10), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= 'Volcano Plot: Control vs mili', pch=19, cex=0.2)
text(res.mili[hits.mili,]$log2FoldChange,-log(res.mili[hits.mili,]$padj,10),labels=names[hits.mili,]$V2,pos=3,cex=0.4)
points(res.mili[hits.mili,'log2FoldChange'],-log(res.mili[hits.mili,'padj'],10),pch=21,cex=0.4,col='turquoise1')
abline(h=-log10(0.05), lty=3)
abline(v=-1, lty=3)
abline(v=1, lty=3)

plot(res.miwi2$log2FoldChange,-log(res.miwi2$padj,10), ylab='-log10(Adjusted P)', xlab='Log2 FoldChange', main= 'Volcano Plot: Control vs miwi2', pch=19, cex=0.2)
text(res.miwi2[hits.miwi2,]$log2FoldChange,-log(res.miwi2[hits.miwi2,]$padj,10),labels=names[rownames(res.miwi2[hits.miwi2,]),'V2'],pos=3,cex=0.4)
points(res.miwi2[hits.miwi2,'log2FoldChange'],-log(res.miwi2[hits.miwi2,'padj'],10),pch=21,cex=0.4,col='turquoise1')
abline(h=-log10(0.05), lty=3)
abline(v=-1, lty=3)
abline(v=1, lty=3)

#Heatmaps
quartz()
heatmap.2(vstMat[hits.dnmt3l,],trace='none',col=hmcol,labRow=names[hits.dnmt3l,'V2'],cexRow=0.4,cexCol=0.6,las=2,Colv=FALSE,dendrogram='row',main='Heatmap Sig. Hits for dnmt3l')

quartz()
heatmap.2(vstMat[hits.mili,],trace='none',col=hmcol,labRow=names[hits.mili,'V2'],cexRow=0.4,cexCol=0.6,las=2,Colv=FALSE,dendrogram='row',main='Heatmap Sig. Hits for mili')

quartz()
heatmap.2(vstMat[hits.miwi2,],trace='none',col=hmcol,labRow=names[hits.miwi2,'V2'],cexRow=0.4,cexCol=0.6,las=2,Colv=FALSE,dendrogram='row',main='Heatmap Sig. Hits for miwi2')
