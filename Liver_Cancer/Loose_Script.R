
##Save Plots in PDF for TCGA Analysis

pdf('TCGA_Data_Plots.pdf', paper = 'a4')

barplot(colSums(counts(dds, normalized=FALSE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Pre Normalised Counts', names.arg =conds)
legend("topleft",levels(conds),cex=0.5,fill=cond_colours)

barplot(colSums(counts(dds, normalized=TRUE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Post Normalised Counts', names.arg =conds)
legend("topleft",levels(conds),cex=0.5,fill=cond_colours)


plotDispEsts(dds,main='Dispersion Plot')


plot(pca$x, main='PCA Variance Stabilised', pch=21, col='black', bg=cond_colours[conds],cex=1)
text(pca$x, as.character(conds), pos=1, cex=0.4)


pca2d(pca, group=conds, col=cond_colours, legend='bottomright')


heatmap.2(cor(vstcounts),trace="none",labRow = conds, cexRow = 0.6, labCol = pdata$Sample_Type, cexCol = 0.6, col=hmcol, main = 'Heatmap of Correlation of the VST Counts')


ggplot(data = tsne.df, aes(tsne.1,  tsne.2, colour=Type, shape=Type))+
  geom_point(size=4)+
  scale_color_manual(values=unique(cond_colours[tsne.df$Type]))

for (i in 2:(length(condition_list))){
  
  normal=condition_list[1]
  compare=condition_list[i]
  
  print(paste('Comparing', normal, 'vs.', compare))
  
  res=results(dds, contrast = c('Sample_Type', normal, compare))
  res=res[order(res$padj),]
  
  hits=(rownames(res[((res$padj<=p_threshold) & (abs(res$log2FoldChange)>=lfc_threshold) & (!is.na(res$padj))),]))
  
  print(paste("Found: ",length(hits)," Statistical hits"))
  
  if (length(hits) != 0){
    if (length(hits)>max_display_hits){
      hits=hits[1:max_display_hits]
    }
    
    plot(res$log2FoldChange,-log10(res$padj), ylab='-log10(Adjusted P-value)', xlab='log2 FoldChange', main=paste('Volcano Plot:', normal, 'vs.', compare), pch=19, cex=0.2)
    points(res[hits,'log2FoldChange'],-log(res[hits,'padj'],10),pch=21,cex=0.4,col='turquoise1')
    text(res[hits,]$log2FoldChange,-log10(res[hits,]$padj), labels = gene_names[hits,]$GeneName, pos=3, cex = 0.4)
    abline(h=-log10(0.05), lty=3, col='red')
    abline(v=-1, lty=3, col='red')
    abline(v=1, lty=3, col='red')
    
    heatmap.2(vstMat[hits,],trace='none',col=hmcol,Colv=FALSE,dendrogram='row',main=paste("Heatmap Sig. Hits for", compare),labCol = conds, cexCol = 0.5) 
  } else {
    print(paste("No Results for ",compare))
  }
  
}


plot(median.normal,median.primary, main="Normal vs. Primary_Tumor", pch=19, col='darkblue', cex=0.4)
points(median.normal[hits],median.primary[hits], pch=21, col='turquoise1')
text(median.normal[hits],median.primary[hits],cex=0.4,pos=4,labels = gene_names[hits,]$GeneName)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)


plot(median.normal,median.recurrent, main="Normal vs. Recurrent_Tumor", pch=19, col='darkblue', cex=0.4)
points(median.normal[hits],median.recurrent[hits], pch=21, col='turquoise1')
text(median.normal[hits],median.recurrent[hits],cex=0.4,pos=4,labels = gene_names[hits,]$GeneName)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)

dev.off()

--------------------------------------
  
#Princomp PCA
pca2=princomp(assay(vsd))
quartz()
plot(pca2$loadings, main='PCA Variance Stabilised for Merged Data', pch=21, col='black', bg=cond_colours[conds], cex=1)
text(pca2$loadings, as.character(conds), pos=1, cex=0.4)

summary(pca2)


--------------------------------------

##Old data labelling
#OLD CODE
colnames(mergedcounts)[1:20]='GTEX_Normal'
colnames(mergedcounts)[21:32]=paste(pdata$Sample_Type)

#SHORT WAY
conds_merged=as.factor(colnames(mergedcounts))
conds_merged_colours=brewer.pal(length(unique(conds_merged)), 'Set1')
names(conds_merged_colours)=unique(conds_merged)

#LONG WAY
conds_merged2=as.factor(c('GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','GTEX_Normal','Primary_Tumor','Primary_Tumor','Primary_Tumor','Primary_Tumor','Recurrent_Tumor', 'Recurrent_Tumor', 'Recurrent_Tumor', 'Solid_Tissue_Normal','Solid_Tissue_Normal','Solid_Tissue_Normal','Solid_Tissue_Normal', 'Solid_Tissue_Normal'))

conds_merged_colours2=brewer.pal(length(unique(conds_merged2)), 'Set1')
names(conds_merged_colours2)=unique(conds_merged2)

mergeddata=read.table('mergeddata.txt', header = T, row.names = 1)



------------------------------------
  

##TEST MATRIX CODE
##IGNORE TEST MATRIX CODE
NEW APPROACH TAKEN FROM HERE - Decided that expression levels in TCGA are too low to directly compare with GTEX and DESeq is insufficient to normalise. Alternative normalisation appraoches taken

A new matrix (testmatrix) is created excluding coutns that are very low - reduces the number of genes)

#TEST WITH NON ZERO GENE COUNTS
#New fucntion created to test merged coutns matrix for non zero values
```{r}
nonzero <- function(val) {
  return(val>=1)
}

#Test matrix created containing only nonzero gene counts based on above function
testmatrix=mergedcounts[apply(mergedcounts,1,min) > 1,]
```


#DESeq normalisation - not including batch effects
Basic code used in all appraoches outlined below
```{r}
testHTSeq<-DESeqDataSetFromMatrix(countData=testmatrix,colData=coldata,design=~Sample_Type)

testdds=estimateSizeFactors(testHTSeq)
testdds=estimateDispersions(testdds)

testnormcounts <- counts(testdds, normalized=TRUE)
testrawcounts=counts(testdds,normalized=FALSE)
testlog2counts=log2(testnormcounts+1)

#Plot to check there is ok distribution of gene counts after removing low gene counts
quartz()
plot(testmatrix[,1],testmatrix[,2])
quartz()
boxplot(log2(testmatrix+1), las=2, col=cond_colours[conds])

#Plot normalised counts - Test matrix produces same bar plot of (pre and post) normalised counts
quartz()
par(mfrow=c(1,2))
barplot(colSums(counts(testdds, normalized=FALSE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Pre Normalised Counts')
legend("topright",levels(conds),cex=0.5,fill=cond_colours[levels(conds)])

barplot(colSums(counts(testdds, normalized=TRUE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Post Normalised Counts')
legend("topright",levels(conds),cex=0.5,fill=cond_colours[levels(conds)])


#Further Normalisation & Dispersion Plot
testdds=nbinomWaldTest(testdds)
test_counts=counts(testdds, normalized=TRUE)

quartz()
plotDispEsts(testdds,main='Dispersion Plot')

#VSD
testvsd <- varianceStabilizingTransformation(testdds)
testvstMat <- assay(testvsd)
testvstcounts <- testvstMat[order(apply(testvstMat,1,sum),decreasing =TRUE),]
testpca=prcomp(t(assay(testvsd)))


The normalisation is generally poor using DESeq functions
Decided to look at other normalisation options to see whether we can compare the GTEX/TCGA batch effect


#Normal DESeq Normalisation and Analysis using test DESeq objects above - USING TEST MATRIX
```{r}
#VSD
testvsd <- varianceStabilizingTransformation(testdds)
testvstMat <- assay(testvsd)
testvstcounts <- testvstMat[order(apply(testvstMat,1,sum),decreasing =TRUE),]

#PCA Plot - All same point return same PC value
testpca=prcomp(t(assay(testvsd)))

quartz()
plot(testpca$x, main='PCA Variance Stabilised for Merged Data', pch=21, col='black', bg=cond_colours[conds], cex=1)
text(testpca$x, as.character(conds), pos=1, cex=0.4)

summary(testpca)

#Heatmap
quartz()
heatmap.2(cor(test_counts),trace="none",labRow = coldata$Sample_Type, cexRow = 0.6, labCol = coldata$Sample_Type, cexCol = 0.6, col=hmcol)

#Rtsne
set.seed(72)
testtsne = Rtsne(t(testvstMat), perplexity = 3, check_duplicates=F)
testtsne.df=data.frame(testtsne.1=testtsne$Y[,1], testtsne.2=testtsne$Y[,2], Type=as.factor(testdds$Sample_Type))

quartz()
ggplot(data = testtsne.df, aes(testtsne.1, testtsne.2, colour=Type, shape=Type))+
  geom_point(size=4)+
  scale_color_manual(values=unique(cond_colours[testtsne.df$Type]))



-------------------------------------------------
  
##LIMMA

#Limma function removal of batch effect
assay(mergevsd)=limma::removeBatchEffect(assay(mergevsd), mergevsd$batch)
quartz()
plotPCA(mergevsd, 'Batch')

assay(mergevsd)=limma::removeBatchEffect(assay(mergevsd), mergevsd$Sample_Type)
quartz()
plotPCA(mergevsd, 'Sample_Type')


filter1=function(x)(IQR(x)>0.5)

filteredcounts=mergevstcounts[genefilter(mergevstcounts,filter1),]

pca=prcomp(t(filteredcounts))
quartz()
plot(pca, type='l')

quartz()
plot(pca$x, main='PCA Variance Stabilised for Merged Data', pch=21, col='black', bg=cond_colours[conds], cex=1)
text(pca$x, as.character(conds), pos=1, cex=0.4)

-------------------------------------------------
  
#SVA Normalisation
  
#Need to create pdata-esque file for SVA
pdata=coldata

mod=model.matrix(~as.factor(conds), data = mergedcounts)
mod0=model.matrix(~1, data = pdata)
nsv=num.sv(rawcounts, mod, method = 'leek')
sva=sva(rawcounts, mod, mod0)

#Adjusting for Surrogate Variables
pValues=f.pvalue(test_counts,mod,mod0)
qValues=p.adjust(pValues, method='BH')

modSv = cbind(mod,sva$sv)
mod0Sv = cbind(mod0,sva$sv)
pValuesSv = f.pvalue(test_counts,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

#ComBat Function to adjust for batch effect between GTEX and TCGA - TO RETURN TO
batch
modcombat=model.matrix(~1, data = testmatrix)
combat_test=ComBat(dat = test_counts, batch = batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)

#Simplified ComBat approach - TO RETURN TO 
modBatch=model.matrix(~as.factor(conds) + as.factor(batch), data = testmatrix)
mod0Batch=model.matrix(~as.factor(batch), data = testmatrix)

pValuesBatch = f.pvalue(test_counts,modBatch,mod0Batch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")

#sva for sequencing
mod1=model.matrix(~conds)
modnull=cbind(mod1[,1])
svseq=svaseq(test_counts, mod1, modnull)$sv

plot(svseq, pch=19, col=cond_colours[conds])

Old TMM Normalisation Code

Type=factor(substring(colnames(mergedcounts),1))
Time=factor(substring(colnames(mergedcounts),1,))

dge=DGEList(counts=mergedcounts)

dge$samples$Type=conds
dge$samples$Batch=batch

dge=calcNormFactors(dge, method = 'TMM')

quartz()
boxplot(log2(dge$counts+1), las=2, col=cond_colours[conds], main='log2 counts', cex.axis=0.6)

quartz()
plotMDS(dge, pch = 19, col= cond_colours[conds], main='Multi Dimensional Scaling (MDS) Plot')
legend('bottomright', legend=levels(conds), col=cond_colours[unique(conds)], pch = 20)


design=model.matrix(~conds+batch, data = dge$samples)
rownames(design)=colnames(dge)

dge=estimateDisp(dge, design, robust=T)
dge$common.dispersion
quartz()
plotBCV(dge)

dge=estimateGLMCommonDisp(dge)
dge=estimateTagwiseDisp(dge)
fit = glmFit(dge, design, robust = T)

quartz()

fit2=glmQLFit(dge, design)
Compare GTEX normal to Primary Tumour (2), Recurrent Tumour (3), Solid 
qlf=glmQLFTest(fit2, coef = 2)
qlf=glmQLFTest(fit2, coef = 3)
qlf=glmQLFTest(fit2, coef = 4)

-----------------------------------
  
#Edge R analysis with only conditions 
  #Design only include conditions
  ```{r}

#Design
design=model.matrix(~0+conds)
contr.matrix=makeContrasts(condsPrimary_Tumor-condsNormal, condsRecurrent_Tumor-condsNormal, levels = colnames(design))

quartz()
v=voom(dge, design, plot = T)

vfit=lmFit(v, design)
vfit=contrasts.fit(vfit, contrasts = contr.matrix)
efit=eBayes(vfit)
contrast1=efit


#PCA
voompca=prcomp(t(v$E))

quartz()
plot(voompca$x, pch=c(21,24)[batch], bg=cond_colours[conds], cex=1)
text(voompca$x,as.character(coldata$Batch),cex=0.3,pos=1)
text(voompca$x,as.character(colnames(v$E)),cex=0.3,pos = 3)


#SA Plot
quartz()
plotSA(efit, main='Mean-Variance Trend')


#MDS Plot
quartz()
plotMDS(v, pch = 19, col= cond_colours[conds], main='Multi Dimensional Scaling (MDS) Plot')
legend('bottomleft', legend=levels(conds), col=cond_colours[unique(conds)], pch = 20)


#Summary
summary(decideTests(efit))


#Volcano Plot
quartz()
volcanoplot(efit, coef=1,names=as.character(gene_names[rownames(efit),1]), highlight = 50, main='Primary Tumour vs. Normal')
abline(h=-log10(0.05), lty=3, col='red')
abline(v=-1, lty=3, col='red')
abline(v=1, lty=3, col='red')

quartz()
volcanoplot(efit, coef=2, names=as.character(gene_names[rownames(efit),1]), highlight = 50, main='Recurrent Tumour vs. Normal') 
abline(h=-log10(0.05), lty=3, col='red')
abline(v=-1, lty=3, col='red')
abline(v=1, lty=3, col='red')


#Further Analysis
tfit = treat(vfit, lfc=1)
dt = decideTests(tfit)

summary(dt)

-------------------------------------------------------------------------
  
#Visualising top gene count distribution
topGene=rownames(rltres)[which.min(rltres$padj)]

quartz()
plotCounts(mergedds, gene=topGene, intgroup=c("Sample_Type"), col=cond_colours[conds], pch=c(21,24)[batch])

TopSigGenes=rownames(rltres, 100)
mat <- assay(rlt)[TopSigGenes, ]
mat <- mat - rowMeans(mat)

-----------------------------------------------------------------------

  #Rtsne SVA
  
  ```{r}
set.seed(72)
svatsne = Rtsne(t(assay(ddssva)), perplexity = 3)
svatsne.df=data.frame(svatsne.1=svatsne$Y[,1], svatsne.2=svatsne$Y[,2], Type=as.factor(ddssva$Sample_Type), Batch=as.factor(ddssva$Batch))

quartz()
ggplot(data = svatsne.df, aes(svatsne.1,  svatsne.2, colour=Type, shape=Batch))+
  geom_point(size=4)+
  scale_color_manual(values=unique(cond_colours[svatsne.df$Type]))+
  ggtitle('Rtsne for SVA')


------------------------------------------------------------------------

  #Rstne - rlog data
  set.seed(72)
rlttsne = Rtsne(t(assay(rlt)), perplexity = 3)
rlttsne.df=data.frame(rlttsne.1=rlttsne$Y[,1], rlttsne.2=rlttsne$Y[,2], Type=as.factor(mergedds$Sample_Type), Batch=as.factor(mergedds$Batch))

quartz()
ggplot(data = rlttsne.df, aes(rlttsne.1,  rlttsne.2, colour=Type, shape=Batch))+
  geom_point(size=4)+
  scale_color_manual(values=unique(cond_colours[rlttsne.df$Type]))

------------------------------------------------------------------------

#Key TCGA code
  Working Directory (Lab) & Libraries
``` {r echo=FALSE, r eval=F}
setwd("~/Desktop/part2_project/Liver_Cancer/TCGA_Counts_mRNA_genelevel")

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(DESeq2)
library(Rtsne)
library(pca3d)
library(limma)
library(genefilter)
library(sva)
library(edgeR)
library(statmod)
```

#Upload pdata (TCGA)
```{r}
pdata=read.table('pdata_rnaseq_genelevel_small.txt', header=T)
rownames(pdata)=pdata$Sample
```

#Read in count data
Data taken from TCGA data repository - Liver Cancer mRNA
datatype=HTSeq count data

```{r}

##TCGA DATA COUNTS
```{r}
ddsHTSeq=DESeqDataSetFromHTSeqCount(sampleTable = pdata, directory= '.', design = ~ Sample_Type)

dds=estimateSizeFactors(ddsHTSeq)
dds=estimateDispersions(dds)

tcgarawcounts=counts(dds,normalized=FALSE)
```

# Load in Gene Names and match up with rownames in ddsHTSeq
```{r}
gene_names = read.table("gene_names.txt",row.names=1,header = F)
gene_names= gene_names[rownames(ddsHTSeq),]
colnames(gene_names)=c("GeneName","Source","BioType")
```

---------------------------------------------------------------
#For scaling up to full TCGA data set
  
pdata=read.table('pdata_rnaseq_genelevel_small.txt', header=T)
rownames(pdata)=pdata$Sample

----------------------------------------------------------------
  #Medians - SVA 
  
  ```{r}

median.primary=apply(mergevstMat[,21:24],1,median)
median.recurrent=apply(mergevstMat[,25:27],1,median)
median.normal=apply(mergevstMat[,c(1:20, 28:32)], 1, median)

quartz()
plot(median.normal,median.primary, main="Normal vs. Primary_Tumor (SVA)", pch=19, col='darkblue', cex=0.4)
points(median.normal[hits],median.primary[hits], pch=21, col='turquoise1')
text(median.normal[hits],median.primary[hits],cex=0.4,pos=4,labels = gene_names[hits,]$GeneName)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)

quartz()
plot(median.normal,median.recurrent, main="Normal vs. Recurrent_Tumor (SVA)", pch=19, col='darkblue', cex=0.4)
points(median.normal[hits],median.recurrent[hits], pch=21, col='turquoise1')
text(median.normal[hits],median.recurrent[hits],cex=0.4,pos=4,labels = gene_names[hits,]$GeneName)
abline(a=0,b=1, col='red', lwd=2, cex=0.4, lty=2)

```
----------------------------------------------------------------
  
  ##Rlog Transformation (Instead of VST)
  
  #Rlog
  
  #Starting from merged counts matrix
  #Convert to DESeq object
  mergedds=DESeqDataSetFromMatrix(mergedcounts, colData=coldata, design = ~Sample_Type + Batch)


#Remove low counts
nrow(mergedds)
mergedds=mergedds[rowSums(counts(mergedds))>1,]
nrow(mergedds)


#rlog transformation
rlt=rlog(mergedds, blind=F)

quartz()
par(mfrow=c(1,2))
plot(log2(counts(mergedds, normalized=FALSE)+ 1), pch=16, cex=0.3, main = 'Pre Normalised Counts')
plot(assay(rlt), pch=16, cex=0.3, main='rlog Normalised Counts')

quartz()
par(mfrow=c(2,2))
boxplot(log2(counts(mergedds, normalized=FALSE)+ 1), pch=16, cex=0.3, main = 'Pre Normalised Counts', col=cond_colours[conds], las=2, cex.axis=0.5)
barplot(log2(counts(mergedds, normalized=FALSE)+ 1), pch=16, cex=0.3, main = 'Pre Normalised Counts', col=cond_colours[conds], las=2, cex.axis=0.5)
boxplot(assay(rlt), col=cond_colours[conds], las=2, cex.axis=0.5, main='rlog Normalised Counts', pch=16, cex=0.3)
barplot(assay(rlt), col=cond_colours[conds], las=2, cex.axis=0.5, main='rlog Normalised Counts', pch=16, cex=0.3)

#Sample Distances - rlog
sampleDists=dist(t(assay(rlt)))
sampleDistMatrix <- as.matrix(sampleDists)

quartz()
heatmap.2(sampleDistMatrix, cexRow = 0.5, cexCol = 0.5, col = hmcol, trace = 'none', main = 'Sample-to-Sample Distances Using Rlog Values')


#PCA Plot
rltpca=prcomp(t(assay(rlt)))
quartz()
plot(rltpca$x, pch=c(21,24)[batch], bg=cond_colours[conds], main='PCA Plot: rlog values')
text(rltpca$x,as.character(coldata$Batch),cex=0.3,pos=1)
text(rltpca$x,as.character(coldata$Sample_Type),cex=0.3,pos = 3)

summary(rltpca)

#MDS Plot
mdsData=data.frame(cmdscale(sampleDistMatrix))
mds=cbind(mdsData, as.data.frame(colData(rlt)))

quartz()
ggplot(mds, aes(X1,X2,color=Sample_Type,shape=Batch))+
  geom_point(size=3)+ 
  ggtitle('MDS Plot - Rlog')+
  scale_color_manual(values=unique(cond_colours[mds$Sample_Type]))

#Differential Expression Analysis
DESeqdds=DESeq(mergedds)

rltres=results(DESeqdds)
summary(rltres)

sum(rltres$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(rltres$pvalue))

sum(rltres$padj < 0.1, na.rm=TRUE)

rltres=rltres[order((rltres$padj)),]

heatmap.2(mat, trace='none', cexRow = 0.5, cexCol = 0.5, col=hmcol)
heatmap.2(assay(rlt)[rownames(rltres[1:100,]),],col=hmcol,trace="none",cexRow = 0.5,cexCol = 0.5, main='Heatmap of top 100 genes - Rlog')


```{r}
save.image("~/Desktop/part2_project/Liver_Cancer/merged_rlog.RData")
```

----------------------------------------------------------------
  
  # Extract the best 100 (by p-value) using the toptable command over all tests (i.e. coef 1 and 2)
  #primarypval=as.data.frame(efit$p.value[,1])
  #primarysig=as.data.frame(primarypval[primarypval<=0.05])
  #rownames(primarysig)=rownames(primarypval)[primarypval<=0.05]
  
  #primarysig=apply(primarysig,2,sort)
  #primarysig=as.data.frame(primarysig[1:100,])
  #colnames(primarysig)=c('p-value')
  
  #rownames(primarysig)=gene_names[rownames(primarysig),]$GeneName

#recurrentpval=as.data.frame(efit$p.value[,2])
#recurrentsig=as.data.frame(recurrentpval[recurrentpval<=0.05])
#rownames(recurrentsig)=rownames(recurrentpval)[recurrentpval<=0.05]

#recurrentsig=apply(recurrentsig,2,sort)
#recurrentsig=as.data.frame(recurrentsig[1:100,])
#colnames(recurrentsig)=c('p-value')

#rownames(recurrentsig)=gene_names[rownames(recurrentsig),]$GeneName


tfit=treat(efit, lfc = log2(0.7))


#lfcprim=topTreat(tfit,coef = 1, number=100000000000000)

#primsig=lfcprim[lfcprim$adj.P.Val<=0.05,]
#primsig=primsig[abs(primsig$logFC)>=0.7,]
#primsig$logFC=abs(primsig$logFC)


#primsig=subset(primsig, select=c('logFC'))
#primsig=apply(primsig,2,sort, decreasing=TRUE)

#topprimsig=as.data.frame(primsig[1:100,])
#rownames(topprimsig)=gene_names[rownames(topprimsig),]$GeneName


#lfcrec=topTreat(tfit,coef = 2, number=100000000000000)

#recsig=lfcrec[lfcrec$adj.P.Val<=0.05,]
#recsig=recsig[abs(recsig$logFC)>=0.7,]
#recsig$logFC=abs(recsig$logFC)


#recsig=subset(recsig, select=c('logFC'))
#recsig=apply(recsig,2,sort, decreasing=TRUE)

#toprecsig=as.data.frame(recsig[1:100,])
#rownames(toprecsig)=gene_names[rownames(toprecsig),]$GeneName



#hits=rownames(topTable(efit,number=100))



#Compare Primary Tumour & Recurrent Tumour to TCGA Normals and to GTEX Normals Separately
```{r}

#Readjust Conditions Vectors to distinguish between GTEX and TCGA Normals - This will allow us to compare effect of comparing each individually
conds_origin=as.character(conds)
conds_origin[1:175]=rep("Normal_GTEX",175)
conds_origin[547:596]=rep("Normal_TCGA",58)
conds_origin=as.factor(conds_origin)


#Design changes to reflect changes in conds > conds_origin
design2=model.matrix(~0+conds_origin+batch)
contr.matrix2=makeContrasts(conds_originPrimary_Tumor-conds_originNormal_TCGA, conds_originPrimary_Tumor-conds_originNormal_GTEX, levels = colnames(design2))

quartz()
v2=voom(dge, design, plot = T)

vfit2=lmFit(v2, design2)
vfit2=contrasts.fit(vfit2, contrasts = contr.matrix2)
efit2=eBayes(vfit2)



#PCA Plot
voompca2=prcomp(t(v2$E))
quartz()
plot(voompca2$x, pch=c(21,24)[batch], bg=cond_colours[conds], cex=1, main = 'PCA Plot - EdgreR (GTEX/TCGA normals = Separate)')
text(voompca2$x,as.character(coldata$Batch),cex=0.3,pos=1)
text(voompca2$x,as.character(colnames(v2$E)),cex=0.3,pos = 3)


#SA Plot
quartz()
plotSA(efit2, main='Mean-Variance Trend - EdgreR (GTEX/TCGA normals = Separate)')


#MDS Plot
quartz()
plotMDS(efit2, pch = 19, col= cond_colours[conds], main='Multi Dimensional Scaling (MDS) Plot - EdgreR (GTEX/TCGA normals = Separate)')
legend('bottomright', legend=levels(conds), col=cond_colours[unique(conds)], pch = 20)


#Summary of findings from comparisons
summary(decideTests(efit2))

#Plot similarities between comparisons
plot(topTable(efit2,coef=1,number=1000000000000000)$logFC,topTable(efit2,coef=2,number=1000000000000000)$logFC, title(main='Butterfly Plot for Primary Tumour vs. TCGA Normals'))
abline(h=0)
abline(v=0)
abline(a=0,b=1, lty=1, col='red')

plot(topTable(efit2,coef=2,number=1000000000000000)$logFC,topTable(efit2,coef=4,number=1000000000000000)$logFC, title(main='Butterfly Plot for Primary Tumour vs. GTEX Normals'))
abline(h=0)
abline(v=0)
abline(a=0,b=1, lty=1, col='red')
```