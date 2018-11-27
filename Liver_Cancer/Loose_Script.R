
#Old data labelling
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
