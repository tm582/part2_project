
# Change working directory to the Introduction page
setwd("/Users/tm582/Desktop/part2_project/R Introduction")

#at the end of each session commit and push into Github repository
#insert data using read.table
genomes<-read.table("genomes.txt", row.names=1, header=TRUE)
#to get more info on data
dim(genomes)
#column names
colnames(genomes)
#row names
rownames(genomes)
#summarise data in a table/dataframe
summary(genomes)
#adding column of average transcript lengths of all available genomes.
genomes$Average.transcript.length<-genomes$Transcripts.length/genomes$Transcripts
#plot
barplot(genomes$Average.transcript.length, names.arg=rownames(genomes), col='lightblue',las=2)

 ##Analysis of sRNA-Seq data for miR-210 over expression exercise
#launch R studio and various libraries inc colour brewer, gplots, DeSSeq2
#define colour palette
library(RColorBrewer)
library(gplots)
library(DESeq2)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

#insert table of counts
mircounts<-read.table('mircounts.txt', header=TRUE,row.names = 1)
#load sample data (pdata)
pdata<-read.table('pdata.txt',header=TRUE,row.names = 1)

# Take the rownames from our experimental description (pdata) and make colnames for the count table
colnames(mircounts)=rownames(pdata)
groups=as.factor(pdata$group)

# Make a conditions object for DESeq
conds=pdata$treatment

##Count Loading and normalisation
#creating DESeq object from count table looking at one condiiton (treatment)
coldata=as.data.frame(conds)

# Add sample names to column data object for DESeq
rownames(coldata)=colnames(mircounts)
coldata$group = groups

# Load a DESeq count object using our counts, description and give it a experimental design.
# Simple - > Only consider the 'condition' wt or miR210
dds<-DESeqDataSetFromMatrix(mircounts,coldata,design=~conds)

# Complicated -> Also consider the group who prepared the sample.
dds<-DESeqDataSetFromMatrix(mircounts,coldata,design=~conds+group)

#pre data normalisation

# Build unique colour scheme for each condition
cond_colours = brewer.pal(length(unique(conds)),"Accent")[as.factor(conds)]
names(cond_colours)=conds

# Build unique colour scheme for each group
group_colours = brewer.pal(3,"Dark2")[as.factor(pdata$group)]
names(group_colours)=pdata$group

# Plot Before we Normalise
barplot(apply(counts(dds),2,sum), las=2,col=cond_colours,main="Pre Normalised Counts",cex.names=0.5)
legend("topleft",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])

#estimating dispersion of the data
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

plotDispEsts(dds)

# Plot the results after
barplot(apply(counts(dds,normalized=T),2,sum), las=2,col=cond_colours,main="Normalised Counts")
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])

# Try to do both plots together to compare
par(mfrow=c(2,1)) # Set up two plots
barplot(apply(counts(dds),2,sum), las=2,col=cond_colours,main="Pre Normalised Counts",cex.names=0.5)
legend("topleft",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
barplot(apply(counts(dds,normalized=T),2,sum), las=2,col=cond_colours,main="Normalised Counts")
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
par(mfrow=c(1,1))


#Post Normalisation
normcounts <- counts(dds, normalized=TRUE)
rawcounts  = counts(dds,normalized=FALSE)
log2counts = log2(normcounts+1)



#Variance Stabilising the Counts
#alternative to log2 transformation = VST transformation in DESeq2
#This allows stabilisation of varaince for low counts

vsd <- varianceStabilizingTransformation(dds)
vstcounts <- assay(vsd)
vstcounts <- vstcounts[order(apply(vstcounts,1,sum),decreasing =TRUE),]

#Doing a pairwise Pearson correlation --> plot on heatmap
heatmap.2(cor(vstcounts),trace="none",col=hmcol,main="Sample to Sample Correlation (VST)",cexRow=0.5,cexCol=0.5,RowSideColors=cond_colours, margins=c(9,7))

#PCA plot used to see if overexpression explains variance in results

pca <- princomp(vstcounts)


par(mfrow=c(2,1))
plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.4)
legend("topright",levels(conds),fill=cond_colours[levels(conds)],cex=0.4)

plot(pca$loadings, col=group_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.4)
legend("topright",levels(groups),fill=group_colours[levels(groups)],cex=0.4)

##ANALYSIS
#Lookign at difference between treatments
m210_median = apply(vstcounts[,1:3],1,median)
ctrl_median= apply(vstcounts[,4:6],1,median)


plot(ctrl_median,m210_median,pch=19,cex=0.5,col="darkblue")
points(ctrl_median["hsa-mir-210-3p"],m210_median["hsa-mir-210-3p"],col="red")
points(ctrl_median["scrambled"],m210_median["scrambled"],col="green")
abline(a=0,b=1,lty=2,col="red",cex=0.5)
text(ctrl_median[1:10],m210_median[1:10],labels=names(ctrl_median[1:10]),cex=0.4,pos=2)
legend("topleft",c("scrambled","miR-210"),fill=c("green","red"),cex=0.7)

#The function nbinomWaldTest fits a negative binomial generalized linear model to each gene and then calculates the significance of the estimated coeffients. We contrast the control and the overexpressed miR210 samples.
p_threshold=0.05
lfc_threshold=0.8

cds <- nbinomWaldTest(dds)

# Compare the two conditions of interest
res=results(cds,contrast=c("conds","miR210","Scr"))

# Resort the list
res <- res[order(res$padj),]

# Find Significant genes   PValue < threshold   Lfc > < threshold
sig = rownames(res[(abs(res$log2FoldChange) > lfc_threshold) & (res$padj < p_threshold) & !is.na(res$padj),])

#Volcano Plot allows visualisation of the differentially expressed genes
#plotting the log2 Fold-change of each gene (x-axis) against the -log10 of its p-value (y-axis).
plot(res$log2FoldChange,-log(res$padj,10),ylab="-log10(Adjusted P)",xlab="Log2 FoldChange",main=paste("Volcano Plot","ctrl v mir210"),pch=19,cex=0.4)      
text(res[sig,]$log2FoldChange,-log(res[sig,]$padj,10),labels=rownames(res[sig,]),pos=3,cex=0.6)
points(res[sig,"log2FoldChange"],-log(res[sig,"padj"],10),pch=19,cex=0.4,col="red")
abline(h=-log10(p_threshold),lty=3)
abline(v=-lfc_threshold,lty=3)
abline(v=lfc_threshold,lty=3) 

#Heatmap for signifcant miRNAs to show expression across replicates
heatmap.2(vstcounts[sig,],trace="none",col=hmcol,main="Significant miRs",cexRow=0.5,cexCol=0.5,ColSideColors=cond_colours, margins=c(9,7),Colv=FALSE,dendrogram="row")

