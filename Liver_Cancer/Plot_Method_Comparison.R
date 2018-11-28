##DESeq Plots

#barplot - Pre/post Normalised counts: Using mergedds 
quartz()
par(mfrow=c(1,2))
barplot(colSums(counts(mergedds, normalized=FALSE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Pre Normalised Counts')
legend("topright",levels(conds),cex=0.5,fill=cond_colours[levels(conds)])

barplot(colSums(counts(mergedds, normalized=TRUE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Post Normalised Counts')
legend("topright",levels(conds),cex=0.5,fill=cond_colours[levels(conds)])

#PCA Plot: Using t(assay(mergevsd))
mergepca=prcomp(t(assay(mergevsd)))
quartz()
plot(mergepca$x, main='PCA Variance Stabilised for Merged Data', pch=c(21,24)[batch], col='black', bg=cond_colours[conds], cex=1)
text(mergepca$x, as.character(conds), pos=1, cex=0.4)

#Volcano Plots & Heatmaps: Using res and mergevstMat
p_threshold=0.05
lfc_threshold=0.7
max_display_hits=100

condition_list=levels(conds)

for (i in 2:(length(condition_list))){
  
  normal=condition_list[1]
  compare=condition_list[i]
  
  print(paste('(GTEX/TCGA Merged) Comparing', normal, 'vs.', compare))
  
  res=results(mergedds, contrast = c('Sample_Type', normal, compare))
  res=res[order(res$padj),]
  
  hits=(rownames(res[((res$padj<=p_threshold) & (abs(res$log2FoldChange)>=lfc_threshold) & (!is.na(res$padj))),]))
  
  print(paste("Found: ",length(hits)," Statistical hits"))
  
  if (length(hits) != 0){
    if (length(hits)>max_display_hits){
      hits=hits[1:max_display_hits]
    }
    
    quartz()
    plot(res$log2FoldChange,-log10(res$padj), ylab='-log10(Adjusted P-value)', xlab='log2 FoldChange', main=paste('Volcano Plot:', normal, 'vs.', compare,'(GTEX/TCGA Merged)'), pch=19, cex=0.2)
    points(res[hits,'log2FoldChange'],-log(res[hits,'padj'],10),pch=21,cex=0.4,col='turquoise1')
    text(res[hits,]$log2FoldChange,-log10(res[hits,]$padj), labels = gene_names[hits,]$GeneName, pos=3, cex = 0.4)
    abline(h=-log10(0.05), lty=3, col='red')
    abline(v=-1, lty=3, col='red')
    abline(v=1, lty=3, col='red')
    
    quartz()
    heatmap.2(mergevstMat[hits,],trace='none',col=hmcol,Colv=FALSE,dendrogram='row',main=paste("Heatmap Sig. Hits for", compare),labCol = conds, cexCol = 0.5) 
  } else {
    print(paste("No Results for ",compare))
  }
  
}

##EdgeR Plots

#boxplot - Pre/Post Normalised: use dge object
quartz()
par(mfrow=c(1,2))
boxplot(log2(dge$counts+1), las=2, col=cond_colours[conds], main="Pre Normalised Data", ylab="Log-cpm")

lcpm=cpm(dge, log=T)

boxplot(lcpm, las=2, col=cond_colours[conds], main="Normalised Data",ylab="Log-cpm")

#PCA Plots: using voom object
voompca=prcomp(t(v$E))
quartz()
plot(voompca$x, pch=c(21,24)[batch], bg=cond_colours[conds], cex=1)
text(voompca$x,as.character(coldata$Batch),cex=0.3,pos=1)
text(voompca$x,as.character(colnames(v$E)),cex=0.3,pos = 3)

#Volcano Plots
quartz()
volcanoplot(efit, coef=1, names=as.character(gene_names[rownames(efit),1]), highlight = 100, main='Primary Tumour vs. Normal')
abline(h=-log10(0.05), lty=3, col='red')
abline(v=-1, lty=3, col='red')
abline(v=1, lty=3, col='red')

quartz()
volcanoplot(efit, coef=2, names=as.character(gene_names[rownames(efit),1]), highlight = 100, main='Recurrent Tumour vs. Normal') 
abline(h=-log10(0.05), lty=3, col='red')
abline(v=-1, lty=3, col='red')
abline(v=1, lty=3, col='red')

#Butterfly plot (extra)
quartz()
plot(topTable(efit,coef=1,number=10000000000)$logFC,topTable(efit,coef=2,number=100000000000)$logFC)
abline(h=0)
abline(v=0)
abline(a=0,b=1, lty=3, col='red')


