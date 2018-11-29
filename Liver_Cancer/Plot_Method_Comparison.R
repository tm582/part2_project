##DESeq Plots

#barplot - Pre/post Normalised counts: Using mergedds 
quartz()
barplot(colSums(counts(mergedds, normalized=TRUE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Post Normalised Counts - Regular DESeq')
legend("topright",levels(conds),cex=0.5,fill=cond_colours[levels(conds)])

#PCA Plot: Using t(assay(mergevsd))
quartz()
mergepca=prcomp(t(assay(mergevsd)))
plot(mergepca$x, main='PCA Variance Stabilised - Regular DESeq', pch=c(21,24)[batch], col='black', bg=cond_colours[conds], cex=1)
text(mergepca$x, as.character(conds), pos=1, cex=0.4)

#Rtsne Plot
quartz()
ggplot(data = tsne.df, aes(tsne.1,  tsne.2, colour=Type, shape=Type))+
  geom_point(size=4)+
  scale_color_manual(values=unique(cond_colours[tsne.df$Type]))

#Volcano Plots & Heatmaps: Using res and mergevstMat
p_threshold=0.05
lfc_threshold=0.7
max_display_hits=100

condition_list=levels(conds)
i=2
for (i in 2:(length(condition_list))){
  
  normal=condition_list[1]
  compare=condition_list[i]
  
  print(paste('(GTEX/TCGA Merged (DESeq)) Comparing', normal, 'vs.', compare))
  
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
    heatmap.2(mergevstMat[hits,],trace='none',col=hmcol,Colv=FALSE,dendrogram='row',main=paste("(DESeq) Heatmap Sig. Hits for", compare),labCol = conds, cexCol = 0.5) 
  } else {
    print(paste("No Results for ",compare))
  }
  
}




##EdgeR Plots

#boxplot - Pre/Post Normalised: use dge object
quartz()
boxplot(lcpm, las=2, col=cond_colours[conds], main="Post-Normalised Data - EdgeR",ylab="Log-cpm", cex.axis=0.5)

#PCA Plots: using voom object
voompca=prcomp(t(v$E))
quartz()
plot(voompca$x, pch=c(21,24)[batch], bg=cond_colours[conds], cex=1, main='PCA Plot - EdgeR')
text(voompca$x,as.character(coldata$Batch),cex=0.3,pos=1)
text(voompca$x,as.character(colnames(v$E)),cex=0.3,pos = 3)

#Volcano Plots
quartz()
volcanoplot(efit, coef=1, names=as.character(gene_names[rownames(efit),1]), highlight = 100, main='Primary Tumour vs. Normal (EdgeR)')
abline(h=-log10(0.05), lty=3, col='red')
abline(v=-1, lty=3, col='red')
abline(v=1, lty=3, col='red')

quartz()
volcanoplot(efit, coef=2, names=as.character(gene_names[rownames(efit),1]), highlight = 100, main='Recurrent Tumour vs. Normal (EdgeR)') 
abline(h=-log10(0.05), lty=3, col='red')
abline(v=-1, lty=3, col='red')
abline(v=1, lty=3, col='red')

#Butterfly plot (extra)
quartz()
plot(topTable(efit,coef=1,number=10000000000)$logFC,topTable(efit,coef=2,number=100000000000)$logFC)
abline(h=0)
abline(v=0)
abline(a=0,b=1, lty=3, col='red')




##Rlog Transformation

#Normalised Counts
quartz()
boxplot(assay(rlt), col=cond_colours[conds], las=2, cex.axis=0.5, main='rlog Normalised Counts', pch=16, cex=0.3)

#Sample Distance Matrix
quartz()
heatmap.2(sampleDistMatrix, cexRow = 0.5, cexCol = 0.5, col = hmcol, trace = 'none', main = 'Sample-to-Sample Distances Using Rlog Values')

#PCA Plot
quartz()
plot(rltpca$x, pch=c(21,24)[batch], bg=cond_colours[conds], main='PCA Plot: rlog values')
text(rltpca$x,as.character(coldata$Batch),cex=0.3,pos=1)
text(rltpca$x,as.character(coldata$Sample_Type),cex=0.3,pos = 3)

#Most Significant Gene(s) Plot
quartz()
plotCounts(dds, gene=topGene, intgroup=c("Sample_Type"))

quartz()
heatmap.2(assay(rlt)[rownames(res[1:100,]),],col=hmcol,trace="none",cexRow = 0.5,cexCol = 0.5, main = 'Heatmap of top 100 most significant genes')




##SVA-DESeq Plots

#Normalised Counts
quartz()
barplot(colSums(counts(ddssva, normalized=TRUE)), col=cond_colours[conds], las=2,cex.names=0.5,main='Post Normalised Counts - SVA')
legend("topleft",levels(conds),cex=0.5,fill=cond_colours[levels(conds)])

#PCA Plot
quartz()
plot(svapca$x, main='PCA Variance Stabilised for SVA', pch=c(21,24)[batch], col='black', bg=cond_colours[conds], cex=1)
text(svapca$x, as.character(conds), pos=1, cex=0.4)

#Rtsne
quartz()
ggplot(data = svatsne.df, aes(svatsne.1,  svatsne.2, colour=Type, shape=Batch))+
  geom_point(size=4)+
  scale_color_manual(values=unique(cond_colours[svatsne.df$Type]))

#Volcano Plots & Heatmaps

for (i in 2:(length(condition_list))){
  
  normal=condition_list[1]
  compare=condition_list[i]
  
  print(paste('(GTEX/TCGA (SVA)) Comparing', normal, 'vs.', compare))
  
  svares=results(ddssva, contrast = c('Sample_Type', normal, compare))
  svares=svares[order(svares$padj),]
  
  svahits=(rownames(svares[((svares$padj<=p_threshold) & (abs(svares$log2FoldChange)>=lfc_threshold) & (!is.na(svares$padj))),]))
  
  print(paste("Found: ",length(svahits)," Statistical hits"))
  
  if (length(svahits) != 0){
    if (length(svahits)>max_display_hits){
      svahits=svahits[1:max_display_hits]
    }
    
    quartz()
    plot(svares$log2FoldChange,-log10(svares$padj), ylab='-log10(Adjusted P-value)', xlab='log2 FoldChange', main=paste('Volcano Plot:', normal, 'vs.', compare,'(GTEX/TCGA (SVA)'), pch=19, cex=0.2)
    points(svares[hits,'log2FoldChange'],-log(svares[svahits,'padj'],10),pch=21,cex=0.4,col='turquoise1')
    text(svares[hits,]$log2FoldChange,-log10(svares[svahits,]$padj), labels = gene_names[svahits,]$GeneName, pos=3, cex = 0.4)
    abline(h=-log10(0.05), lty=3, col='red')
    abline(v=-1, lty=3, col='red')
    abline(v=1, lty=3, col='red')
    
    quartz()
    heatmap.2(mergevstMat[svahits,],trace='none',col=hmcol,Colv=FALSE,dendrogram='row',main=paste("(SVA) Heatmap Sig. Hits for", compare),labCol = conds, cexCol = 0.5) 
  } else {
    print(paste("No Results for ",compare))
  }
  
}



