
# Change working directory to the Introduction page
setwd("/Users/tm582/Desktop/part2_project/R Introduction")

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

