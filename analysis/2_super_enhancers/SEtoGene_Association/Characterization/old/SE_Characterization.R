#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler", version = "3.8")

library(ggplot2)
library(clusterProfiler)
library(stringr)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers")

SE.Gene <- read.table('./SEtoGene_Association/Cvv_SE_NearestTSS.bed',sep='\t',stringsAsFactors = F)

associatedSEs <- unique(SE.Gene$V4)
length(associatedSEs)

hist(SE.Gene$V9,breaks=20,ylim=c(0,1000),main='',xlab='Distance from Nearest TSS to SE',col='grey80')

sum(SE.Gene$V9 == 0)/length(SE.Gene$V9)
sum(SE.Gene$V9 > 0)/length(SE.Gene$V9)

slices <- c(sum(SE.Gene$V9 == 0),sum(SE.Gene$V9 > 0))
pie(slices,labels=c('Overlapping (85.9%)','Adjacent (14.1%)'),main='SE-Associated Gene Locations\n Relative to SE')

hist(SE.Gene[which(SE.Gene$V9 > 0),9],breaks=40,ylim=c(0,40),xlim=c(0,250000),main='Non-Overlapping SE-Associated Genes',xlab='Distance from Nearest TSS to SE',col='grey80')
boxplot(SE.Gene[which(SE.Gene$V9 > 0),9],xlab='Non-Overlapping SE-Associated Genes',ylim=c(0,250000),ylab='Distance from Nearest TSS to SE',col='grey80')


Venom.SE <- read.table('./SEtoGene_Association/Cvv_VGene_NearestSE.bed',sep='\t',stringsAsFactors = F)
hist(Venom.SE$V9,breaks=50)
