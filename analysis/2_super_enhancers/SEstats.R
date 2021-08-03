
#install.packages('ggrepel')
#install.packages('cowplot')

library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers")

SEs <- read.table('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',sep='\t',stringsAsFactors = F)

chromLengths <- read.table('/Users/perryb/Desktop/Desktop_dir/Cvv_ChrLengths.txt',sep='\t',stringsAsFactors = F)
chromLengths$V1 <- gsub("scaffold-","",chromLengths$V1) 

SEs$length <- abs(SEs$V3-SEs$V2)


#hist(SEs$length,breaks=20,col='grey50',main='',xlab='Super-Enhancer Length (bp)')

length_hist <- ggplot(SEs,aes(x=length)) + 
  geom_histogram(color='black',binwidth=5000) + 
  ylab(label='Count') + xlab(label='Super-Enhancer Length (bp)') +
  scale_y_continuous(expand = c(0,0),limits=c(0,75)) +
  geom_vline(aes(xintercept=mean(length)),linetype='dashed',color='red',size=1)

length_hist


mean(SEs$length)
median(SEs$length)



chrom_SECounts <- as.data.frame(table(unlist(SEs$V1)),stringsAsFactors = F)
chrom_SECounts <- rbind(chrom_SECounts,c('scaffold-UNAS',22))
chrom_SECounts$Freq <- as.numeric(chrom_SECounts$Freq)
chrom_SECounts <- chrom_SECounts[c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,9,31,32),]
chrom_SECounts$Var1 <- gsub("scaffold-","",chrom_SECounts$Var1) 

chromCount <- ggplot(chrom_SECounts,aes(x=Var1,y=Freq)) + 
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0),limits=c(0,110)) +
  ylab(label='# of SEs') +
  xlab(label='') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

chromCount  

length_by_count <- merge(chromLengths,chrom_SECounts,by.x=1,by.y=1,all.x=T)
#length_by_count$V1 <- as.factor(gsub("scaffold-", "", length_by_count$V1))

lengthCount_Scatter <- ggplot(length_by_count,aes(x=V2,y=Freq,label=V1)) + 
  geom_point() +
  geom_text_repel() +
  ylab(label='Number of SEs') + xlab(label='Length of Chromosome')
  
lengthCount_Scatter

SEchromLength.lm <- lm(Freq ~ V2,data=length_by_count)
summary(SEchromLength.lm)

lengthCount_Scatter_regLine <- lengthCount_Scatter + geom_smooth(method='lm',formula = y~x,se=T,lty=2,lwd=.5,fill='grey90')  +
  labs(title = paste("Adj R2 = ",signif(summary(SEchromLength.lm)$adj.r.squared, 5),
                     "Intercept =",signif(SEchromLength.lm$coef[[1]],5 ),
                     " Slope =",signif(SEchromLength.lm$coef[[2]], 5),
                     " P =",signif(summary(SEchromLength.lm)$coef[2,4], 5)))


plot_grid(length_hist,lengthCount_Scatter_regLine,ncol=2,align='v')

#plot_grid(length_hist,lengthCount_Scatter,chromCount,se.plot,ncol=2,labels=c('A','B','C','D'),align='v')



#### Hist of genes per SE

SEgene.pairs <- read.table('SEtoGene_Association/Cvv_SEassocPairs_Genes.bed',sep='\t',stringsAsFactors = F)
SEgene.pairs[,4:5] <- str_split_fixed(SEgene.pairs$V4,pattern = '@',2)


genePerSE <- as.data.frame(str_split_fixed(SEgene.pairs$V4,pattern = '@',2))

genePerSE <- as.data.frame(table(genePerSE$V1))

genePerSE_hist <- ggplot(genePerSE,aes(x=Freq)) + 
  geom_histogram(color='black') + 
  ylab(label='Count') + xlab(label='# Genes per SE') +
  scale_y_continuous(expand = c(0,0),limits=c(0,320))

genePerSE_hist  

min5geneSEs <- genePerSE[which(genePerSE$Freq >= 5),]
min5geneSEs <- SEgene.pairs[which(SEgene.pairs$V4 %in% min5geneSEs$Var1),]

max5geneSEs <- genePerSE[which(genePerSE$Freq < 5),]
max5geneSEs <- SEgene.pairs[which(SEgene.pairs$V4 %in% max5geneSEs$Var1),]

#write.table(min5geneSEs,'SEassocGenes_min5geneSEs.txt',sep='\t',quote=F,col.names=F,row.names = F)
#write.table(max5geneSEs,'SEassocGenes_max5geneSEs.txt',sep='\t',quote=F,col.names=F,row.names = F)

