#source("https://bioconductor.org/biocLite.R")
#biocLite("tximport")
#install.packages('pheatmap')
#install.packages('viridis')

#install.packages('matrixStats')

library(pheatmap)
library(viridis)
library(stringr)
library(matrixStats)


setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/SEtoGene_Association/Characterization")

normcounts <- read.csv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/Transcription_Factors/Analyses/Gene_expression/Salmon_preSTAR_10.13.18/Cvv_totRNA_NormCounts_10.29.18.csv',header=T,row.names=1)
rownames(normcounts) <- str_split_fixed(row.names(normcounts),'_',2)[,2]

SEgenes <- read.table('SEassocGenes_TranscIDs.txt')

SEgene.normcounts <-as.matrix(normcounts[which(row.names(normcounts) %in% SEgenes$V1),])

SEgene.avg.normcounts <- as.data.frame(SEgene.normcounts[,1:2])
SEgene.avg.normcounts$Panc.Avg <- rowMeans2(SEgene.normcounts,cols=3:5)
SEgene.avg.normcounts$SI.Avg <- rowMeans2(SEgene.normcounts,cols=6:8)
SEgene.avg.normcounts$Skin.Avg <- rowMeans2(SEgene.normcounts,cols=9:11)
SEgene.avg.normcounts$Stomach.Avg <- rowMeans2(SEgene.normcounts,cols=12:14)
SEgene.avg.normcounts$UnextVG.Avg <- rowMeans2(SEgene.normcounts,cols=15:17)
SEgene.avg.normcounts$VG1DPE.Avg <- rowMeans2(SEgene.normcounts,cols=18:20)
SEgene.avg.normcounts$VG3DPE.Avg <- rowMeans2(SEgene.normcounts,cols=21:23)

SEgene.avg.normcounts <- SEgene.avg.normcounts[,c(7,8,9,1,2,3,4,5,6)]
SEgene.avg.normcounts <- SEgene.avg.normcounts[order(SEgene.avg.normcounts$VG1DPE.Avg,decreasing = T),]


pheatmap(SEgene.avg.normcounts,scale='row',cluster_cols = F,cluster_rows = T,col=magma(50),show_rownames = F,border_color = NA,cellwidth = 10)
pheatmap(log2(SEgene.avg.normcounts+1),scale='none',cluster_cols = T,cluster_rows = T,clustering_method = 'mcquitty',col=magma(50),show_rownames = T,fontsize_row = 7,border_color = NA,cellwidth = 10,treeheight_row = 20,treeheight_col = 15,main = 'log2(norm. count  + 1)')

de.genes <- read.csv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/Transcription_Factors/Analyses/Gene_expression/Salmon_preSTAR_10.13.18/AllDEgenes_Unext.vs.1DPE_12.05.18.csv',header=T,row.names=1)
rownames(de.genes) <- str_split_fixed(row.names(de.genes),'_',2)[,2]
de.genes <- de.genes[which(de.genes$IHW_pvalue < 0.05),]


SEgene.normcounts <-as.matrix(normcounts[which(row.names(normcounts) %in% SEgenes$V1),])

NonSEgene.normcounts <- as.matrix(normcounts[which(!(row.names(normcounts) %in% SEgenes$V1)),])


