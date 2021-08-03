
library(stringr)
library(matrixStats)
library(pheatmap)
library(viridis)
library(grid)

draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

####################################################################################################

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/")

SEgene.pairs <- read.table('SEtoGene_Association/Cvv_SEassocPairs_Genes.bed',sep='\t',stringsAsFactors = F)
SEgene.pairs[,4:5] <- str_split_fixed(SEgene.pairs$V4,pattern = '@',2)
SEgene.pairs$V5 <- str_split_fixed(sub(".*Crovir_Transcript_ID=", "", SEgene.pairs$V5),';',2)[,1]

IDconvert <- read.csv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/Gene_expression/_newRNAseq_STAR_03.13.19/tx2gene_v3.csv',stringsAsFactors = F,header = F)

SEgene.pairs.ID <- merge(SEgene.pairs,IDconvert,by.x='V5',by.y='V1')

SEgene.pairs.ID <- as.data.frame(SEgene.pairs.ID$V2.y)

tf.list <- read.table('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/Transcription_Factors/TF_lists/Cvv_TF.DNABinding.txt')[4]


TF.SEgene.pairs <- as.data.frame(SEgene.pairs.ID[which(SEgene.pairs.ID$`SEgene.pairs.ID$V2.y` %in% tf.list$V4),])

write.table(TF.SEgene.pairs,'SEassocTFsList_03.20.19.txt',sep='\t',quote = F,row.names = F,col.names = F)

TFSEs <- 68
TFnotSEs <- length(tf.list$V4)-68

GeneSEs <- 946
GeneNotSEs <-  18538 - 946

TFSE_FisherMatrix <- matrix(c(TFSEs,GeneSEs,TFnotSEs,GeneNotSEs),nrow=2)
fisher.test(TFSE_FisherMatrix)

normcounts <- read.csv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/Gene_expression/z_old/Salmon_preSTAR_10.13.18/Cvv_totRNA_NormCounts_10.29.18.csv',header = T,row.names = 1)

tf.VGnormcounts <- as.matrix(normcounts[which(row.names(normcounts) %in% TF.SEgene.pairs[,1]),])

tf.avg.VGnormcounts <- as.data.frame(tf.VGnormcounts[,1:2])
tf.avg.VGnormcounts$Panc.Avg <- rowMeans2(tf.VGnormcounts,cols=3:5)
tf.avg.VGnormcounts$SI.Avg <- rowMeans2(tf.VGnormcounts,cols=6:8)
tf.avg.VGnormcounts$Skin.Avg <- rowMeans2(tf.VGnormcounts,cols=9:11)
tf.avg.VGnormcounts$Stomach.Avg <- rowMeans2(tf.VGnormcounts,cols=12:14)
tf.avg.VGnormcounts$UnextVG.Avg <- rowMeans2(tf.VGnormcounts,cols=15:17)
tf.avg.VGnormcounts$VG1DPE.Avg <- rowMeans2(tf.VGnormcounts,cols=18:20)
tf.avg.VGnormcounts$VG3DPE.Avg <- rowMeans2(tf.VGnormcounts,cols=21:23)

tf.avg.VGnormcounts <- tf.avg.VGnormcounts[,c(7,8,9,1,2,3,4,5,6)]
tf.avg.VGnormcounts <- tf.avg.VGnormcounts[order(tf.avg.VGnormcounts$VG1DPE.Avg,decreasing = T),]

row.names(tf.avg.VGnormcounts) <- str_split_fixed(row.names(tf.avg.VGnormcounts),'_',2)[,1]

allGenes.Unext.vs.1DPE <- read.csv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/Gene_expression/z_old/Salmon_preSTAR_10.13.18/AllDEgenes_Unext.vs.1DPE_12.05.18.csv',header = T,row.names = 1)
TF.Unext.vs.1DPE <- allGenes.Unext.vs.1DPE[which(row.names(allGenes.Unext.vs.1DPE) %in% row.names(tf.avg.VGnormcounts)),]

TFs.Ext.vs.NonVenom <- read.csv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation//Gene_expression/z_old/Salmon_preSTAR_10.13.18/TF_ExtVGvsNonVG_UpregulatedPWresults_12.17.18.csv',header = T,row.names = 1)
rownames(TFs.Ext.vs.NonVenom) <- str_split_fixed(row.names(TFs.Ext.vs.NonVenom),'_',2)[,1]

TFs.Ext.vs.NonVenom <- tf.avg.VGnormcounts[which(row.names(tf.avg.VGnormcounts) %in% row.names(TFs.Ext.vs.NonVenom)),]

#write.table(row.names(TFs.Ext.vs.NonVenom),'UpregSEassocTFs_extVGvsNonVG.txt',quote=F)

pheatmap(tf.avg.VGnormcounts,scale='row',cluster_cols = F,cluster_rows = T,col=magma(50),show_rownames = T,border_color = NA,cellwidth = 15,cellheight = 15)
pheatmap(log2(TFs.Ext.vs.NonVenom+0.1),scale='none',cluster_cols = F,cluster_rows = T,col=magma(50),show_rownames = T,border_color = NA,cellwidth = 15,cellheight = 15)


## NFKB1

nfkb1 <- normcounts[which(row.names(normcounts)=='NFKB1_crovir-transcript-4379'),]
nfkb1 <- nfkb1[,-c(9:11)]
pheatmap(nfkb1,scale='row',cluster_cols = F,cluster_rows = F,col=magma(50),show_rownames = T,border_color = NA,cellwidth = 15,cellheight = 15)
