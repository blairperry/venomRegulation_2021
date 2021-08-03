
data <- read.table('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/SEtoGene_Association/Cvv_SE_NearestTSS.bed',sep='\t',stringsAsFactors = F)
gff <- read.table('/Users/perryb/Downloads/CroVir_rnd1.all.maker.final.homologIDs.gff',sep='\t',stringsAsFactors = F)

data2 <- merge(data,gff,by.x='V8',by.y='V9',all.x=T)

data2 <- data2[,c(2,3,4,5,10,13,14,1,9)]

write.table(data2,'/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/SEtoGene_Association/Cvv_SE_NearestTSS_GeneCoords.bed',sep='\t',quote = F, row.names=F, col.names = F)
