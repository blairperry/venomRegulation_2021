

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/SEtoGene_Association/Characterization/GO_analysis")

BiolProc <- read.table('BiolProc/enrichment_results_wg_result1546959400.txt',header = 1,sep='\t',stringsAsFactors = F)
CellComp <- read.table('CellComp/enrichment_results_wg_result1546959304.txt',header = 1,sep='\t',stringsAsFactors = F)
MolFunc <- read.table('MolFun/enrichment_results_wg_result1546959173.txt',header = 1,sep='\t',stringsAsFactors = F)

BiolProc <- BiolProc[which(BiolProc$FDR < 0.05),]
CellComp <- CellComp[which(CellComp$FDR < 0.05),]
MolFunc <- MolFunc[which(MolFunc$FDR < 0.05),]

allGO <- rbind(MolFunc,BiolProc)
allGO <- allGO[order(allGO$R,decreasing = T),]

par(las=2)
par(mar=c(20,4,1,1)) # adjust as needed

barplot(allGO$R,names.arg = allGO$description,ylim=c(0,6),ylab='Ratio of enrichment')

