
library(Sushi)
library(viridis)
library(tidyr)

#Load venom gene annotations
all.venomGenes <- read.table('/Users/perryb/Desktop/Desktop_dir/CroVir_IGB/__Cvv_IGBtracks_Updated-01.08.19/VenomGene_regions_full.bed',sep='\t',stringsAsFactors = F)

### SVMP

# Subset out SVMPs
svmp.genes <- all.venomGenes[c(1:11),]
svmp.genes2 <- read.table('../../../MultiTrackPlotting_Sushi/misc_data/VenomGene_regions-mi1.bed',sep='\t',stringsAsFactors = F)


#Load TADs
TADs <- read.table('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/HiC/Analyses/TAD_Identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hier1_01.18.19.bed',sep='\t',stringsAsFactors = F)

#Load and format HiC matrix
svmp.hic.long <- read.table('./misc_data/mi1_intra_10kb_contact_matrix_01.16.19.txt',sep='\t',stringsAsFactors = F)
svmp.hic.long$V3 <- log10(svmp.hic.long$V3+.1)
colnames(svmp.hic.long) <- c('b1','b2','count')

svmp.hic.wide <- spread(svmp.hic.long,key = b2,value = count, fill=0)
row.names(svmp.hic.wide) <- svmp.hic.wide$b1
svmp.hic.wide <- svmp.hic.wide[,-1]
svmp.hic.wide[lower.tri(svmp.hic.wide)] = t(svmp.hic.wide)[lower.tri(svmp.hic.wide)]


chrom_svmp = 'scaffold-mi1'
chromstart_svmp = min(svmp.genes$V2-105000)
chromend_svmp =  max(svmp.genes$V3+105000)

par(mfrow=c(3,1),mar=c(1,2,1,1))

plotHic(svmp.hic.wide,chrom_svmp,chromstart_svmp,chromend_svmp,max_y=20,zrange = c(0,1.5),palette = viridis)
labelgenome(chrom_svmp,chromstart_svmp,chromend_svmp,n=4,scale="Mb")
mtext("HiC",side=2,line=0.5,font=2)

plotBed(TADs,chrom_svmp,chromstart_svmp,chromend_svmp,row = 'auto',color='seagreen')
labelgenome(chrom_svsp,chromstart_svsp,chromend_svsp,n=4,scale="Mb")
mtext("TADs",side=2,line=0.5,font=2)

plotGenes(svmp.genes2,chrom_svmp,chromstart_svmp,chromend_svmp,row = 'supplied',bentline = F,plotgenetype = 'box')
mtext("Gene",side=2,line=0.5,font=2)
#plotBed(svmp.genes,chrom_svmp,chromstart_svmp,chromend_svmp,row=1)



### SVSP

# Subset out svsps
svsp.genes <- read.table('../../../MultiTrackPlotting_Sushi/misc_data/VenomGene_regions-mi2.bed',sep='\t',stringsAsFactors = F)

#Load and format HiC matrix
svsp.hic.long <- read.table('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/HiC/HiC_Matrices/10kb/mi2_intra_10kb_contact_matrix_01.16.19.txt',sep='\t',stringsAsFactors = F)
svsp.hic.long$V3 <- log10(svsp.hic.long$V3+.1)
colnames(svsp.hic.long) <- c('b1','b2','count')

svsp.hic.wide <- spread(svsp.hic.long,key = b2,value = count, fill=0)
row.names(svsp.hic.wide) <- svsp.hic.wide$b1
svsp.hic.wide <- svsp.hic.wide[,-1]
svsp.hic.wide[lower.tri(svsp.hic.wide)] = t(svsp.hic.wide)[lower.tri(svsp.hic.wide)]


chrom_svsp = 'scaffold-mi2'
chromstart_svsp = min(svsp.genes$V2-105000)
chromend_svsp =  max(svsp.genes$V3+105000)

par(mfrow=c(3,1),mar=c(1,2,1,1))

plotHic(svsp.hic.wide,chrom_svsp,chromstart_svsp,chromend_svsp,max_y=20,zrange = c(0,1.5),palette = viridis)
labelgenome(chrom_svsp,chromstart_svsp,chromend_svsp,n=4,scale="Mb")
mtext("HiC",side=2,line=0.5,font=2)

plotBed(TADs,chrom_svsp,chromstart_svsp,chromend_svsp,row = 'auto',color='seagreen')
mtext("TADs",side=2,line=0.5,font=2)
labelgenome(chrom_svsp,chromstart_svsp,chromend_svsp,n=4,scale="Mb")

plotGenes(svsp.genes,chrom_svsp,chromstart_svsp,chromend_svsp,row = 'supplied',bentline = F,plotgenetype = 'box')
mtext("Gene",side=2,line=0.5,font=2)



### pla2

# Subset out pla2s
pla2.genes <- read.table('../../../MultiTrackPlotting_Sushi/misc_data/VenomGene_regions-mi7ExonsOnly.Sushi.bed',sep='\t',stringsAsFactors = F)

#Load and format HiC matrix
pla2.hic.long <- read.table('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/HiC/HiC_Matrices/10kb/mi7_intra_10kb_contact_matrix_01.16.19.txt',sep='\t',stringsAsFactors = F)
pla2.hic.long$V3 <- log10(pla2.hic.long$V3+.1)
colnames(pla2.hic.long) <- c('b1','b2','count')

pla2.hic.wide <- spread(pla2.hic.long,key = b2,value = count, fill=0)
row.names(pla2.hic.wide) <- pla2.hic.wide$b1
pla2.hic.wide <- pla2.hic.wide[,-1]
pla2.hic.wide[lower.tri(pla2.hic.wide)] = t(pla2.hic.wide)[lower.tri(pla2.hic.wide)]


chrom_pla2 = 'scaffold-mi7'
chromstart_pla2 = min(pla2.genes$V2-300000)
chromend_pla2 =  max(pla2.genes$V3+300000)

par(mfrow=c(3,1),mar=c(1,2,1,1))

plotHic(pla2.hic.wide,chrom_pla2,chromstart_pla2,chromend_pla2,max_y=20,zrange = c(0,1.5),palette = viridis)
labelgenome(chrom_pla2,chromstart_pla2,chromend_pla2,n=4,scale="Mb")
mtext("HiC",side=2,line=0.5,font=2)

plotBed(TADs,chrom_pla2,chromstart_pla2,chromend_pla2,row = 'auto',color='seagreen')
mtext("TADs",side=2,line=0.5,font=2)
labelgenome(chrom_pla2,chromstart_pla2,chromend_pla2,n=4,scale="Mb")

plotGenes(pla2.genes,chrom_pla2,chromstart_pla2,chromend_pla2,row = 'supplied',bentline = F,plotgenetype = 'box')
mtext("Gene",side=2,line=0.5,font=2)





### PLOT ALL
par(mfrow=c(3,3),mar=c(1,2,1,1))

plotHic(svmp.hic.wide,chrom_svmp,chromstart_svmp,chromend_svmp,max_y=20,zrange = c(0,1.5),palette = viridis)
labelgenome(chrom_svmp,chromstart_svmp,chromend_svmp,n=4,scale="Mb")
mtext("HiC",side=2,line=0.5,font=2)

plotHic(svsp.hic.wide,chrom_svsp,chromstart_svsp,chromend_svsp,max_y=20,zrange = c(0,1.5),palette = viridis)
labelgenome(chrom_svsp,chromstart_svsp,chromend_svsp,n=4,scale="Mb")

plotHic(pla2.hic.wide,chrom_pla2,chromstart_pla2,chromend_pla2,max_y=20,zrange = c(0,1.5),palette = viridis)
labelgenome(chrom_pla2,chromstart_pla2,chromend_pla2,n=4,scale="Mb")

plotBed(TADs,chrom_svmp,chromstart_svmp,chromend_svmp,row = 'auto',color='seagreen')
labelgenome(chrom_svsp,chromstart_svsp,chromend_svsp,n=4,scale="Mb")
mtext("TADs",side=2,line=0.5,font=2)

plotBed(TADs,chrom_svsp,chromstart_svsp,chromend_svsp,row = 'auto',color='seagreen')
labelgenome(chrom_svsp,chromstart_svsp,chromend_svsp,n=4,scale="Mb")

plotBed(TADs,chrom_pla2,chromstart_pla2,chromend_pla2,row = 'auto',color='seagreen')
labelgenome(chrom_pla2,chromstart_pla2,chromend_pla2,n=4,scale="Mb")

plotGenes(svmp.genes2,chrom_svmp,chromstart_svmp,chromend_svmp,row = 'supplied',bentline = F,plotgenetype = 'box')
mtext("Gene",side=2,line=0.5,font=2)

plotGenes(svsp.genes,chrom_svsp,chromstart_svsp,chromend_svsp,row = 'supplied',bentline = F,plotgenetype = 'box')

plotGenes(pla2.genes,chrom_pla2,chromstart_pla2,chromend_pla2,row = 'supplied',bentline = F,plotgenetype = 'box')











