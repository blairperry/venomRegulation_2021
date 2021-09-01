
#library(devtools)
#install_github("wbaopaul/rGMAP")

library(rGMAP)


## 10kb
ma1 <- read.table('input_files/10kb/ma1_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
ma2 <- read.table('input_files/10kb/ma2_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
ma3 <- read.table('input_files/10kb/ma3_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
ma4 <- read.table('input_files/10kb/ma4_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
ma5 <- read.table('input_files/10kb/ma5_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
ma6 <- read.table('input_files/10kb/ma6_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
ma7 <- read.table('input_files/10kb/ma7_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
maZ <- read.table('input_files/10kb/maZ_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi1 <- read.table('input_files/10kb/mi1_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi2 <- read.table('input_files/10kb/mi2_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi3 <- read.table('input_files/10kb/mi3_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi4 <- read.table('input_files/10kb/mi4_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi5 <- read.table('input_files/10kb/mi5_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi6 <- read.table('input_files/10kb/mi6_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi7 <- read.table('input_files/10kb/mi7_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi8 <- read.table('input_files/10kb/mi8_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi9 <- read.table('input_files/10kb/mi9_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)
mi10 <-read.table('input_files/10kb/mi10_intra_10kb_rGMAPinput_01.16.19.txt',sep='\t',col.names=c('n1','n2','count'),stringsAsFactors = F)

res.ma1 = rGMAP(ma1,resl = 10000)
res.ma2 = rGMAP(ma2,resl = 10000)
res.ma3 = rGMAP(ma3,resl = 10000)
res.ma4 = rGMAP(ma4,resl = 10000)
res.ma5 = rGMAP(ma5,resl = 10000)
res.ma6 = rGMAP(ma6,resl = 10000)
res.ma7 = rGMAP(ma7,resl = 10000)
res.maZ = rGMAP(maZ,resl = 10000)
res.mi1 = rGMAP(mi1,resl = 10000)
res.mi2 = rGMAP(mi2,resl = 10000)
res.mi3 = rGMAP(mi3,resl = 10000)
res.mi4 = rGMAP(mi4,resl = 10000)
res.mi5 = rGMAP(mi5,resl = 10000)
res.mi6 = rGMAP(mi6,resl = 10000)
res.mi7 = rGMAP(mi7,resl = 10000)
res.mi8 = rGMAP(mi8,resl = 10000)
res.mi9 = rGMAP(mi9,resl = 10000)
res.mi10 = rGMAP(mi10,resl = 10000)


write.table(res.ma1$hierTads,'TADs/10kb/ma1_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.ma2$hierTads,'TADs/10kb/ma2_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.ma3$hierTads,'TADs/10kb/ma3_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.ma4$hierTads,'TADs/10kb/ma4_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.ma5$hierTads,'TADs/10kb/ma5_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.ma6$hierTads,'TADs/10kb/ma6_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.ma7$hierTads,'TADs/10kb/ma7_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.maZ$hierTads,'TADs/10kb/maZ_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi1$hierTads,'TADs/10kb/mi1_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi2$hierTads,'TADs/10kb/mi2_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi3$hierTads,'TADs/10kb/mi3_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi4$hierTads,'TADs/10kb/mi4_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi5$hierTads,'TADs/10kb/mi5_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi6$hierTads,'TADs/10kb/mi6_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi7$hierTads,'TADs/10kb/mi7_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi8$hierTads,'TADs/10kb/mi8_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi9$hierTads,'TADs/10kb/mi9_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)
write.table(res.mi10$hierTads,'TADs/10kb/mi10_10kb_HierTADS_01.16.19.txt',quote=F,sep='\t',row.names = F)



###
### Misc. TAD Stats
###

all.h1 <- read.table('./TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hier1_01.18.19.bed',sep='\t',stringsAsFactors = F)
all.h2 <- read.table('./TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hier2_01.18.19.bed',sep='\t',stringsAsFactors = F)


# histogram of TAD size

len.h1 <- all.h1$V3-all.h1$V2
len.h2 <- all.h2$V3-all.h2$V2

median(len.h1)/1000
median(len.h2)/1000

hist(len.h1,breaks=50,col='grey',ylim=c(0,200),xlab='Length (bp)',main='TADs',xlim=c(0,3500000))
hist(len.h2,breaks=50,col='grey',ylim=c(0,200),xlab='Length (bp)',main='subTADs',xlim=c(0,3500000))

# TADs per chromosome
chrom_TADCounts.h1 <- as.data.frame(table(unlist(all.h1$V1)),stringsAsFactors = F)
chrom_TADCounts.h2 <- as.data.frame(table(unlist(all.h2$V1)),stringsAsFactors = F)

par(las=3,mar=c(8,3,1,1))
barplot(chrom_TADCounts.h1$Freq,names.arg = chrom_TADCounts.h1$Var1,ylim=c(0,600))
barplot(chrom_TADCounts.h2$Freq,names.arg = chrom_TADCounts.h2$Var1,ylim=c(0,600))
par(las=1,mar=c(8,3,1,1))




