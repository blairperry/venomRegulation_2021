
library(tidyverse)
library(circlize)

###
### NOTE: NOT YET FINISHED. As of 8/5/21 at 8:40am, this script does not work properly. 
###



cvv.chroms <- read_tsv('./data/misc/Cvv_ChromSizes.txt') %>% 
  mutate(color = ifelse(str_detect(old_id,'ma'),'grey','goldenrod')) %>% 
  mutate(color = ifelse(old_id == 'scaffold-Z','seagreen',color))

cvv.exp <- read_csv('analysis/1_gene_expression/norm_counts/CvvVenomReg_RepRNAseq_wNonVen_VSTNormCounts_08.02.21.csv') %>% 
  mutate(avg1DPE = rowMeans(.[,c(4,8,9)])) %>% 
  select(txid = 1, avg1DPE) %>% 
  mutate(txid = str_split_fixed(txid,'[_]',2)[,2])

cvv.gff <- read_tsv('data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
  select(chr=1,start=4,end=5,desc=9) %>% 
  filter(str_detect(desc,'trnascan',negate = T)) %>% 
  mutate(txid = str_split_fixed(desc,'[;]',4)[,3]) %>% 
  mutate(txid = str_remove(txid,'Crovir_Transcript_ID=')) %>% 
  select(1,2,3,5)

cvv.exp.gff <- cvv.exp %>% 
  left_join(cvv.gff) %>% 
  left_join(cvv.chroms,by=c('chr'='old_id')) %>% 
  select(chr = 6, start=4, end=5, avg1DPE)



###### Run this whole block together each time
circos.par(start.degree=90,gap.degree=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,10))
circos.genomicInitialize(cvv.chroms,plotType = NULL)

circos.track(ylim = c(0, 1), 
             bg.col = cvv.chroms$color, 
             bg.border = 'black', track.height = 0.08,panel.fun = function(x, y) {
               circos.genomicAxis(h = "bottom", direction = "inside",labels.cex = 0.5)}
             )

circos.genomicTrackPlotRegion(cvv.exp.gff, panel.fun = function(region, value, ...) {
  circos.genomicLines(cvv.exp.gff, cvv.exp.gff$avg1DPE)
})

circos.genomicLabels(cvv.chroms, labels.column = 1, side = "inside",niceFacing = T,cex = 1,line_lwd = 0)



# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   chr = CELL_META$sector.index
#   xlim = CELL_META$xlim
#   ylim = CELL_META$ylim
#   circos.rect(xlim[1], 0, xlim[2], 1, col = c('green','white'))
#   circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "black",
#               facing = "inside", niceFacing = TRUE)
# }, track.height = 0.15, bg.border = NA)

circos.clear()

######