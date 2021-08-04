
library(tidyverse)
library(ggsignif)
library(ggbeeswarm)
library(patchwork)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/promoter_chipseq")

ovrlp <- read_tsv('PromoterH3K4me3_Overlap_12.08.20.bed',col_names = F)

all.prom <- read_tsv('input_files/CvvGeneCoords_v2.bed',col_names = F)

no.ovrlp <- all.prom %>% 
  mutate(txid = str_remove_all(X4,'ID=')) %>% 
  filter(str_detect(txid,'trnascan',negate = T)) %>% 
  filter(!(txid %in% ovrlp$X4))

length(ovrlp$X1)
length(ovrlp$X1) / length(all.prom$X1)

exp <- read_csv('../../data/rnaseq/norm_counts/cvv_VenomRNAseq_NormCounts_updatedGFF_FullVenomAnnot_11.19.19.csv') %>% 
  mutate(avg1DPE = rowMeans(.[,2:4])) %>% 
  mutate(geneID = str_split_fixed(X1,'_',2)[,1]) %>% 
  mutate(txID = str_split_fixed(X1,'_',2)[,2]) %>% 
  mutate(ChIPsignal = ifelse(txID %in% ovrlp$X4,'TRUE','FALSE')) %>% 
  mutate(ChIPsignal = factor(ChIPsignal,levels = c('TRUE','FALSE')))

venPriorityGenes <- read_tsv('../../data/venom_annotations/___PriorityVenomGenes_Sept2020/PriorityVenomGeneList_withNVPs_02.25.21.txt',col_names = F) %>% 
  select(txid = 6)


allGenes <- ggplot(exp,aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5 ) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  ggtitle('d) All Genes') +
  scale_y_continuous(limits=c(0,max(log10(exp$avg1DPE+1)*1.2))) +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

venom.exp <- exp %>% 
  filter(exp$txID %in% venPriorityGenes$txid)

allVen <- ggplot(venom.exp,aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5 ) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('e) All Venom Genes') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))


svmp <- ggplot(subset(venom.exp,str_detect(txID,'SVMP')),aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('f) SVMPs') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

svsp <- ggplot(subset(venom.exp,str_detect(txID,'SVSP')),aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('g) SVSPs') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))

svsp

svsp.2removed <- ggplot(venom.exp %>% filter(str_detect(txID,'SVSP') & txID != 'crovir-transcript-SVSP_2'),aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('g) SVSPs') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))

svsp.2removed

pla2 <- ggplot(subset(venom.exp,str_detect(txID,'PLA2') & txID != 'crovir-transcript-PLA2_gIIE'),aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = .5) +
  scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  ggtitle('h) PLA2s') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

other <- ggplot(subset(venom.exp,str_detect(txID,'SVSP|SVSP|PLA2|128',negate = T)),aes(x=ChIPsignal,y=log10(avg1DPE+1),fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  xlab('H3K4me3 Peak Near Promoter') +
  ggtitle('i) Other VGs') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))

#(allGenes + allVen) / (svmp | svsp | pla2 | other) 

#(allGenes | allVen | svmp | svsp | pla2 | other) + plot_annotation(tag_levels = 'a' )


(allGenes | allVen) / (svmp | svsp ) / (pla2 | other) 

# ### Overlap with promoters near ATACseq peak only in VG3
# 
# atac.ovrlp <- read_tsv('./z_old/promoters_H3K4me3_overlap_genomeWide_extended_1000bp_VG3_not_in_VGU.bed',col_names = F) %>% 
#   mutate(X1 = str_remove_all(X1,'ID='))
# 
# exp.atac <- exp %>% 
#   mutate(atac_overlap = ifelse(txID %in% atac.ovrlp$X1,T,F)) %>% 
#   mutate(both_overlap = ifelse(atac_overlap == T & ChIPsignal == T,T,NA)) %>% 
#   mutate(both_overlap = ifelse(atac_overlap == F & ChIPsignal== F,F,both_overlap))
# 
# ggplot(subset(exp.atac,!is.na(both_overlap)),aes(x=both_overlap,y=log10(avg1DPE+1),color=both_overlap)) +
#   geom_jitter(alpha=0.25,size=0.3) +
#   geom_boxplot(show.legend = F,color='black',fill='NA') +
#   # geom_quasirandom(alpha=0.4) +
#   geom_signif(comparisons = list(c("FALSE", "TRUE")),
#               map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5 ) +
#   scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
#   ggtitle('All Genes') +
#   xlab('H3K4me3 Peak Near Promoter & Unique ATAC-seq peak at VG3') +
#   theme_bw()
