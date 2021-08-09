
library(tidyverse)
library(ggsignif)
library(ggbeeswarm)
library(patchwork)

ovrlp <- read_tsv('./analysis/5_promoters/PromoterH3K4me3_Overlap_08.09.21.bed',col_names = F)

all.prom <- read_tsv('./analysis/5_promoters/CvvGeneCoords_forABC_08.09.21.bed',col_names = F)

no.ovrlp <- all.prom %>% 
  mutate(txid = str_remove_all(X4,'ID=')) %>% 
  filter(str_detect(txid,'trnascan',negate = T)) %>% 
  filter(!(txid %in% ovrlp$X4))

length(ovrlp$X1)
length(ovrlp$X1) / length(all.prom$X1)

exp <- read_csv('./analysis/1_gene_expression/norm_counts/CvvVenomReg_RepRNAseq_wNonVen_VSTNormCounts_08.02.21.csv') %>% 
  mutate(avg1DPE = rowMeans(.[,c(4,8,9)])) %>% 
  mutate(geneID = str_split_fixed(X1,'_',2)[,1]) %>% 
  mutate(txID = str_split_fixed(X1,'_',2)[,2]) %>% 
  mutate(ChIPsignal = ifelse(txID %in% ovrlp$X4,'TRUE','FALSE')) %>% 
  mutate(ChIPsignal = factor(ChIPsignal,levels = c('TRUE','FALSE')))

venPriorityGenes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F) %>% 
  select(txid = 6)


allGenes <- ggplot(exp,aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5 ) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  ggtitle('d) All Genes') +
  # scale_y_continuous(limits=c(0,max(log10(exp$avg1DPE+1)*1.2))) +
  xlab('H3K4me3 Peak Near Promoter') +
  ylab('VST Normalized Gene Expression')+
  theme_bw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

allGenes

venom.exp <- exp %>% 
  filter(exp$txID %in% venPriorityGenes$txid)

allVen <- ggplot(venom.exp,aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5 ) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  # scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('e) All Venom Genes') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))

allVen

svmp <- ggplot(subset(venom.exp,str_detect(txID,'SVMP')),aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  # scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('f) SVMPs') +
  xlab('H3K4me3 Peak Near Promoter') +
  ylab('VST Normalized Gene Expression')+
  theme_bw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

svsp <- ggplot(subset(venom.exp,str_detect(txID,'SVSP')),aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  # scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  ggtitle('g) SVSPs') +
  xlab('H3K4me3 Peak Near Promoter') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))

svsp

# svsp.2removed <- ggplot(venom.exp %>% filter(str_detect(txID,'SVSP') & txID != 'crovir-transcript-SVSP_2'),aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
#   geom_boxplot(show.legend = F) +
#   geom_jitter(width = 0.1,size=1,show.legend = F) +
#   geom_signif(comparisons = list(c("FALSE", "TRUE")),
#               map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = 0.5) +
#   scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
#   scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
#   ggtitle('g) SVSPs') +
#   xlab('H3K4me3 Peak Near Promoter') +
#   theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))
# 
# svsp.2removed

pla2 <- ggplot(subset(venom.exp,str_detect(txID,'PLA2') & txID != 'crovir-transcript-PLA2_gIIE'),aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0 ,vjust = .5) +
  # scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  ggtitle('h) PLA2s') +
  xlab('H3K4me3 Peak Near Promoter') +
  ylab('VST Normalized Gene Expression')+
  theme_bw() + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

other <- ggplot(subset(venom.exp,str_detect(txID,'SVSP|SVSP|PLA2|128',negate = T)),aes(x=ChIPsignal,y=avg1DPE,fill=ChIPsignal)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.1,size=1,show.legend = F) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")),
              map_signif_level = TRUE, textsize=4, tip_length = 0,vjust = 0.5) +
  scale_fill_manual(values=c('FALSE'='grey','TRUE'='skyblue2')) +
  # scale_y_continuous(limits=c(0,max(log10(venom.exp$avg1DPE+1)*1.2))) +
  xlab('H3K4me3 Peak Near Promoter') +
  ggtitle('i) Other VGs') +
  theme_bw() + theme(axis.title.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold'))

#(allGenes + allVen) / (svmp | svsp | pla2 | other) 

#(allGenes | allVen | svmp | svsp | pla2 | other) + plot_annotation(tag_levels = 'a' )


(allGenes | allVen) / (svmp | svsp ) / (pla2 | other) 

