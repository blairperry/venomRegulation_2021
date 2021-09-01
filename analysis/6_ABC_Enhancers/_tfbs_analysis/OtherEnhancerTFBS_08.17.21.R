
library(tidyverse)
library(ggvenn)
library(patchwork)
library(scico)

tf.families <- read_tsv('analysis/3_transcription_factors/Motif_TFfamilies_08.17.21.txt',col_names = c('MA_ID','symbol','family')) %>% 
  mutate(symbol = str_split_fixed(symbol,'[(]',2)[,1]) %>% 
  separate_rows(symbol,sep = '::') %>% 
  mutate(symbol = str_to_upper(symbol)) %>% 
  select(symbol,family)

vg.info <- read_tsv('data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = c('chr','start','end','desc','messy_id','txid','gene_id'))

exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  select(txid,Median1DPE)

boundThresholds <- read_tsv('analysis/4_atacseq/tobias_footprinting/BoundThresholds_08.11.21.tsv')
rvg11.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG11') %>% select(threshold))
rvg1.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG1') %>% select(threshold))
rvg4.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG4') %>% select(threshold))


tf.characterization <- read_tsv('analysis/3_transcription_factors/_functional_characterization/CandTFs_InterestingSubset_FULL_08.27.21.tsv') %>% 
  filter(value == T) 


other.tfscan <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/Other.scanOnly/Other_Enhancers_TFBS.gencoords.footprints.relPos.bed',col_names = F) %>% 
  select(PER_ID=1,tf_id = 3, tf_start = 5, RVG11_Footprint=11,RVG1_Footprint=12,RVG4_Footprint=13) %>% 
  mutate(avg_footprint = rowMeans(.[,3:5])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,8:10])) %>% 
  mutate(vg_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(vg_id = str_split_fixed(vg_id,pattern = '[-]',n=2)[,2]) %>% 
  left_join(vg.info,by=c('vg_id'='gene_id')) %>% 
  left_join(exp)

other.bound <- other.tfscan %>% 
  filter(tot.bound > 1) %>% 
  unique()

all.boundTFcounts <- other.bound %>% 
  select(PER_ID,txid,tf_id,avg_footprint,tot.bound) %>% 
  group_by(PER_ID,txid,tf_id) %>% 
  tally() %>% 
  mutate(symbol = str_to_upper(tf_id)) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'[::]'),str_split_fixed(symbol,'[::]',2)[,1],symbol)) %>% 
  left_join(tf.families) %>% 
  left_join(exp) %>% 
  arrange(-Median1DPE) %>% 
  mutate(PER_ID = factor(PER_ID,levels = unique(.$PER_ID))) %>% 
  left_join(tf.characterization,by=c('symbol'='id')) %>% 
  mutate(feature = factor(feature,levels=c("Upregulated (RNA-seq)","SE-Associated (ChIP-seq)","Differentially Bound (ATAC-seq)","Previously implicated in venom regulation","Adrenoceptor Signaling Pathway","AP-1 Complex Member","Interacts with ERK","Unfolded Protein Response Pathway")))




p1.alt <- ggplot(all.boundTFcounts,aes(x=reorder(PER_ID,-Median1DPE),y=tf_id,size=n)) +
  geom_point(pch=23,show.legend = F,fill='grey') +
  # scale_size_continuous(range = c(2,8)) +
  scale_size_area(max_size = 5) +
  facet_grid(rows = vars(str_replace_all(family,'_',' ')),scales = 'free',space='free') +
  xlab('Venom Gene') + ylab('TFBS') +
  # ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',
                           plot.title.position = 'plot',plot.title = element_text(face='bold'),strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank() )

p1.alt

p2.alt <- ggplot(all.boundTFcounts %>% select(PER_ID,Median1DPE) %>% unique,aes(y=log10(Median1DPE+1),x=reorder(PER_ID,-Median1DPE))) +
  geom_bar(stat='identity',fill='firebrick4') +
  ylab('Log10(TPM)') + 
  ggtitle('"Other" Venom Gene vPERs - Putatively-Bound TFBS') +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.title.x=element_blank(),
                           plot.title.position = 'plot',plot.title = element_text(face='bold'),strip.background = element_rect(fill='grey20'))

p3.alt <- ggplot(all.boundTFcounts,aes(x=feature,y=tf_id,fill=feature)) +
  geom_point(size=3,show.legend = F,pch=21) +
  # scale_size_continuous(range = c(2,8)) +
  # scale_size_area(max_size = 4) +
  facet_grid(rows = vars(str_replace_all(family,'_',' ')),scales = 'free',space='free') +
  # scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  scale_fill_scico_d() +
  xlab('') +
  ylab('') +
  # ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',
                           axis.text.y=element_blank(),
                           plot.title.position = 'plot',plot.title = element_text(face='bold'),strip.text.x = element_blank(),
                           strip.text.y = element_text(angle = 0,hjust = 0,colour = 'black'),strip.background = element_rect(fill='grey95',color='NA') )

p3.alt

p2.alt + plot_spacer() + p1.alt + p3.alt + plot_layout(heights = c(1,8),widths = c(3,1),ncol = 2)

