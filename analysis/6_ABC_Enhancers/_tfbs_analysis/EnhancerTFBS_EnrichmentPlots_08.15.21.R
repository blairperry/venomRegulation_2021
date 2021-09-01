library(scico)
library(tidyverse)
library(ggvenn)
library(patchwork)

tf.families <- read_tsv('analysis/3_transcription_factors/Motif_TFfamilies_08.17.21.txt',col_names = c('MA_ID','symbol','family')) %>% 
  mutate(symbol = str_split_fixed(symbol,'[(]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[::]',2)[,1]) %>% 
  mutate(symbol = str_to_upper(symbol)) %>% 
  select(symbol,family)


svmp.enrich <- read_csv('analysis/6_ABC_Enhancers/_tfbs_analysis/SVMP.vs.NonVen/Enrichment: Text8356234261001362229_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up') %>% 
  mutate(vgfamily = 'SVMP')%>% 
  mutate(symbol = str_to_upper(transcription_factor_name)) %>% 
  left_join(tf.families)


svsp.enrich <- read_csv('analysis/6_ABC_Enhancers/_tfbs_analysis/SVSP.vs.NonVen/Enrichment: Text8832593824906450786_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')%>% 
  mutate(vgfamily = 'SVSP')%>% 
  mutate(symbol = str_to_upper(transcription_factor_name)) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'[::]'),str_split_fixed(symbol,'[::]',2)[,1],symbol)) %>% 
  left_join(tf.families)


pla2.enrich <- read_csv('analysis/6_ABC_Enhancers/_tfbs_analysis/PLA2.vs.NonVen/Enrichment: Text7970206761508930896_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')%>% 
  mutate(vgfamily = 'PLA2')%>% 
  mutate(symbol = str_to_upper(transcription_factor_name)) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'[::]'),str_split_fixed(symbol,'[::]',2)[,1],symbol)) %>% 
  left_join(tf.families)


forVenn <- list('SVMP'=svmp.enrich$transcription_factor_name,
                  'SVSP'=svsp.enrich$transcription_factor_name,
                  'PLA2'=pla2.enrich$transcription_factor_name)

ggvenn(forVenn,fill_color = c('#A3C0CF','#7DCBAC','#A23B5F'),stroke_size = 0.5)


all.enrich <- svmp.enrich %>% 
 bind_rows(svsp.enrich,pla2.enrich) %>% 
 mutate(vgfamily = factor(vgfamily,levels = c('SVMP','SVSP','PLA2'))) 

p.enrich <- ggplot(all.enrich,aes(x=vgfamily,y=reorder(transcription_factor_name,desc(transcription_factor_name)), fill=vgfamily)) +
  geom_point(pch=21,size=3,show.legend = F) + 
  labs(x='Venom Gene Family',y='TFBS',title = 'Enriched TFBS') +
  facet_grid(rows = vars(str_replace_all(family,'_',' ')),scales = 'free_y',space='free') +
  scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  theme_linedraw() + theme(plot.title.position = 'plot',strip.text.y = element_text(angle = 0,hjust = 0,colour = 'black'),strip.background = element_rect(fill='grey95',color='NA'),axis.text.x = element_text(angle=45,hjust=1))
p.enrich



# Bound/Unbound TFBS Plots ------------------------------------------------


exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv')
vg.info <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = c('chr','start','end','desc','messy_name','txid','symbol')) %>% 
  left_join(exp)

vg.exp <- vg.info %>% select(symbol,Median1DPE)


boundThresholds <- read_tsv('analysis/4_atacseq/tobias_footprinting/BoundThresholds_08.11.21.tsv')
rvg11.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG11') %>% select(threshold))
rvg1.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG1') %>% select(threshold))
rvg4.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG4') %>% select(threshold))

# SVMP --------------------------------------------------------------------

svmp.tf.ali <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_EnhancerPredictionsFull_08.18.21.aln.mafft.singleLine_alignTFBS.txt',
                         col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[_]',2)[,2]) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'R_'),str_split_fixed(symbol,'[_]',3)[,3],symbol)) %>% 
  left_join(vg.exp) %>% 
  mutate(IDsplit = ifelse(is.na(Median1DPE),str_split(symbol,'[.]'),'none')) %>% 
  unnest(IDsplit) %>% 
  left_join(vg.exp,by=c('IDsplit'='symbol')) %>% 
  group_by(PER_ID,start) %>% 
  mutate(meanExp = mean(Median1DPE.y)) %>% 
  mutate(Median1DPE = ifelse(is.na(Median1DPE.x),meanExp,Median1DPE.x)) %>% 
  mutate(meanExp = ifelse(is.na(meanExp),F,T)) %>% 
  select(-Median1DPE.y,-Median1DPE.x,-IDsplit) %>% 
  ungroup() %>% 
  arrange(-Median1DPE) %>% 
  mutate(PER_ID = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,14:16])) %>% 
  mutate(vgfamily = str_split_fixed(symbol,'[_]',2)[,1])



svsp.tf.ali <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVSP.vs.NonVen/tfbs_aligned/SVSP_EnhancerPredictionsFull_08.18.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[_]',2)[,2]) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'R_'),str_split_fixed(symbol,'[_]',3)[,3],symbol)) %>% 
  left_join(vg.exp) %>% 
  mutate(IDsplit = ifelse(is.na(Median1DPE),str_split(symbol,'[.]'),'none')) %>% 
  unnest(IDsplit) %>% 
  left_join(vg.exp,by=c('IDsplit'='symbol')) %>% 
  group_by(PER_ID,start) %>% 
  mutate(meanExp = mean(Median1DPE.y)) %>% 
  mutate(Median1DPE = ifelse(is.na(Median1DPE.x),meanExp,Median1DPE.x)) %>% 
  mutate(meanExp = ifelse(is.na(meanExp),F,T)) %>% 
  select(-Median1DPE.y,-Median1DPE.x,-IDsplit) %>% 
  ungroup() %>% 
  arrange(-Median1DPE) %>% 
  mutate(PER_ID = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,14:16])) %>% 
  mutate(vgfamily = str_split_fixed(symbol,'[_]',2)[,1])



pla2.tf.ali <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/PLA2_EnhancerPredictionsFull_08.18.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[_]',2)[,2]) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'R_'),str_split_fixed(symbol,'[_]',3)[,3],symbol)) %>% 
  left_join(vg.exp) %>% 
  mutate(IDsplit = ifelse(is.na(Median1DPE),str_split(symbol,'[.]'),'none')) %>% 
  unnest(IDsplit) %>% 
  left_join(vg.exp,by=c('IDsplit'='symbol')) %>% 
  group_by(PER_ID,start) %>% 
  mutate(meanExp = mean(Median1DPE.y)) %>% 
  mutate(Median1DPE = ifelse(is.na(Median1DPE.x),meanExp,Median1DPE.x)) %>% 
  mutate(meanExp = ifelse(is.na(meanExp),F,T)) %>% 
  select(-Median1DPE.y,-Median1DPE.x,-IDsplit) %>% 
  ungroup() %>% 
  arrange(-Median1DPE) %>% 
  mutate(PER_ID = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,14:16])) %>% 
  mutate(vgfamily = str_split_fixed(symbol,'[_]',2)[,1])

## Experimenting with other plot formats

tf.characterization <- read_tsv('analysis/3_transcription_factors/_functional_characterization/CandTFs_InterestingSubset_FULL_08.27.21.tsv') %>% 
  filter(value == T) 

all.boundTFcounts <- svmp.tf.ali %>% 
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
  # mutate(vgfamily = str_split_fixed(gene_id,'[ ]',2)[,1]) %>% 
  select(vgfamily,PER_ID,tfbs_id,avg_footprint,tot.bound,Median1DPE) %>% 
  unique() %>% 
  filter(tot.bound > 1) %>% 
  group_by(vgfamily,PER_ID,tfbs_id,Median1DPE) %>% 
  tally() %>% 
  arrange(Median1DPE) %>% 
  mutate(PER_ID = str_remove_all(PER_ID,'_R_')) %>%
  mutate(PER_ID = factor(PER_ID,levels=unique(.$PER_ID)))%>% 
  mutate(symbol = str_to_upper(tfbs_id)) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'[::]'),str_split_fixed(symbol,'[::]',2)[,1],symbol)) %>% 
  left_join(tf.families) %>% 
  mutate(vgfamily = factor(vgfamily,levels=c('SVMP','SVSP','PLA2'))) %>% 
  left_join(tf.characterization,by=c('symbol'='id')) %>% 
  mutate(feature = factor(feature,levels=c("Upregulated (RNA-seq)","SE-Associated (ChIP-seq)","Differentially Bound (ATAC-seq)","Previously implicated in venom regulation","Adrenoceptor Signaling Pathway","AP-1 Complex Member","Interacts with ERK","Unfolded Protein Response Pathway")))

#write_tsv(all.boundTFcounts,'analysis/6_ABC_Enhancers/_tfbs_analysis/_allBoundTFBS_vPERs_08.23.21.tsv')


all.exp <- svmp.tf.ali %>% 
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
  # mutate(vgfamily = str_split_fixed(gene_id,'[ ]',2)[,1]) %>% 
  select(vgfamily,PER_ID,Median1DPE,meanExp) %>% 
  unique() %>% 
  mutate(PER_ID = str_remove_all(PER_ID,'_R_')) %>% 
  filter(PER_ID %in% all.boundTFcounts$PER_ID) %>% 
  arrange(Median1DPE) %>% 
  mutate(PER_ID = factor(PER_ID,levels=unique(.$PER_ID)))%>% 
  mutate(vgfamily = factor(vgfamily,levels=c('SVMP','SVSP','PLA2')))


p1.alt <- ggplot(all.boundTFcounts,aes(x=reorder(PER_ID,-Median1DPE),y=tfbs_id,size=n,fill=vgfamily)) +
  geom_point(pch=23,show.legend = F) +
  scale_size_continuous(range = c(2,4)) +
  # scale_size_area(max_size = 4) +
  facet_grid(cols=vars(vgfamily),rows = vars(str_replace_all(family,'_',' ')),scales = 'free',space='free') +
  xlab('Venom Gene') + ylab('TFBS') +
  scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  # ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',
                           plot.title.position = 'plot',plot.title = element_text(face='bold'),strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank() )

p2.alt <- ggplot(all.exp,aes(y=log10(Median1DPE+1),x=reorder(PER_ID,-Median1DPE),fill=meanExp)) +
  geom_bar(stat='identity',show.legend = F) +
  facet_grid(cols=vars(vgfamily),scales = 'free_x',space='free') +
  ylab('Log10(TPM)') + 
  scale_fill_manual(values=c('FALSE'='firebrick4','TRUE'='lightpink')) +
  ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.title.x=element_blank(),
                           plot.title.position = 'plot',plot.title = element_text(face='bold'),strip.background = element_rect(fill='grey20'))


p3.alt <- ggplot(all.boundTFcounts,aes(x=feature,y=tfbs_id,fill=feature)) +
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


# p3.alt

p2.alt + plot_spacer() + p1.alt + p3.alt + plot_layout(heights = c(1,5),widths = c(3,1),ncol = 2)


# Enrichment plot from above but colored by whether there is evidence for bound sites

svmp.bound <- svmp.tf.ali %>% 
  filter(tot.bound>1) %>% 
  select(tfbs_id) %>% 
  mutate(tfbs_id = str_to_upper(tfbs_id)) %>% 
  unique()

svsp.bound <- svsp.tf.ali %>% 
  filter(tot.bound>1) %>% 
  select(tfbs_id) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_id,'[::]',2)[,1]) %>% 
  mutate(tfbs_id = str_to_upper(tfbs_id)) %>% 
  unique()

pla2.bound <- pla2.tf.ali %>% 
  filter(tot.bound>1) %>% 
  select(tfbs_id) %>% 
  mutate(tfbs_id = str_to_upper(tfbs_id)) %>% 
  unique()

all.enrich_bound <- all.enrich %>% 
  mutate(bound = ifelse(vgfamily == 'SVMP' & symbol %in% svmp.bound$tfbs_id,T,
                        ifelse(vgfamily == 'SVSP' & symbol %in% svsp.bound$tfbs_id,T,
                               ifelse(vgfamily == 'PLA2' & symbol %in% pla2.bound$tfbs_id,T,F))))

ggplot(all.enrich_bound,aes(x=vgfamily,y=reorder(transcription_factor_name,desc(transcription_factor_name)), fill=bound)) +
  geom_point(pch=21,size=3,show.legend = T) + 
  labs(x='Venom Gene Family',y='TFBS',title = 'Enriched TFBS') +
  facet_grid(rows = vars(str_replace_all(family,'_',' ')),scales = 'free_y',space='free') +
  # scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  theme_linedraw() + theme(plot.title.position = 'plot',strip.text.y = element_text(angle = 0,hjust = 0,colour = 'black'),strip.background = element_rect(fill='grey95',color='NA'),axis.text.x = element_text(angle=45,hjust=1))



# 
# # Similar to above but with unbound TFs too
# 
# 
# all.boundUnboundTFcounts <- svmp.tf.ali %>% 
#   bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
#   # mutate(vgfamily = str_split_fixed(PER_ID,'[ ]',2)[,1]) %>% 
#   select(vgfamily,PER_ID,tfbs_id,avg_footprint,tot.bound,Median1DPE) %>% 
#   unique() %>% 
#   mutate(bound = ifelse(tot.bound > 1,T,F)) %>% 
#   group_by(vgfamily,PER_ID,tfbs_id,Median1DPE,bound) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   arrange(Median1DPE) %>% 
#   mutate(PER_ID = factor(PER_ID,levels=unique(.$PER_ID)))
# 
# all.exp.v2 <- svmp.tf.ali %>% 
#   bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
#   # mutate(vgfamily = str_split_fixed(PER_ID,'[ ]',2)[,1]) %>% 
#   select(vgfamily,PER_ID,Median1DPE,meanExp) %>% 
#   unique() %>% 
#   filter(PER_ID %in% all.boundUnboundTFcounts$PER_ID) %>% 
#   arrange(Median1DPE) %>% 
#   mutate(PER_ID = factor(PER_ID,levels=unique(all.boundUnboundTFcounts$PER_ID)))
# 
# p1.v2 <- ggplot(all.boundUnboundTFcounts,aes(x=tfbs_id,y=PER_ID,size=n,color=tfbs_id)) +
#   geom_point(data=subset(all.boundUnboundTFcounts, bound == F),show.legend = F,pch=1,color='black',size=6) +
#   geom_point(data=subset(all.boundUnboundTFcounts, bound == T),show.legend = F,size=6) +
#   # scale_size_continuous(range = c(2,8)) +
#   scale_size_area(max_size = max(all.boundUnboundTFcounts$n)) +
#   facet_grid(rows=vars(vgfamily),scales = 'free_y',space='free') +
#   xlab('TFBS') + ylab('') +
#   ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
#   theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',strip.text.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold') )
# 
# p2.v2 <- ggplot(all.exp.v2,aes(x=log10(Median1DPE+1),y=PER_ID,fill=meanExp)) +
#   geom_bar(stat='identity') +
#   facet_grid(rows=vars(vgfamily),scales = 'free_y',space='free') +
#   ylab('') +
#   scale_fill_manual(values=c('FALSE'='firebrick4','TRUE'='lightpink')) +
#   theme_linedraw() + theme(panel.grid = element_blank(),
#                            axis.text.y = element_blank(),
#                            axis.title.y = element_blank())
# 
# p1.v2 + p2.v2 + plot_layout(widths = c(3,1))
# 
