
library(tidyverse)
library(ggvenn)
library(patchwork)

svmp.enrich <- read_csv('analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/Enrichment: Text5563405300722031706_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up') %>% 
  mutate(family = 'SVMP')

svsp.enrich <- read_csv('analysis/5_promoters/_tfbs_analysis/SVSP.vs.NonVen/Enrichment: Text1420245787213203463_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')%>% 
  mutate(family = 'SVSP')

pla2.enrich <- read_csv('analysis/5_promoters/_tfbs_analysis/PLA2.vs.NonVen/Enrichment: Text1784563758267191801_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')%>% 
  mutate(family = 'PLA2')

forVenn <- list('SVMP'=svmp.enrich$transcription_factor_name,
                  'SVSP'=svsp.enrich$transcription_factor_name,
                  'PLA2'=pla2.enrich$transcription_factor_name)

ggvenn(forVenn,fill_color = c('#A3C0CF','#7DCBAC','#A23B5F'),stroke_size = 0.5)

all.enrich <- svmp.enrich %>% 
 bind_rows(svsp.enrich,pla2.enrich) %>% 
 mutate(family = factor(family,levels = c('SVMP','SVSP','PLA2')))

ggplot(all.enrich,aes(x=family,y=reorder(transcription_factor_name,desc(transcription_factor_name)), fill=family)) +
  geom_point(pch=21,size=6,show.legend = F) + 
  labs(x='Venom Gene Family',y='TFBS',title = 'Enriched TFBS') +
  scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  theme_linedraw() + theme(plot.title.position = 'plot')



# Bound/Unbound TFBS Plots ------------------------------------------------

exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  select(txid,Median1DPE)

boundThresholds <- read_tsv('analysis/4_atacseq/tobias_footprinting/BoundThresholds_08.11.21.tsv')
rvg11.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG11') %>% select(threshold))
rvg1.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG1') %>% select(threshold))
rvg4.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG4') %>% select(threshold))

# SVMP --------------------------------------------------------------------

svmp.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_PromPeaks_seqs_08.14.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint','MA_ID')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(tx_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  left_join(exp,by=c('tx_id'='txid')) %>% 
  arrange(-Median1DPE) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' ')) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,15:17]))


svsp.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/svsp.vs.NonVen/tfbs_aligned/svsp_PromPeaks_seqs_08.14.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint','MA_ID')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(tx_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  left_join(exp,by=c('tx_id'='txid')) %>% 
  arrange(-Median1DPE) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' ')) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,15:17]))


pla2.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_PromPeaks_seqs_08.14.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint','MA_ID')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(tx_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  left_join(exp,by=c('tx_id'='txid')) %>% 
  arrange(-Median1DPE) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' ')) %>% 
  filter(gene_id != 'PLA2 gIIE') %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,15:17]))

## Experimenting with other plot formats

all.boundTFcounts <- svmp.tf.ali %>% 
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
  mutate(family = str_split_fixed(gene_id,'[ ]',2)[,1]) %>% 
  select(family,gene_id,tfbs_id,avg_footprint,tot.bound,Median1DPE) %>% 
  unique() %>% 
  filter(tot.bound > 1) %>% 
  group_by(family,gene_id,tfbs_id,Median1DPE) %>% 
  tally() %>% 
  arrange(Median1DPE) %>% 
  mutate(gene_id = factor(gene_id,levels=unique(.$gene_id)))

all.exp <- svmp.tf.ali %>% 
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
  mutate(family = str_split_fixed(gene_id,'[ ]',2)[,1]) %>% 
  select(family,gene_id,Median1DPE) %>% 
  unique() %>% 
  filter(gene_id %in% all.boundTFcounts$gene_id) %>% 
  arrange(Median1DPE) %>% 
  mutate(gene_id = factor(gene_id,levels=unique(.$gene_id)))

p1 <- ggplot(all.boundTFcounts,aes(x=tfbs_id,y=gene_id,size=n,color=tfbs_id)) +
  geom_point(show.legend = F) +
  # scale_size_continuous(range = c(2,8)) +
  scale_size_area(max_size = 8) +
  facet_grid(rows=vars(family),scales = 'free_y',space='free') +
  xlab('TFBS') + ylab('') +
  ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',strip.text.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold') )

p2 <- ggplot(all.exp,aes(x=log10(Median1DPE+1),y=gene_id)) +
  geom_bar(stat='identity',fill='firebrick4') +
  facet_grid(rows=vars(family),scales = 'free_y',space='free') +
  ylab('') +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.title.y = element_blank())

p1 + p2 + plot_layout(widths = c(3,1))


# Similar to above but with unbound TFs too


all.boundUnboundTFcounts <- svmp.tf.ali %>% 
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
  mutate(family = str_split_fixed(gene_id,'[ ]',2)[,1]) %>% 
  select(family,gene_id,tfbs_id,avg_footprint,tot.bound,Median1DPE) %>% 
  unique() %>% 
  mutate(bound = ifelse(tot.bound > 1,T,F)) %>% 
  group_by(family,gene_id,tfbs_id,Median1DPE,bound) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(Median1DPE) %>% 
  mutate(gene_id = factor(gene_id,levels=unique(.$gene_id)))

all.exp.v2 <- svmp.tf.ali %>% 
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>% 
  mutate(family = str_split_fixed(gene_id,'[ ]',2)[,1]) %>% 
  select(family,gene_id,Median1DPE) %>% 
  unique() %>% 
  filter(gene_id %in% all.boundUnboundTFcounts$gene_id) %>% 
  arrange(Median1DPE) %>% 
  mutate(gene_id = factor(gene_id,levels=unique(.$gene_id)))

p1.v2 <- ggplot(all.boundUnboundTFcounts,aes(x=tfbs_id,y=gene_id,size=n,color=tfbs_id)) +
  geom_point(data=subset(all.boundUnboundTFcounts, bound == F),show.legend = F,pch=1,color='black',size=6) +
  geom_point(data=subset(all.boundUnboundTFcounts, bound == T),show.legend = F,size=6) +
  # scale_size_continuous(range = c(2,8)) +
  scale_size_area(max_size = max(all.boundUnboundTFcounts$n)) +
  facet_grid(rows=vars(family),scales = 'free_y',space='free') +
  xlab('TFBS') + ylab('') +
  ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',strip.text.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(face='bold') )

p2.v2 <- ggplot(all.exp.v2,aes(x=log10(Median1DPE+1),y=gene_id)) +
  geom_bar(stat='identity',fill='firebrick4') +
  facet_grid(rows=vars(family),scales = 'free_y',space='free') +
  ylab('') +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.title.y = element_blank())

p1.v2 + p2.v2 + plot_layout(widths = c(3,1))

