
library(tidyverse)
library(patchwork)

exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  select(txid,Median1DPE)

boundThresholds <- read_tsv('analysis/4_atacseq/tobias_footprinting/BoundThresholds_08.11.21.tsv')
rvg11.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG11') %>% select(threshold))
rvg1.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG1') %>% select(threshold))
rvg4.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG4') %>% select(threshold))
  
# SVMP --------------------------------------------------------------------

svmp.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_alignTFBS.txt',
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
  mutate(tot.bound = rowSums(.[,15:17])) %>% 
  mutate(gene_id = factor(gene_id,levels = unique(.$gene_id)))

svmp.cons_score <- read_tsv('./analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_consScore.txt',col_names = c('position','score'))

svmp.gaps <- read_tsv('./analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(tx_id = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  # mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  filter(tx_id %in% svmp.tf.ali$tx_id) %>% 
  mutate(tx_id = factor(tx_id,levels=unique(svmp.tf.ali$tx_id))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  select(4,2,3) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' '))%>% 
  mutate(gene_id = factor(gene_id,levels = unique(svmp.tf.ali$gene_id)))

svmp.splits <- svmp.tf.ali %>%
  filter(str_detect(tfbs_split,'split')) %>%
  mutate(split_num = str_split_fixed(tfbs_split,'_',2)[,2]) %>%
  group_by(PER_ID,tfbs_id,avg_footprint,strand) %>%
  mutate(split_start = min(end), split_end = max(start)) %>%
  select(tfbs=tfbs_id,split_start,split_end) %>%
  unique() %>% 
  mutate(tx_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' '))%>% 
  mutate(gene_id = factor(gene_id,levels = unique(svmp.tf.ali$gene_id)))

svmp.p1 <- ggplot(svmp.tf.ali,aes(xmin=start,xmax=end,ymin=0,ymax=avg_footprint*strand,fill=tfbs_id)) +
  geom_segment(inherit.aes = F,data=svmp.splits,aes(x=split_start,xend=split_end,y=0.5*avg_footprint*strand,yend=0.5*avg_footprint*strand),lwd=0.75,alpha=0.75) +
  geom_rect(data=subset(svmp.tf.ali,tot.bound > 1),alpha=0.95,col='black',lwd=0.25) +
  geom_rect(data=subset(svmp.tf.ali,tot.bound <= 1),alpha=0.5,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(svmp.cons_score$position+20)),lwd=0.75,color='grey40') +
  geom_segment(data=svmp.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=1,alpha=0.85) +
  facet_wrap(~gene_id,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  # coord_cartesian(xlim =c(0,max(svmp.tf.ali$end + 15)),clip = 'on',expand = F)+
  coord_cartesian(xlim =c(150,400),clip = 'on',expand = F)+
  # scale_fill_identity() +
  # scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        strip.placement = 'outside',
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='NA',size = 0.5,color = NA),
        panel.spacing = unit(0.25, "lines"))
svmp.p1

svmp.p2 <- ggplot(svmp.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  # coord_cartesian(xlim =c(0,max(svmp.tf.ali$end + 15)),clip = 'on',expand = F)+
  coord_cartesian(xlim =c(150,400),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

svmp.p3 <- ggplot(svmp.tf.ali %>% select(gene_id,Median1DPE) %>% unique(),aes(x=reorder(gene_id,desc(gene_id)),y=Median1DPE)) +
  geom_bar(stat='identity',fill='firebrick4') +
  coord_flip() +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(svmp.tf.ali$Median1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.title.y = element_blank(),axis.title.x = element_blank())

svmp.all <- svmp.p1 + svmp.p3 + svmp.p2 + plot_spacer() + plot_layout(heights = c(4,1),widths = c(5,1),nrow = 2)
svmp.all


# SVSP --------------------------------------------------------------------

svsp.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/svsp.vs.NonVen/tfbs_aligned/svsp_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_alignTFBS.txt',
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
  mutate(tot.bound = rowSums(.[,15:17]))%>% 
  mutate(gene_id = factor(gene_id,levels = unique(.$gene_id)))

svsp.cons_score <- read_tsv('./analysis/5_promoters/_tfbs_analysis/svsp.vs.NonVen/tfbs_aligned/svsp_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_consScore.txt',col_names = c('position','score'))

svsp.gaps <- read_tsv('./analysis/5_promoters/_tfbs_analysis/svsp.vs.NonVen/tfbs_aligned/svsp_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(tx_id = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  # mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  filter(tx_id %in% svsp.tf.ali$tx_id) %>% 
  mutate(tx_id = factor(tx_id,levels=unique(svsp.tf.ali$tx_id))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  select(4,2,3) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' '))%>% 
  mutate(gene_id = factor(gene_id,levels = unique(svsp.tf.ali$gene_id)))

svsp.splits <- svsp.tf.ali %>%
  filter(str_detect(tfbs_split,'split')) %>%
  mutate(split_num = str_split_fixed(tfbs_split,'_',2)[,2]) %>%
  group_by(PER_ID,tfbs_id,avg_footprint,strand) %>%
  mutate(split_start = min(end), split_end = max(start)) %>%
  select(tfbs=tfbs_id,split_start,split_end) %>%
  unique() %>% 
  mutate(tx_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' '))%>% 
  mutate(gene_id = factor(gene_id,levels = unique(svsp.tf.ali$gene_id)))

svsp.p1 <- ggplot(svsp.tf.ali,aes(xmin=start,xmax=end,ymin=0,ymax=avg_footprint*strand,fill=tfbs_id)) +
  geom_segment(inherit.aes = F,data=svsp.splits,aes(x=split_start,xend=split_end,y=0.5*avg_footprint*strand,yend=0.5*avg_footprint*strand),lwd=0.75,alpha=0.75) +
  geom_rect(data=subset(svsp.tf.ali,tot.bound > 1),alpha=0.65,col='black',lwd=0.25) +
  geom_rect(data=subset(svsp.tf.ali,tot.bound <= 1),alpha=0.25,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(svsp.cons_score$position+20)),lwd=0.75,color='grey40') +
  geom_segment(data=svsp.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=1,alpha=0.85) +
  facet_wrap(~gene_id,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(svsp.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_identity() +
  # scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        strip.placement = 'outside',
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='NA',size = 0.5,color = NA),
        panel.spacing = unit(0.25, "lines"))
svsp.p1

svsp.p2 <- ggplot(svsp.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(svsp.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

svsp.p3 <- ggplot(svsp.tf.ali %>% select(gene_id,Median1DPE) %>% unique(),aes(x=reorder(gene_id,desc(gene_id)),y=Median1DPE)) +
  geom_bar(stat='identity',fill='firebrick4') +
  coord_flip() +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(svsp.tf.ali$Median1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

svsp.all <- svsp.p1 + svsp.p3 + svsp.p2 + plot_spacer() + plot_layout(heights = c(4,1),widths = c(5,1),nrow = 2)
svsp.all


# pla2 --------------------------------------------------------------------

pla2.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_alignTFBS.txt',
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
  mutate(tot.bound = rowSums(.[,15:17]))%>% 
  mutate(gene_id = factor(gene_id,levels = unique(.$gene_id)))

pla2.cons_score <- read_tsv('./analysis/5_promoters/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_consScore.txt',col_names = c('position','score'))

pla2.gaps <- read_tsv('./analysis/5_promoters/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_PromPeaks_seqs_08.16.21.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(tx_id = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  # mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  filter(tx_id %in% pla2.tf.ali$tx_id) %>% 
  mutate(tx_id = factor(tx_id,levels=unique(pla2.tf.ali$tx_id))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  select(4,2,3) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' ')) %>% 
  filter(gene_id != 'PLA2 gIIE') %>% 
  mutate(gene_id = factor(gene_id,levels = unique(pla2.tf.ali$gene_id)))

pla2.splits <- pla2.tf.ali %>%
  filter(str_detect(tfbs_split,'split')) %>%
  mutate(split_num = str_split_fixed(tfbs_split,'_',2)[,2]) %>%
  group_by(PER_ID,tfbs_id,avg_footprint,strand) %>%
  mutate(split_start = min(end), split_end = max(start)) %>%
  select(tfbs=tfbs_id,split_start,split_end) %>%
  unique() %>% 
  mutate(tx_id = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(gene_id = str_remove(tx_id,'crovir-transcript-')) %>% 
  mutate(gene_id = str_replace_all(gene_id,'_',' ')) %>% 
  filter(gene_id != 'PLA2 gIIE') %>% 
  mutate(gene_id = factor(gene_id,levels = unique(pla2.tf.ali$gene_id)))

pla2.p1 <- ggplot(pla2.tf.ali,aes(xmin=start,xmax=end,ymin=0,ymax=avg_footprint*strand,fill=tfbs_id)) +
  geom_segment(inherit.aes = F,data=pla2.splits,aes(x=split_start,xend=split_end,y=0.5*avg_footprint*strand,yend=0.5*avg_footprint*strand),lwd=0.75,alpha=0.75) +
  geom_rect(data=subset(pla2.tf.ali,tot.bound > 1),alpha=0.65,col='black',lwd=0.25) +
  geom_rect(data=subset(pla2.tf.ali,tot.bound <= 1),alpha=0.25,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(pla2.cons_score$position+20)),lwd=0.75,color='grey40') +
  geom_segment(data=pla2.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=1,alpha=0.85) +
  facet_wrap(~gene_id,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(pla2.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_identity() +
  # scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        strip.placement = 'outside',
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='NA',size = 0.5,color = NA),
        panel.spacing = unit(0.25, "lines"))
pla2.p1

pla2.p2 <- ggplot(pla2.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(pla2.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

pla2.p3 <- ggplot(pla2.tf.ali %>% select(gene_id,Median1DPE) %>% unique(),aes(x=reorder(gene_id,desc(gene_id)),y=Median1DPE)) +
  geom_bar(stat='identity',fill='firebrick4') +
  coord_flip() +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(pla2.tf.ali$Median1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

pla2.all <- pla2.p1 + pla2.p3 + pla2.p2 + plot_spacer() + plot_layout(heights = c(4,1),widths = c(5,1),nrow = 2)
pla2.all




