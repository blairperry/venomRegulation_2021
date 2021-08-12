
library(tidyverse)
library(patchwork)

exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_1DPEAvgExpression_08.08.21.tsv',col_names = c('txid','avg1DPE'))
vg.info <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = c('chr','start','end','desc','messy_name','txid','symbol')) %>% 
  left_join(exp)

vg.exp <- vg.info %>% select(symbol,avg1DPE)

boundThresholds <- read_tsv('analysis/4_atacseq/tobias_footprinting/BoundThresholds_08.11.21.tsv')
rvg11.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG11') %>% select(threshold))
rvg1.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG1') %>% select(threshold))
rvg4.thresh <- as.numeric(boundThresholds %>% filter(sample=='RVG4') %>% select(threshold))
  
# SVMP --------------------------------------------------------------------

svmp.tf.ali <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[_]',2)[,2]) %>% 
  left_join(vg.exp) %>% 
  mutate(IDsplit = ifelse(is.na(avg1DPE),str_split(symbol,'[.]'),'none')) %>% 
  unnest(IDsplit) %>% 
  left_join(vg.exp,by=c('IDsplit'='symbol')) %>% 
  group_by(PER_ID,start) %>% 
  mutate(meanExp = mean(avg1DPE.y)) %>% 
  mutate(avg1DPE = ifelse(is.na(avg1DPE.x),meanExp,avg1DPE.x)) %>% 
  mutate(meanExp = ifelse(is.na(meanExp),F,T)) %>% 
  select(-avg1DPE.y,-avg1DPE.x,-IDsplit) %>% 
  ungroup() %>% 
  arrange(-avg1DPE) %>% 
  mutate(PER_ID = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,14:16]))
  
svmp.cons_score <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_consScore.txt',col_names = c('position','score'))

svmp.gaps <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(PER_ID = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  filter(PER_ID %in% svmp.tf.ali$PER_ID) %>% 
  mutate(PER_ID = factor(PER_ID,levels=unique(svmp.tf.ali$PER_ID))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  select(4,2,3)

svmp.splits <- svmp.tf.ali %>%
  filter(str_detect(tfbs_split,'split')) %>%
  mutate(split_num = str_split_fixed(tfbs_split,'_',2)[,2]) %>%
  group_by(PER_ID,tfbs_id,avg_footprint,strand) %>%
  mutate(split_start = min(end), split_end = max(start)) %>%
  select(tfbs=tfbs_id,split_start,split_end) %>%
  unique()

svmp.p1 <- ggplot(svmp.tf.ali,aes(xmin=start,xmax=end,ymin=0,ymax=avg_footprint*strand,fill=tfbs_id)) +
  geom_segment(inherit.aes = F,data=svmp.splits,aes(x=split_start,xend=split_end,y=0.5*avg_footprint*strand,yend=0.5*avg_footprint*strand),lwd=0.75,alpha=0.75) +
  geom_rect(data=subset(svmp.tf.ali,tot.bound > 1),alpha=0.75,col='black',lwd=0.25) +
  geom_rect(data=subset(svmp.tf.ali,tot.bound <= 1),alpha=0.25,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(svmp.cons_score$position+20)),lwd=0.75,color='grey40') +
  geom_segment(data=svmp.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=1,alpha=0.85) +
  facet_wrap(~PER_ID,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(svmp.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_identity() +
  # scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
svmp.p1

svmp.p2 <- ggplot(svmp.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(svmp.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

svmp.p3 <- ggplot(svmp.tf.ali %>% select(PER_ID,avg1DPE,meanExp) %>% unique(),aes(x=reorder(PER_ID,desc(PER_ID)),y=avg1DPE,fill=meanExp)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=c('FALSE'='firebrick4','TRUE'='lightpink')) +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(svmp.tf.ali$avg1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

svmp.all <- svmp.p1 + svmp.p3 + svmp.p2 + plot_spacer() + plot_layout(heights = c(4,1),widths = c(5,1),nrow = 2)
svmp.all


# SVSP --------------------------------------------------------------------

SVSP.tf.ali <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVSP.vs.NonVen/tfbs_aligned/SVSP_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[_]',2)[,2]) %>% 
  left_join(vg.exp) %>% 
  mutate(IDsplit = ifelse(is.na(avg1DPE),str_split(symbol,'[.]'),'none')) %>% 
  unnest(IDsplit) %>% 
  left_join(vg.exp,by=c('IDsplit'='symbol')) %>% 
  group_by(PER_ID,start) %>% 
  mutate(meanExp = mean(avg1DPE.y)) %>% 
  mutate(avg1DPE = ifelse(is.na(avg1DPE.x),meanExp,avg1DPE.x)) %>% 
  mutate(meanExp = ifelse(is.na(meanExp),F,T)) %>% 
  select(-avg1DPE.y,-avg1DPE.x,-IDsplit) %>% 
  ungroup() %>% 
  arrange(-avg1DPE) %>% 
  mutate(PER_ID = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,14:16]))

SVSP.cons_score <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVSP.vs.NonVen/tfbs_aligned/SVSP_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_consScore.txt',col_names = c('position','score'))

SVSP.gaps <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/SVSP.vs.NonVen/tfbs_aligned/SVSP_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(PER_ID = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  filter(PER_ID %in% SVSP.tf.ali$PER_ID) %>% 
  mutate(PER_ID = factor(PER_ID,levels=unique(SVSP.tf.ali$PER_ID))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  select(4,2,3)

SVSP.splits <- SVSP.tf.ali %>%
  filter(str_detect(tfbs_split,'split')) %>%
  mutate(split_num = str_split_fixed(tfbs_split,'_',2)[,2]) %>%
  group_by(PER_ID,tfbs_id,avg_footprint,strand) %>%
  mutate(split_start = min(end), split_end = max(start)) %>%
  select(tfbs=tfbs_id,split_start,split_end) %>%
  unique()

SVSP.p1 <- ggplot(SVSP.tf.ali,aes(xmin=start,xmax=end,ymin=0,ymax=avg_footprint*strand,fill=tfbs_id)) +
  geom_segment(inherit.aes = F,data=SVSP.splits,aes(x=split_start,xend=split_end,y=0.5*avg_footprint*strand,yend=0.5*avg_footprint*strand),lwd=0.75,alpha=0.75) +
  geom_rect(data=subset(SVSP.tf.ali,tot.bound > 1),alpha=0.75,col='black',lwd=0.25) +
  geom_rect(data=subset(SVSP.tf.ali,tot.bound <= 1),alpha=0.25,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(SVSP.cons_score$position+20)),lwd=0.75,color='grey40') +
  geom_segment(data=SVSP.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=1,alpha=0.85) +
  facet_wrap(~PER_ID,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(SVSP.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_identity() +
  # scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
SVSP.p1

SVSP.p2 <- ggplot(SVSP.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(SVSP.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

SVSP.p3 <- ggplot(SVSP.tf.ali %>% select(PER_ID,avg1DPE,meanExp) %>% unique(),aes(x=reorder(PER_ID,desc(PER_ID)),y=avg1DPE,fill=meanExp)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=c('FALSE'='firebrick4','TRUE'='lightpink')) +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(SVSP.tf.ali$avg1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

SVSP.all <- SVSP.p1 + SVSP.p3 + SVSP.p2 + plot_spacer() + plot_layout(heights = c(4,1),widths = c(5,1),nrow = 2)
SVSP.all

# pla2 --------------------------------------------------------------------

pla2.tf.ali <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_alignTFBS.txt',
                        col_names = c('PER_ID','start','end','tfbs_split','strand','RVG11_Footprint','RVG1_Footprint','RVG4_Footprint')) %>% 
  mutate(tfbs_id = str_split_fixed(tfbs_split,'_',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(symbol = str_split_fixed(symbol,'[_]',2)[,2]) %>% 
  left_join(vg.exp) %>% 
  mutate(IDsplit = ifelse(is.na(avg1DPE),str_split(symbol,'[.]'),'none')) %>% 
  unnest(IDsplit) %>% 
  left_join(vg.exp,by=c('IDsplit'='symbol')) %>% 
  group_by(PER_ID,start) %>% 
  mutate(meanExp = mean(avg1DPE.y)) %>% 
  mutate(avg1DPE = ifelse(is.na(avg1DPE.x),meanExp,avg1DPE.x)) %>% 
  mutate(meanExp = ifelse(is.na(meanExp),F,T)) %>% 
  select(-avg1DPE.y,-avg1DPE.x,-IDsplit) %>% 
  ungroup() %>% 
  arrange(-avg1DPE) %>% 
  mutate(PER_ID = str_split_fixed(PER_ID,'[::]',2)[,1]) %>% 
  mutate(avg_footprint = rowMeans(.[,6:8])) %>% 
  mutate(RVG11.bound = ifelse(RVG11_Footprint >= rvg11.thresh,1,0),
         RVG1.bound = ifelse(RVG1_Footprint >= rvg1.thresh,1,0),
         RVG4.bound = ifelse(RVG4_Footprint >= rvg4.thresh,1,0)) %>% 
  mutate(tot.bound = rowSums(.[,14:16]))

pla2.cons_score <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_consScore.txt',col_names = c('position','score'))

pla2.gaps <- read_tsv('./analysis/6_ABC_Enhancers/_tfbs_analysis/pla2.vs.NonVen/tfbs_aligned/pla2_EnhancerPredictionsFull_08.10.21.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(PER_ID = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  filter(PER_ID %in% pla2.tf.ali$PER_ID) %>% 
  mutate(PER_ID = factor(PER_ID,levels=unique(pla2.tf.ali$PER_ID))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  select(4,2,3)

pla2.splits <- pla2.tf.ali %>%
  filter(str_detect(tfbs_split,'split')) %>%
  mutate(split_num = str_split_fixed(tfbs_split,'_',2)[,2]) %>%
  group_by(PER_ID,tfbs_id,avg_footprint,strand) %>%
  mutate(split_start = min(end), split_end = max(start)) %>%
  select(tfbs=tfbs_id,split_start,split_end) %>%
  unique()

pla2.p1 <- ggplot(pla2.tf.ali,aes(xmin=start,xmax=end,ymin=0,ymax=avg_footprint*strand,fill=tfbs_id)) +
  geom_segment(inherit.aes = F,data=pla2.splits,aes(x=split_start,xend=split_end,y=0.5*avg_footprint*strand,yend=0.5*avg_footprint*strand),lwd=0.75,alpha=0.75) +
  geom_rect(data=subset(pla2.tf.ali,tot.bound > 1),alpha=0.75,col='black',lwd=0.25) +
  geom_rect(data=subset(pla2.tf.ali,tot.bound <= 1),alpha=0.25,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(pla2.cons_score$position+20)),lwd=0.75,color='grey40') +
  geom_segment(data=pla2.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=1,alpha=0.85) +
  facet_wrap(~PER_ID,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(pla2.tf.ali$end + 15)),clip = 'on',expand = F)+
  # scale_fill_identity() +
  # scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
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

pla2.p3 <- ggplot(pla2.tf.ali %>% select(PER_ID,avg1DPE,meanExp) %>% unique(),aes(x=reorder(PER_ID,desc(PER_ID)),y=avg1DPE,fill=meanExp)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=c('FALSE'='firebrick4','TRUE'='lightpink')) +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(pla2.tf.ali$avg1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

pla2.all <- pla2.p1 + pla2.p3 + pla2.p2 + plot_spacer() + plot_layout(heights = c(4,1),widths = c(5,1),nrow = 2)
pla2.all



