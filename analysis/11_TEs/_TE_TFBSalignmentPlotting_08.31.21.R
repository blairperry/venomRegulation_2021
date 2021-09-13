
library(tidyverse)
library(patchwork)
library(colorspace)
library(ggrepel)
library(janitor)
library(randomcoloR)
library(ggtree)

tf.fams <- read_tsv('analysis/3_transcription_factors/Motif_TFfamilies_08.17.21.txt',col_names = c('MA_ID','tfbs','family')) %>% 
  mutate(tfbs = str_split_fixed(tfbs,'[(]',2)[,1]) %>% 
  mutate(tfbs = str_to_upper(tfbs))


te.ms.tf.ali <- read_tsv('./analysis/11_TEs/alignments/SVSP_FunctionalRegionTEOverlaps_wConsensus_seqs.mafftalign.singleLine_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1])%>% 
  mutate(X1 = str_replace_all(X1,'[.:;()]','_')) %>% 
  mutate(X1 = str_remove(X1,'hAT-Tip100__')) %>% 
  # mutate(X1 = factor(X1,levels=te.ms.upgma$id)) %>%
  unique() %>% 
  mutate(type = str_split_fixed(X1,'_',2)[,1]) %>% 
  mutate(type = ifelse(type == 'croMit','Consensus',type)) %>% 
  mutate(type = factor(type,levels = c('Consensus','Enhancer','Promoter','Other'))) %>% 
  mutate(X1 = str_remove(X1,'Cmitchellii_rnd-5_family-368#DNA/')) %>% 
  mutate(X1 = ifelse(str_detect(X1,'newConsensus'),'Cv1-hAT-Tip100 Consensus',X1)) %>% 
  mutate(tfbs_id = str_to_upper(tfbs_id)) %>% 
  left_join(tf.fams,by=c('tfbs_id'='tfbs')) %>% 
  mutate(clean_id = str_split_fixed(X1,'__',2)[,1] %>% str_remove_all('_R_')) %>% 
  mutate(group = str_split_fixed(clean_id,'_',2)[,1]) %>% 
  mutate(group = factor(group,levels = c('Cv1-hAT-Tip100 Consensus','Promoter','Enhancer','Other'))) %>% 
  mutate(region_start = as.numeric(str_split_fixed(clean_id,'[_-]',5)[,4])) %>% 
  arrange(desc(group),desc(region_start)) %>% 
  mutate(clean_id = factor(clean_id,levels = unique(.$clean_id))) 
  
  
te.ms.tf.ali



te.ms.cons_score <- read_tsv('./analysis/11_TEs/alignments/SVSP_FunctionalRegionTEOverlaps_wConsensus_seqs.mafftalign.singleLine_consScore.txt',col_names = F)


te.ms.gaps <- read_tsv('./analysis/11_TEs/alignments/SVSP_FunctionalRegionTEOverlaps_wConsensus_seqs.mafftalign.singleLine_gaps.txt',col_names = F) %>%
  # mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>%
  mutate(X1 = str_replace_all(X1,'[.:;()]','_')) %>% 
  mutate(X1 = str_remove(X1,'hAT-Tip100__')) %>% 
  # mutate(X1 = factor(X1,levels=te.ms.upgma$id)) %>%
  unique() %>%
  mutate(X2 = X2 - 1, X3 = X3 + 1) %>% 
  mutate(X1 = str_remove(X1,'Cmitchellii_rnd-5_family-368#DNA/')) %>% 
  mutate(X1 = ifelse(str_detect(X1,'newConsensus'),'Cv1-hAT-Tip100 Consensus',X1)) %>% 
  mutate(X1 = factor(X1,levels = unique(te.ms.tf.ali$X1))) %>% 
  mutate(clean_id = str_split_fixed(X1,'__',2)[,1] %>% str_remove_all('_R_')) %>% 
  mutate(clean_id = factor(clean_id,levels = unique(te.ms.tf.ali$clean_id))) 




te.ms.splits <- te.ms.tf.ali %>%
  filter(str_detect(X4,'split')) %>%
  mutate(split_num = str_split_fixed(X4,'_',2)[,2]) %>%
  group_by(X1,X5,X6) %>%
  mutate(split_start = min(X3), split_end = max(X2)) %>%
  select(tfbs=tfbs_id,split_start,split_end,X8,clean_id,group) %>%
  unique() 

 
# # Full Plot
# ggplot(te.ms.tf.ali %>% filter(str_detect(X1,'8924100',negate = T)),aes(x=X2,xend=X3,y=X1,yend=X1,color=tfbs_id)) +
#   # geom_segment(inherit.aes = F,data=te.ms.splits,aes(x=split_start,xend=split_end,y=X1,yend=X1,color=tfbs),lwd=0.25,alpha=0.75) +
#   geom_segment(alpha=0.65,lwd=3) +
#   # geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(te.ms.cons_score$X1+20)),lwd=0.75,color='grey60') +
#   # geom_segment(data=te.ms.tf.ali,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
#   ylab('')+
#   xlab('Relative Position') +
#   coord_cartesian(xlim =c(0,max(te.ms.tf.ali$X3 + 15)),clip = 'on',expand = T)+
#   # scale_y_continuous(position = 'right',limits = c(-abs(max(te.ms.tf.ali$X6)*1.2),abs(max(te.ms.tf.ali$X6)*1.2))) +
#   theme_classic() +
#   theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
#         legend.position = 'left',
#         strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
#         axis.title.x = element_blank(),
#         axis.text.y=element_text(size=8),
#         strip.background = element_rect(fill='grey90',size = 0.5),
#         panel.spacing = unit(0.25, "lines"))
# 


## Filtered plot based on TFs with enrichment in SVMP enhancers

te.enh.cand_tfs <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVSP/SVSP.boundTFs.tsv') %>% 
  select(tfbs_id) %>% 
  unique()

te.prom.cand_tfs <- read_tsv('analysis/5_promoters/_tfbs_analysis/allBoundTFBS_Promoters_08.24.21.tsv') %>% 
  filter(vgfamily == 'SVSP') %>% 
  select(tfbs_id) %>% 
  unique()

te.all.cand_tfs <- te.enh.cand_tfs %>% 
  bind_rows(te.prom.cand_tfs) %>% 
  unique()



focal.tf.ali <- te.ms.tf.ali %>% 
  filter(tfbs_id %in% str_to_upper(te.all.cand_tfs$tfbs_id))

focal.splits <- te.ms.splits %>% 
  filter(tfbs %in% te.all.cand_tfs$tfbs_id)

focal.gaps <- te.ms.gaps %>% 
  filter(X1 %in% focal.tf.ali$X1)



p.align <- ggplot(focal.tf.ali,aes(x=X2,xend=X3,y=clean_id,yend=clean_id,color=family)) +
  geom_segment(inherit.aes = F,aes(y=clean_id,yend=clean_id,x=0,xend=max(te.ms.cons_score$X1+20)),lwd=0.75,color='grey70') +
  geom_segment(data=focal.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=clean_id,yend=clean_id),color='white',lwd=2,alpha=0.85) +
  geom_segment(alpha=0.35,lwd=4) +
  ylab('')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(te.ms.tf.ali$X3)),clip = 'on',expand = T)+
  # scale_y_continuous(position = 'right',limits = c(-abs(max(te.ms.tf.ali$X6)*1.2),abs(max(te.ms.tf.ali$X6)*1.2))) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))

p.align



p.cons <- ggplot(te.ms.cons_score,aes(x=X1,y=X2)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(te.ms.tf.ali$X3)),clip = 'on',expand = T)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw() + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left', axis.title.x = element_blank())



p.cons + p.align +plot_layout(heights = c(1,10),nrow = 2,ncol=1)



### Footprint plots
# ggplot(focal.tf.ali.color,aes(xmin=X2,xmax=X3,ymin=0,ymax=X5*X6,fill=fam_col)) +
#   geom_rect(data=subset(focal.tf.ali.color,X6 > X7),color='black',lwd=0.25,alpha=0.5) +
#   geom_rect(data=subset(focal.tf.ali.color,str_detect(X1,'Consensus')),aes(ymin=-1,ymax=1),color='black',lwd=0.25,alpha=0.5) +
#   geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(te.ms.cons_score$X1+20)),lwd=0.5,color='grey60') +
#   geom_segment(data=focal.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
#   # geom_segment(inherit.aes = F,data=svmp.splits,aes(x=split_start,xend=split_end,y=0.5*X5*X6,yend=0.5*X5*X6,color=tfbs),lwd=1,alpha=0.5) +
#   facet_wrap(~X1,ncol = 1,strip.position = 'left',scales="free_y") +
#   ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
#   xlab('Relative Position') +
#   coord_cartesian(xlim =c(0,max(focal.tf.ali.color$X3 + 15)),clip = 'on',expand = F)+
#   scale_fill_identity() +
#   scale_y_continuous(position = 'right',limits = c(-abs(max(focal.tf.ali.color$X6)*1.2),abs(max(focal.tf.ali.color$X6)*1.2))) +
#   theme_classic() +
#   theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
#         legend.position = 'left',
#         strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
#         axis.title.x = element_blank(),
#         axis.text.y=element_text(size=8),
#         strip.background = element_rect(fill='white',size = 0),
#         strip.text = element_text(hjust=1),
#         panel.spacing = unit(0.25, "lines"))

ggplot(focal.tf.ali.color_noOutlier,
       aes(xmin=X2,xmax=X3,ymin=0,ymax=X5*X6,fill=fam_col)) +
  geom_rect(data=subset(focal.tf.ali.color_noOutlier,X6 > X7),color='black',lwd=0.25,alpha=0.5) +
  geom_rect(data=subset(focal.tf.ali.color_noOutlier,str_detect(X1,'Consensus')),aes(ymin=-4,ymax=4),color='black',lwd=0.25,alpha=0.5) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(te.ms.cons_score$X1+20)),lwd=0.5,color='grey60') +
  geom_segment(data=focal.gaps%>% filter(X1 != 'Enhancer_scaffold-mi2_8535725-8536125_+_'),
               inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
  # geom_segment(inherit.aes = F,data=svmp.splits,aes(x=split_start,xend=split_end,y=0.5*X5*X6,yend=0.5*X5*X6,color=tfbs),lwd=1,alpha=0.5) +
  facet_wrap(~X1,ncol = 1,strip.position = 'left',scales="fixed") +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(focal.tf.ali.color_noOutlier$X3 + 15)),clip = 'on',expand = F)+
  scale_fill_identity() +
  scale_y_continuous(position = 'right',limits = c(-abs(max(focal.tf.ali.color_noOutlier$X6)*1.2),abs(max(focal.tf.ali.color_noOutlier$X6)*1.2))) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='white',size = 0),
        strip.text = element_text(hjust=1),
        panel.spacing = unit(0.25, "lines"))


ggplot(focal.tf.ali.color_noOutlier %>% select(name,family,fam_col) %>% unique(),aes(x=1,y=name,fill=fam_col)) +
  geom_point(pch=22,size=5,alpha=0.5) +
  geom_point(pch=22,size=5,fill=NA) +
  facet_grid(rows=vars(family),scales='free_y',space = 'free_y')+
  scale_fill_identity() + 
  ggtitle('Legend') +
  theme_linedraw() + theme(panel.grid = element_blank(),strip.text.y = element_text(angle = 0,color='black',face='bold'),strip.background = element_rect(fill='white'),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())
