
library(tidyverse)
library(patchwork)
library(nationalparkcolors)

# setwd('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/transcription_factors/_tfbs_analysis/scripts/tfbs_aligner/')
# 
# 
# tf.raw <- read_tsv('example_files/tfbs.gff',col_names = F,skip = 1) %>% 
#   mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])
# 
# 
# tf.ali <- read_tsv('./example_outs/aligned_input2_alignTFBS.txt',col_names = F) %>% 
#   mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
#   mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])
# 
# cons_score <- read_tsv('example_outs/aligned_input2_consScore.txt',col_names = F)
# 
# gaps <- read_tsv('example_outs/aligned_input2_gaps.txt',col_names = F) %>% 
#   mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])
# 
# 
# p1 <- ggplot(tf.raw,aes(x=X4,xend=X5,y=X1,yend=X1,color=X3)) +
#   geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(X5+20)),lwd=0.5,color='grey') +
#   geom_segment(lwd=10,alpha=0.5,show.legend = F) +
#   ggtitle('TFBS Predictions on Raw Enhancer Sequences') +
#   ylab('vPER') +
#   theme_classic() + theme(plot.title.position = 'plot')
# 
# p2 <- ggplot(tf.ali,aes(x=X2,xend=X3,y=X1,yend=X1,color=tfbs_id)) +
#   geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(X3+20)),lwd=0.5,color='grey') +
#   geom_segment(inherit.aes = F,data=gaps,aes(x=X2,xend=X3,y=X1,yend=X1),color='white',lwd=2,alpha=0.85) +
#   geom_segment(lwd=10,alpha=0.5) +
#   ggtitle('TFBS Positions After Sequence Alignment') +
#   ylab('vPER') +
#   theme_classic() + theme(plot.title.position = 'plot')
# 
# 
# p3 <- ggplot(cons_score,aes(x=X1,y=X2)) +
#   geom_line() +
#   scale_y_continuous(limits=c(0,1))+
#   scale_x_continuous(limits=c(0,max(tf.ali$X3+20)))+
#   ylab('Consensus Score') +
#   theme_classic() + theme(panel.grid.major.y = element_line())
# 
# 
# # p1 / p2  / p3
# 
# p1 / p2 / p3 + plot_layout(heights = c(3,3,1))
# 
# 


setwd('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/transcription_factors/_tfbs_analysis/_Enhancer_CIIIDER_TFBS_EnrichmentAnalyses_12.11.20/SVSP.vs.NonVen/mafft_aligned')


# tf.raw <- read_tsv('./hSVMP.vs.paralogEnriched_allSeqs.gff',col_names = F,skip = 1) %>% 
  # mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])


tf.ali <- read_tsv('~/Desktop/test/SVSP_PERs_seqs.aln.mafft.singleLine_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])

cons_score <- read_tsv('~/Desktop/test/SVSP_PERs_seqs.aln.mafft.singleLine_consScore.txt',col_names = F)

gaps <- read_tsv('~/Desktop/test/SVSP_PERs_seqs.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])


# p1 <- ggplot(tf.raw,aes(x=X4,xend=X5,y=X1,yend=X1,color=X3)) +
#   geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(X5+20)),lwd=0.5,color='grey') +
#   geom_segment(lwd=5,alpha=0.5,show.legend = F) +
#   scale_x_continuous(expand=c(0,0))+
#   ggtitle('TFBS Predictions on Raw Sequences') +
#   ylab('Promoter') +
#   xlab('')+
#   theme_classic() + theme(plot.title.position = 'plot')

p2 <- ggplot(tf.ali,aes(x=X2,xend=X3,y=X1,yend=X1,color=tfbs_id)) +
  geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(cons_score$X1+20)),lwd=0.75,color='grey') +
  geom_segment(inherit.aes = F,data=gaps,aes(x=X2,xend=X3,y=X1,yend=X1),color='white',lwd=2,alpha=0.85) +
  geom_segment(lwd=6,alpha=0.5) +
  ggtitle('TFBS Positions After Sequence Alignment') +
  scale_x_continuous(expand=c(0,0))+
  # scale_color_manual(values=c('#33561F','#39B481','#F88427','#BF301B')) +
  ylab('Promoter') +
  xlab('')+
  theme_classic() + theme(plot.title.position = 'plot')
p2

p3 <- ggplot(cons_score,aes(x=X1,y=X2)) +
  geom_line() +
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(expand=c(0,0),limits=c(0,max(cons_score$X1+20)))+
  ylab('Consensus Score') +
  xlab('Position')+
  theme_classic() + theme(panel.grid.major.y = element_line())


# p1 / p2  / p3

#p1 / p2 / p3 + plot_layout(heights = c(3,3,1))

p2 / p3 + plot_layout(heights = c(4,1))




#### SVSP Promoters


setwd('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/transcription_factors/_tfbs_analysis/_Promoter_CIIIDER_TFBS_HierarchicalEnrichmentAnalyses_09.29.20/7_HighestExpressed.vs.Paralogs/hSVSP.vs.paralogs/')


tf.raw <- read_tsv('./hSVSP.vs.paralogEnriched_allSeqs.gff',col_names = F,skip = 1) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])


tf.ali <- read_tsv('./SVSPandParalogs_promSeqs.align.fixed_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])

cons_score <- read_tsv('./SVSPandParalogs_promSeqs.align.fixed_consScore.txt',col_names = F)

gaps <- read_tsv('SVSPandParalogs_promSeqs.align.fixed_gaps.txt',col_names = F) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])


p1 <- ggplot(tf.raw,aes(x=X4,xend=X5,y=X1,yend=X1,color=X3)) +
  geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(X5+20)),lwd=0.5,color='grey') +
  geom_segment(lwd=5,alpha=0.5,show.legend = F) +
  scale_x_continuous(expand=c(0,0))+
  ggtitle('TFBS Predictions on Raw Sequences') +
  ylab('Promoter') +
  xlab('')+
  theme_classic() + theme(plot.title.position = 'plot')

p2 <- ggplot(tf.ali,aes(x=X2,xend=X3,y=X1,yend=X1,color=tfbs_id)) +
  geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(cons_score$X1+20)),lwd=0.75,color='grey') +
  geom_segment(inherit.aes = F,data=gaps,aes(x=X2,xend=X3,y=X1,yend=X1),color='white',lwd=2,alpha=0.85) +
  geom_segment(lwd=6,alpha=0.5) +
  ggtitle('TFBS Positions After Sequence Alignment') +
  scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=c('#33561F','#39B481','#F88427','#BF301B','#443020')) +
  ylab('Promoter') +
  xlab('')+
  theme_classic() + theme(plot.title.position = 'plot')


p3 <- ggplot(cons_score,aes(x=X1,y=X2)) +
  geom_line() +
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(expand=c(0,0),limits=c(0,max(cons_score$X1+20)))+
  ylab('Consensus Score') +
  xlab('Position')+
  theme_classic() + theme(panel.grid.major.y = element_line())


# p1 / p2  / p3

#p1 / p2 / p3 + plot_layout(heights = c(3,3,1))

p2 / p3 + plot_layout(heights = c(4,1))




#### PLA2 Promoters


setwd('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/transcription_factors/_tfbs_analysis/_Promoter_CIIIDER_TFBS_HierarchicalEnrichmentAnalyses_09.29.20/4_VenomFamilies.vs.NonVenGenes/PLA2/')


tf.raw <- read_tsv('./PLA2.vs.NonVenEnriched_allSeqs.gff',col_names = F,skip = 1) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])


tf.ali <- read_tsv('./PLA2andParalogs_promSeqs.align.fixed_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])

cons_score <- read_tsv('./PLA2andParalogs_promSeqs.align.fixed_consScore.txt',col_names = F)

gaps <- read_tsv('PLA2andParalogs_promSeqs.align.fixed_gaps.txt',col_names = F) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1])


p1 <- ggplot(tf.raw,aes(x=X4,xend=X5,y=X1,yend=X1,color=X3)) +
  geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(X5+20)),lwd=0.5,color='grey') +
  geom_segment(lwd=5,alpha=0.5,show.legend = F) +
  scale_x_continuous(expand=c(0,0))+
  ggtitle('TFBS Predictions on Raw Sequences') +
  ylab('Promoter') +
  xlab('')+
  theme_classic() + theme(plot.title.position = 'plot')

p2 <- ggplot(tf.ali,aes(x=X2,xend=X3,y=X1,yend=X1,color=tfbs_id)) +
  geom_segment(inherit.aes = F,aes(y=X1,yend=X1,x=0,xend=max(cons_score$X1+20)),lwd=0.75,color='grey') +
  geom_segment(inherit.aes = F,data=gaps,aes(x=X2,xend=X3,y=X1,yend=X1),color='white',lwd=2,alpha=0.85) +
  geom_segment(lwd=6,alpha=0.5) +
  ggtitle('TFBS Positions After Sequence Alignment') +
  scale_x_continuous(expand=c(0,0))+
  # scale_color_manual(values=c('#33561F','#39B481','#F88427','#BF301B','#443020')) +
  ylab('Promoter') +
  xlab('')+
  theme_classic() + theme(plot.title.position = 'plot')


p3 <- ggplot(cons_score,aes(x=X1,y=X2)) +
  geom_line() +
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(expand=c(0,0),limits=c(0,max(cons_score$X1+20)))+
  ylab('Consensus Score') +
  xlab('Position')+
  theme_classic() + theme(panel.grid.major.y = element_line())


# p1 / p2  / p3

#p1 / p2 / p3 + plot_layout(heights = c(3,3,1))

p2 / p3 + plot_layout(heights = c(4,1))
