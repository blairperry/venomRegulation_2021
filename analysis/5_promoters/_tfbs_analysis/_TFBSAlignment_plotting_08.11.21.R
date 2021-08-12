
library(tidyverse)
library(patchwork)
library(colorspace)
library(ggrepel)
library(janitor)
library(randomcoloR)
library(ggtree)
library(ggVennDiagram)


# Load in average expression ----------------------------------------------

exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_1DPEAvgExpression_08.08.21.tsv',col_names = c('txid','Avg1DPE'))

# Bound threshold (minimum score for any bound TFs in BINDetect results overview file) ---------------------------------------------------

# bound_thresh = 4.97183

# # Read in clusters 
# 
# clusters <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/transcription_factors/clustermotifs_output/_Motif_to_Clusters_12.18.20.tsv',col_names = c('tfbs_id','M_id','cluster')) 
# 
# focal_clusters <- read_tsv('focalTFBS_Clusters_01.07.21.txt')
# 
# set.seed(18)
# 
# cluster_colors <- clusters %>% 
#   filter(cluster %in% focal_clusters$cluster) %>% 
#   select(cluster) %>% 
#   unique() %>% 
#   mutate(cluster_col = distinctColorPalette(length(cluster)))
# 
# clusters <- clusters %>% 
#   left_join(cluster_colors) %>% 
#   mutate(cluster_col = ifelse(is.na(cluster_col),'black',cluster_col)) # set non-focal colors to black, mostly to make it obvious if something unexpected (i.e. non-focal TFBS) is being plotted
# 
# length(unique(clusters$cluster_col))

# Read in families

# tf_fams <- read_tsv('../../transcription_factors/_functional_characterization/Motif_TFfamilies_03.08.21.txt',col_names = c('id','name','family'))

# focal_fams <- read_tsv('focalTfFamilies_03.08.21.txt')

# fam_colors.fromEnhancers <- read_tsv('../../transcription_factors/_tfbs_analysis/_Enhancer_CIIIDER_TFBS_EnrichmentAnalyses_12.11.20/focalTfFamilies_COLORS_Enhancers_03.08.21.txt')
# 
# fam_colors <- focal_fams %>% 
#   left_join(fam_colors.fromEnhancers) %>% 
#   mutate(fam_col = ifelse(family == 'PAS_domain_factors','firebrick4',ifelse(family == 'TBX2-related_factors','lightsalmon4',fam_col)))

# write_tsv(fam_colors,'focalTfFamilies_COLORS_promoters_03.08.21.txt')

# ggplot(fam_colors,aes(x=1,y=family,fill=fam_col)) +
#   geom_point(pch=21,size=14,alpha=0.65) +
#   scale_fill_identity() +
#   theme_classic()

# tf_fams <- tf_fams %>%
#   left_join(fam_colors) %>%
#   mutate(fam_col = ifelse(is.na(fam_col),'black',fam_col)) # set non-focal colors to black, mostly to make it obvious if something unexpected (i.e. non-focal TFBS) is being plotted

# Enrichment Results Plot -------------------------------------------------

svmp.enrich <- read_csv('_tfbs_enrichment_binding/SVMP_vs_nonVenProms/Enrichment: Text4699459774371289374_MostSigDeficit.csv') %>% 
  clean_names() %>% 
  filter(gene_p_value < 0.05 & site_representation == 'Up') %>% 
  mutate(group = 'SVMP')

svsp.enrich <- read_csv('_tfbs_enrichment_binding/SVSP_vs_nonVenProms/Enrichment: Text3187430031164419848_MostSigDeficit.csv') %>% 
  clean_names() %>% 
  filter(gene_p_value < 0.05 & site_representation == 'Up') %>% 
  mutate(group = 'SVSP')

pla2.enrich <- read_csv('_tfbs_enrichment_binding/PLA2_vs_nonVenProms/Enrichment: Text5033974188740469399_MostSigDeficit.csv') %>% 
  clean_names() %>% 
  filter(gene_p_value < 0.05 & site_representation == 'Up') %>% 
  mutate(group = 'PLA2')

all.enrich <- svmp.enrich %>% 
  bind_rows(svsp.enrich,pla2.enrich) %>% 
  mutate(group = factor(group,levels=c('SVMP','SVSP','PLA2'))) %>% 
  add_count(transcription_factor_name) %>% 
  arrange(-n,transcription_factor_name) %>% 
  mutate(transcription_factor_name = factor(transcription_factor_name,levels = unique(.$transcription_factor_name)))

write_tsv(all.enrich,'allPromoterEnrichedTFBS_forSuppTable_04.06.21.tsv')
  
  


all.enrich %>% 
  group_by(group) %>% 
  tally()


all.enrich.venn <- list(SVMP = svmp.enrich$transcription_factor_name, SVSP = svsp.enrich$transcription_factor_name, PLA2 = pla2.enrich$transcription_factor_name)

ggVennDiagram(all.enrich.venn,color='white') + scale_fill_gradient(low='white',high='dodgerblue4')

ggplot(all.enrich,aes(x=group,y=reorder(transcription_factor_name,desc(transcription_factor_name)),size=significance_score,fill=significance_score)) +
  geom_point(pch=21) +
  scale_size_binned(name = 'Significance Score',range=c(3,8)) +
  scale_fill_binned_sequential(name= 'Significance Score','Mint') +
  ylab('Transcription Factor') +
  xlab('Venom Gene Family') +
  theme_linedraw()
#


# SVMP --------------------------------------------------------------------

svmp.upgma <- fortify(read.tree('promoterPeakSeqs/aligned/mafft/upgma_trees/SVMP_proms_nj.21010803115799.on.nh')) %>% 
  filter(isTip) %>% 
  mutate(id = str_split_fixed(label,'\\_\\_',2)[,1]) %>% 
  mutate(id = str_split_fixed(id, '\\-',3)[,3]) %>% 
  arrange(-y)

svmp.tf.ali <- read_tsv('./analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/tfbs_aligned/SVMP_PromPeaks_seqs_08.10.21.aln.mafft.singleLine_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  left_join(exp,by=c('X1'='gene')) %>% 
  arrange(-avg1DPE) %>% 
  mutate(X1 = factor(X1,levels=svmp.upgma$id)) %>% 
  mutate(bound = ifelse(X5 >= bound_thresh,'Bound','Unbound')) %>% 
  unique() %>% 
  left_join(tf_fams,by=c('tfbs_id'='name'))

svmp.cons_score <- read_tsv('./_tfbs_enrichment_binding/SVMP_vs_nonVenProms/mafft_aligned_TFBS/SVMP_PromPeaks_seqs.aln.mafft.singleLine_consScore.txt',col_names = F)

svmp.gaps <- read_tsv('./_tfbs_enrichment_binding/SVMP_vs_nonVenProms/mafft_aligned_TFBS/SVMP_PromPeaks_seqs.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  filter(X1 %in% svmp.tf.ali$X1) %>% 
  mutate(X1 = factor(X1,levels=unique(svmp.tf.ali$X1))) %>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1)

svmp.posScores <- read_tsv('./_tfbs_enrichment_binding/_TFBS_scoring/SVMP_PositionScores_12.22.20.tsv') %>% mutate(group = 'svmp') %>% 
  left_join(tf_fams,by=c('tfbs'='name')) %>% 
  mutate(high = ifelse(posScore >= mean(.$posScore),T,F))

svmp.posScores.mean <- mean(svmp.posScores$posScore)

# svmp.splits <- svmp.tf.ali %>% 
#   filter(str_detect(X4,'split')) %>% 
#   mutate(split_num = str_split_fixed(X4,'_',2)[,2]) %>% 
#   group_by(X1,X5,X6) %>% 
#   mutate(split_start = min(X3), split_end = max(X2)) %>% 
#   select(tfbs=tfbs_id,split_start,split_end,cluster_col) %>% 
#   unique()



svmp.p1 <- ggplot(svmp.tf.ali,aes(xmin=X2,xmax=X3,ymin=0,ymax=X5*X6,fill=fam_col)) +
  # geom_segment(inherit.aes = F,data=svmp.splits,aes(x=split_start,xend=split_end,y=0.5*X5*X6,yend=0.5*X5*X6,color=cluster_col),lwd=1,alpha=0.75) +
  geom_rect(data=subset(svmp.tf.ali,X6 > X7),alpha=0.65,col='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(svmp.cons_score$X1+20)),lwd=0.75,color='grey60') +
  geom_segment(data=svmp.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
  facet_wrap(~X1,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(svmp.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_y_continuous(position = 'right',limits = c(-abs(max(svmp.tf.ali$X6)*1.2),abs(max(svmp.tf.ali$X6)*1.2))) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
svmp.p1

svmp.p2 <- ggplot(svmp.cons_score,aes(x=X1,y=X2)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(svmp.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

svmp.p3 <- ggplot(svmp.tf.ali %>% select(X1,avg1DPE) %>% unique(),aes(x=reorder(X1,desc(X1)),y=avg1DPE)) +
  geom_bar(stat='identity',fill='firebrick4') +
  coord_flip() +
  # ylab('Gene Expression at 1DPE\n(Average Norm. Counts)') +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(svmp.tf.ali$avg1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

svmp.p4 <- ggplot(svmp.posScores,aes(xmin=start,xmax=stop,ymin=0,ymax=posScore,fill=fam_col)) +
  geom_rect(alpha=0.65,show.legend = F,col='black',lwd=0.25) +
  # geom_point(data=subset(svmp.posScores,high == T),aes(x=start,y=4)) +
  geom_hline(yintercept = mean(svmp.posScores$posScore),lty=2) +
  coord_cartesian(xlim =c(0,max(svmp.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_y_continuous(expand = c(0,0),limits=c(0,max(svmp.posScores$posScore)*1.1)) +
  geom_text_repel(data=svmp.posScores %>% filter(posScore >= mean(svmp.posScores$posScore)), aes(x=(start+stop)/2,y=posScore,label=tfbs,color=fam_col),max.overlaps = 15) +
  scale_fill_identity() + 
  scale_color_identity() +
  ylab('TFBS Position Score') +
  theme_classic(base_size = 16) + theme(axis.title.x = element_blank())
svmp.p4

svmp.p5 <- ggplot(svmp.tf.ali %>% select(family,fam_col) %>% unique(),aes(x=1,y=family,label=family,fill=fam_col)) +
  scale_fill_identity() +
  geom_label(alpha=0.65,) +
  theme_void()
svmp.p5

svmp.all <- svmp.p4 + plot_spacer() + svmp.p1 + svmp.p3 + svmp.p2 + plot_spacer() + plot_layout(heights = c(1,4,1),widths = c(5,1),nrow = 3)
svmp.all




# SVSP --------------------------------------------------------------------


svsp.upgma <- fortify(read.tree('promoterPeakSeqs/aligned/mafft/upgma_trees/SVSP_proms_nj.21010803157588.on.nh')) %>% 
  filter(isTip) %>% 
  mutate(id = str_split_fixed(label,'\\_\\_',2)[,1]) %>% 
  mutate(id = str_split_fixed(id, '\\-',3)[,3]) %>% 
  arrange(-y)

svsp.tf.ali <- read_tsv('./_tfbs_enrichment_binding/svsp_vs_nonVenProms/mafft_aligned_TFBS/SVSP_PromPeaks_seqs.aln.mafft.singleLine_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  left_join(exp,by=c('X1'='gene')) %>% 
  arrange(-avg1DPE) %>% 
  mutate(X1 = factor(X1,levels=svsp.upgma$id)) %>% 
  mutate(bound = ifelse(X5 >= bound_thresh,'Bound','Unbound')) %>% 
  unique() %>% 
  left_join(tf_fams,by=c('tfbs_id'='name'))

svsp.cons_score <- read_tsv('./_tfbs_enrichment_binding/svsp_vs_nonVenProms/mafft_aligned_TFBS/SVSP_PromPeaks_seqs.aln.mafft.singleLine_consScore.txt',col_names = F)

svsp.gaps <- read_tsv('./_tfbs_enrichment_binding/svsp_vs_nonVenProms/mafft_aligned_TFBS/SVSP_PromPeaks_seqs.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  filter(X1 %in% svsp.tf.ali$X1) %>% 
  mutate(X1 = factor(X1,levels=unique(svsp.tf.ali$X1)))%>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1)

svsp.posScores <- read_tsv('./_tfbs_enrichment_binding/_TFBS_scoring/SVSP_PositionScores_12.22.20.tsv') %>% mutate(group = 'svsp') %>% 
  left_join(tf_fams,by=c('tfbs'='name'))%>% 
  mutate(high = ifelse(posScore >= mean(.$posScore),T,F))

svsp.posScores.mean <- mean(svsp.posScores$posScore)

# Currently no splits in SVSP alignment 
# svsp.splits <- svsp.tf.ali %>%
#   filter(str_detect(X4,'split')) %>%
#   mutate(split_num = str_split_fixed(X4,'_',2)[,2]) %>%
#   group_by(X1,X5,X6) %>%
#   mutate(split_start = min(X3), split_end = max(X2)) %>%
#   select(tfbs=tfbs_id,split_start,split_end) %>%
#   unique()


svsp.p1 <- ggplot(svsp.tf.ali,aes(xmin=X2,xmax=X3,ymin=0,ymax=X5*X6,fill=fam_col)) +
  geom_rect(data=subset(svsp.tf.ali,X6 > X7),alpha=0.65,color='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(svsp.cons_score$X1+20)),lwd=0.75,color='grey60') +
  geom_segment(data=svsp.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
  facet_wrap(~X1,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(min(svsp.tf.ali$X2 - 15),max(svsp.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_y_continuous(position = 'right',limits = c(-abs(max(svsp.tf.ali$X6)*1.2),abs(max(svsp.tf.ali$X6)*1.2))) +
  scale_fill_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
svsp.p1


svsp.p2 <- ggplot(svsp.cons_score,aes(x=X1,y=X2)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(min(svsp.tf.ali$X2 - 15),max(svsp.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

svsp.p3 <- ggplot(svsp.tf.ali %>% select(X1,avg1DPE) %>% unique(),aes(x=reorder(X1,desc(X1)),y=avg1DPE)) +
  geom_bar(stat='identity',fill='firebrick4') +
  coord_flip() +
  # ylab('Gene Expression at 1DPE\n(Average Norm. Counts)') +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(svsp.tf.ali$avg1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())


svsp.p4 <- ggplot(svsp.posScores,aes(xmin=start,xmax=stop,ymin=0,ymax=posScore,fill=fam_col)) +
  geom_rect(alpha=0.65,color='black',lwd=0.25,show.legend = F) +
  geom_hline(yintercept = mean(svsp.posScores$posScore),lty=2) +
  coord_cartesian(xlim =c(min(svsp.tf.ali$X2 - 15),max(svsp.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_y_continuous(expand = c(0,0),limits=c(0,max(svsp.posScores$posScore)*1.1)) +
  geom_text_repel(data=svsp.posScores %>% filter(posScore >= mean(svsp.posScores$posScore)), aes(x=(start+stop)/2,y=posScore,label=tfbs,color=fam_col),max.overlaps = 15) +
  scale_fill_identity() + 
  scale_color_identity() +
  ylab('TFBS Position Score') +
  theme_classic(base_size = 16) + theme(axis.title.x = element_blank())
svsp.p4

svsp.all <- svsp.p4 + plot_spacer() + svsp.p1 + svsp.p3 + svsp.p2 + plot_spacer() + plot_layout(heights = c(1,4,1),widths = c(5,1),nrow = 3)
svsp.all


# pla2 --------------------------------------------------------------------
# 
# pla2.upgma <- fortify(read.tree('promoterPeakSeqs/aligned/mafft/upgma_trees/PLA2_nj.21010804039280.on.nh')) %>% 
#   filter(isTip) %>% 
#   mutate(id = str_split_fixed(label,'\\_\\_',2)[,1]) %>% 
#   mutate(id = str_split_fixed(id, '\\-',3)[,3]) %>% 
#   arrange(-y)

pla2.tf.ali <- read_tsv('./_tfbs_enrichment_binding/pla2_vs_nonVenProms/mafft_aligned_TFBS/PLA2_PromPeaks_noNVP_seqs.aln.mafft.singleLine_alignTFBS.txt',col_names = F) %>% 
  mutate(tfbs_id = str_split_fixed(X4,'_',2)[,1]) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  left_join(exp,by=c('X1'='gene')) %>% 
  arrange(-avg1DPE) %>% 
  mutate(X1 = factor(X1,levels=c('PLA2_A1','PLA2_B1','PLA2_C1'))) %>% 
  mutate(bound = ifelse(X5 >= bound_thresh,'Bound','Unbound')) %>% 
  # mutate(X2 = ifelse(X3-X2 == 0, X2-1,X2)) %>% # to Fix invisible 1bp long split segment of one TFBS
  unique() %>% 
  left_join(tf_fams,by=c('tfbs_id'='name'))

pla2.cons_score <- read_tsv('./_tfbs_enrichment_binding/pla2_vs_nonVenProms/mafft_aligned_TFBS/PLA2_PromPeaks_noNVP_seqs.aln.mafft.singleLine_consScore.txt',col_names = F)

pla2.gaps <- read_tsv('./_tfbs_enrichment_binding/pla2_vs_nonVenProms/mafft_aligned_TFBS/PLA2_PromPeaks_noNVP_seqs.aln.mafft.singleLine_gaps.txt',col_names = F) %>% 
  mutate(X1 = str_split_fixed(X1,'[::]',2)[,1]) %>% 
  mutate(X1 = str_remove(X1,'crovir-transcript-')) %>% 
  filter(X1 %in% pla2.tf.ali$X1) %>% 
  mutate(X1 = factor(X1,levels=unique(pla2.tf.ali$X1)))%>% 
  mutate(X2 = X2 - 1, X3 = X3 + 1)

pla2.posScores <- read_tsv('./_tfbs_enrichment_binding/_TFBS_scoring/PLA2_PositionScores_02.26.21.tsv') %>% mutate(group = 'pla2') %>% 
  left_join(tf_fams,by=c('tfbs'='name'))%>% 
  mutate(high = ifelse(posScore >= mean(.$posScore),T,F))

pla2.posScores.mean <- mean(pla2.posScores$posScore)


pla2.splits <- pla2.tf.ali %>% 
  filter(str_detect(X4,'split')) %>% 
  mutate(split_num = str_split_fixed(X4,'_',2)[,2]) %>% 
  group_by(X1,X5,X6) %>% 
  mutate(split_start = min(X3), split_end = max(X2)) %>% 
  select(tfbs=tfbs_id,split_start,split_end,fam_col) %>% 
  unique()


pla2.p1 <- ggplot(pla2.tf.ali,aes(xmin=X2,xmax=X3,ymin=0,ymax=X5*X6,fill=fam_col)) +
  geom_segment(inherit.aes = F,data=pla2.splits,aes(x=split_start,xend=split_end,y=0.5*X5*X6,yend=0.5*X5*X6,color=fam_col),lwd=1,alpha=0.75) +
  geom_rect(data=subset(pla2.tf.ali,X6 > X7),alpha=0.65,color='black',lwd=0.25) +
  geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(pla2.cons_score$X1+20)),lwd=0.75,color='grey60') +
  geom_segment(data=pla2.gaps,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
  facet_wrap(~X1,ncol = 1,strip.position = 'left') +
  ylab('Strand (Pos = 1, Neg = -1) * Footprint Score')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(pla2.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_y_continuous(position = 'right',limits = c(-abs(max(pla2.tf.ali$X6)*1.2),abs(max(pla2.tf.ali$X6)*1.2))) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_classic(base_size = 16) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
pla2.p1


pla2.p2 <- ggplot(pla2.cons_score,aes(x=X1,y=X2)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(pla2.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw(base_size = 16) + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

pla2.p3 <- ggplot(pla2.tf.ali %>% select(X1,avg1DPE) %>% unique(),aes(x=reorder(X1,avg1DPE),y=avg1DPE)) +
  geom_bar(stat='identity',fill='firebrick4') +
  coord_flip() +
  ylab('Gene Expression at 1DPE\n(Average Norm. Counts)') +
  scale_y_continuous(expand=c(0,0),limits=c(0,max(pla2.tf.ali$avg1DPE)*1.1),position = 'left') +
  theme_classic(base_size = 16) + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())

pla2.p4 <- ggplot(pla2.posScores,aes(xmin=start,xmax=stop,ymin=0,ymax=posScore,fill=fam_col)) +
  geom_rect(alpha=0.5,show.legend = F,color='black',lwd=0.25) +
  geom_hline(yintercept = mean(pla2.posScores$posScore),lty=2) +
  coord_cartesian(xlim =c(min(pla2.tf.ali$X2 - 15),max(pla2.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_y_continuous(expand = c(0,0),limits=c(0,max(pla2.posScores$posScore)*1.1)) +
  geom_text_repel(data=pla2.posScores %>% filter(posScore >= mean(pla2.posScores$posScore)), aes(x=(start+stop)/2,y=posScore,label=tfbs,color=fam_col),max.overlaps = 15) +
  scale_fill_identity() + 
  coord_cartesian(xlim =c(0,max(pla2.tf.ali$X3 + 15)),clip = 'on',expand = F)+
  scale_color_identity() +
  ylab('TFBS Position Score') +
  theme_classic(base_size = 16) + theme(axis.title.x = element_blank())
pla2.p4


pla2.all <- pla2.p4 + plot_spacer() + pla2.p1 + pla2.p3 + pla2.p2 + plot_spacer() + plot_layout(heights = c(1,4,1),widths = c(5,1),nrow = 3)
pla2.all


### Combined data frame for cluster color simplification

all.tf.ali.clusters <- svmp.tf.ali %>%
  bind_rows(svsp.tf.ali,pla2.tf.ali) %>%
  select(family) %>%
  unique()

# write_tsv(all.tf.ali.clusters,'focalTfFamilies_03.08.21.txt')


### Plot all 3 together
svmp.all
svsp.all
pla2.all

svmp.all / svsp.all / pla2.all






### Plot 1 promoter per family with only TFBS w/ Position Score > mean in a given family

highPosScores <- svmp.posScores %>% 
  bind_rows(svsp.posScores,pla2.posScores) %>% 
  mutate(group = str_to_upper(group)) %>% 
  mutate(group = factor(group, levels = c('SVMP','SVSP','PLA2'))) %>% 
  filter(high == T) %>% 
  arrange(start) %>% 
  mutate(tfbs = factor(tfbs,levels=unique(.$tfbs)))

highPosScores %>% 
  select(group,tfbs,posScore) %>% 
  group_by(group) %>% 
  filter(posScore == max(posScore))


# write_tsv(highPosScores,'./_tfbs_enrichment_binding/_TFBS_scoring/HighScoringTFBS_main3VGs_02.26.21.tsv')


highPosScores.splits <- highPosScores %>% 
  filter(str_detect(tfbs_pos_id,'split')) %>% 
  mutate(split_num = str_split_fixed(tfbs_pos_id,'_',3)[,2]) %>% 
  group_by(presentScore,tfbs,strand) %>% 
  mutate(split_start = min(stop), split_end = max(start)) %>% 
  select(tfbs,split_start,split_end,group,fam_col) %>% 
  unique()

ggplot(highPosScores,aes(x=start,xend=stop,y=tfbs,yend=tfbs,color=fam_col)) +
  geom_segment(data=highPosScores.splits,aes(x=split_start,xend=split_end,y=tfbs,yend=tfbs,color=fam_col),lwd=1,alpha=0.5,show.legend = F) +
  geom_segment(alpha=0.5,size=8,show.legend = F) +
  scale_color_identity() +
  facet_wrap(~group,scale='free',ncol = 1) +
  theme_linedraw() + theme(panel.grid.minor.x = element_blank())


### Similar plots but w/ Secondary Motifs as well

tf.exp <- read_tsv('../../transcription_factors/_TF_NormCounts_Avg1DPEandNonVen_01.05.20.tsv') %>% 
  mutate(id = str_to_upper(id)) 

tf.de <- read_tsv('../../transcription_factors/_atacseq_footprinting/combinedTFs_FIXED.txt') %>% 
  select(id,upreg) %>% 
  mutate(id = str_to_upper(id)) %>% 
  mutate(id = str_remove(id, '\\(VAR\\.2\\)')) %>% 
  filter(upreg == 'TRUE')


svmp.1st2nd <- read_tsv('_tfbs_enrichment_binding/SVMP_vs_nonVenProms/secondary_tfbs/SVMP_1stAnd2ndTFBS_PositionScores_02.26.21.tsv') %>% mutate(group = 'SVMP') %>% 
  mutate(highScore = ifelse(posScore > svmp.posScores.mean,T,F))

svsp.1st2nd <- read_tsv('_tfbs_enrichment_binding/SVSP_vs_nonVenProms/secondary_tfbs/SVSP_1stAnd2ndTFBS_PositionScores_02.26.21.tsv')%>% mutate(group = 'SVSP')%>% 
  mutate(highScore = ifelse(posScore > svsp.posScores.mean,T,F))

pla2.1st2nd <- read_tsv('_tfbs_enrichment_binding/PLA2_vs_nonVenProms/secondary_tfbs/PLA2_1stAnd2ndTFBS_PositionScores_02.26.21.tsv')%>% mutate(group = 'PLA2')%>% 
  mutate(highScore = ifelse(posScore > pla2.posScores.mean,T,F))

all.1st2nd <- svmp.1st2nd %>% 
  bind_rows(svsp.1st2nd,pla2.1st2nd) %>% 
  mutate(group = str_to_upper(group)) %>% 
  mutate(group = factor(group, levels = c('SVMP','SVSP','PLA2'))) %>% 
  # filter(posScore > 1) %>%
  filter(highScore == 'TRUE') %>% 
  arrange(start) %>% 
  mutate(tfbs = factor(tfbs,levels=unique(.$tfbs))) %>% 
  unique() %>% 
  group_by(tfbs,group) %>% 
  mutate(maxPosScore = max(posScore)) %>% 
  left_join(tf_fams,by=c('tfbs'='name')) %>% 
  group_by(family) %>% 
  mutate(order = length(unique(group))) %>% 
  arrange(-order) %>% 
  mutate(family = factor(family,levels=unique(.$family))) 
  
#write_tsv(all.1st2nd,'./_tfbs_enrichment_binding/All3_PrimSecHighTFs_02.26.21.tsv')

all.1st2nd.splits <- all.1st2nd %>% 
  filter(str_detect(tfbs_pos_id,'split')) %>% 
  mutate(split_num = str_split_fixed(tfbs_pos_id,'_',3)[,2]) %>% 
  group_by(presentScore,tfbs,strand) %>% 
  mutate(split_start = min(stop), split_end = max(start)) %>% 
  select(tfbs,split_start,split_end,group,round) %>% 
  unique()


all.1st2nd.exp <- all.1st2nd %>% 
  mutate(tfbs = str_to_upper(tfbs)) %>% 
  mutate(tfbs2 = tfbs) %>% 
  mutate(tfbs2 = ifelse(str_detect(tfbs2,'\\:'),str_split_fixed(tfbs2,'\\:\\:',2),tfbs2)) %>% 
  left_join(tf.exp,by=c('tfbs2'='id')) %>% 
  select(family,tfbs,avg1DPE) %>% 
  unique()



# plotA <- ggplot(all.1st2nd,aes(x=group,y=reorder(tfbs,order),fill=round,size=avgPosScore)) +
#   geom_point(pch=21) +
#   scale_size_continuous(range = c(4,8)) +
#   scale_fill_discrete_sequential(palette = "Blues", nmax = 6,order = c(2,6)) +
#   ylab('TFBS') +
#   xlab('Venom Family') +
#   theme_linedraw() + theme(panel.grid.minor.x = element_blank(),legend.position = 'right')
# plotA


p1 <- ggplot(all.1st2nd,aes(x=group,y=tfbs,fill=maxPosScore,shape=round,size=maxPosScore)) +
  geom_point() +
  scale_size_continuous(range = c(3,7)) +
  scale_shape_manual(values = c(21,23)) +
  # scale_fill_discrete_sequential(palette = "Blues", nmax = 6,order = c(2,6)) +
  scale_fill_continuous_sequential('Blues') +
  facet_grid(rows=vars(family),scales='free_y',space = 'free_y')+
  ylab('TFBS') +
  xlab('Venom Family') +
  theme_linedraw(base_size = 16) + theme(panel.grid.minor.x = element_blank(),legend.position = 'bottom',
                           strip.text.y = element_text(angle = 0),strip.background.y = element_rect(fill='grey20'))
p1

# p2 <- ggplot(all.1st2nd.exp %>% mutate(avg1DPE = ifelse(is.na(avg1DPE),0,avg1DPE)),aes(x=(avg1DPE+1),y=tfbs)) +
#   geom_bar(stat='identity',orientation = 'y',fill='firebrick4') +
#   geom_point(data=subset(all.1st2nd.exp,is.na(avg1DPE)),aes(x=3,y=tfbs),pch=4,color='red') +
#   geom_point(data=subset(all.1st2nd.exp,tfbs %in% tf.de$id),aes(x=avg1DPE*3,y=tfbs),pch=8,color='black') +
#   facet_grid(rows=vars(family),scales='free_y',space = 'free_y')+
#   scale_x_log10(expand=c(0,0)) +
#   xlab('Avg. 1DPE Expression\n(Log10 Axis)') +
#   theme_classic() + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),strip.background = element_blank(),strip.text = element_blank(),legend.position = 'bottom')
# 
# p3 <- p1 + p2 + plot_layout(widths = c(2,1))
# 
# p3

# plotB <- ggplot(all.1st2nd,aes(x=start,xend=stop,y=tfbs,yend=tfbs,color=round,size=posScore)) +
#   geom_segment(data=all.1st2nd.splits,aes(x=split_start,xend=split_end,y=tfbs,yend=tfbs,color=round),lwd=1,alpha=0.5,show.legend = T) +
#   # geom_point(inherit.aes = F,aes(x=start-15,y=tfbs,fill=round,group=group),pch=21) +
#   geom_segment(show.legend = T) +
#   scale_size_continuous(range = c(2,4)) +
#   scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#   facet_wrap(~group,scale='free',ncol = 1) +
#   ylab('') +
#   xlab('Relative Position') +
#   theme_linedraw() + theme(panel.grid.minor.x = element_blank())
# plotB

# plotA + plotB + plot_layout(widths = c(1,4))


#### Other VGs

vg.info <- read_tsv('../../../data/venom_annotations/___PriorityVenomGenes_Sept2020/PriorityVenomGeneList_withNVPs_02.25.21.txt',col_names = F) %>% 
  select(txid=X6,gene_id = X7) %>% 
  filter(str_detect(gene_id,'ADAM',negate = T))

other.tfbs.enrich <- read_tsv('_tfbs_enrichment_binding/Other_VGs/otherVG.enriched.genCoords.VG3VGUscores.relPos.bed',col_names = F) %>%
  mutate(highScoring = ifelse(X3 %in% highPosScores$tfbs,T,F))
# 
# other.tfbs.enrich %>% filter(highScoring)

other.tfbs.all <- read_tsv('_tfbs_enrichment_binding/Other_VGs/otherVG.allCiiider.genCoords.VG3VGUscores.relPos.bed',col_names = F) %>% 
  mutate(highScoring = ifelse(X3 %in% highPosScores$tfbs,T,F)) %>% 
  mutate(enriched = ifelse(X3 %in% other.tfbs.enrich$X3,T,F)) %>% 
  mutate(txid = str_split_fixed(X1,'::',2)[,1]) %>% 
  left_join(vg.info)


other.tfbs.all %>% filter(highScoring)



# ggplot(other.tfbs.enrich,aes(xmin=X5,xmax=X6,ymin=0,ymax=X7*X11, fill=X3)) +
#   geom_rect(size=10,alpha=0.5) +
#   facet_wrap(~X1,ncol=1,scales='free_x') +
#   theme_linedraw() 
# 
# ggplot(other.tfbs.enrich %>% filter(highScoring),aes(xmin=X5,xmax=X6,ymin=0,ymax=X7*X11, fill=X3)) +
#   geom_rect(size=10,alpha=0.5) +
#   facet_wrap(~X1,ncol=1,scales='free_x') +
#   theme_linedraw() 


ggplot(other.tfbs.all,aes(xmin=X5,xmax=X6,ymin=0,ymax=X7*X11, fill=X3)) +
  geom_rect(size=10,alpha=0.5) +
  facet_wrap(~gene_id,ncol=1,scales='free_x') +
  theme_linedraw() 

ggplot(other.tfbs.all %>% filter(highScoring),aes(xmin=X5,xmax=X6,ymin=0,ymax=X7*X11, fill=X3)) +
  geom_rect(size=10,alpha=0.5) +
  geom_hline(yintercept = 0,color='grey60',alpha=0.8) +
  geom_text(data=other.tfbs.all %>% filter(enriched & X11 > 0 & highScoring),aes(x=(X5+X6)/2,y=X7*X11*0.5,label="*")) +
  facet_wrap(~gene_id,ncol=1,scales='free_x',strip.position = 'left') +
  ylab('Venom Gene (Other)') +
  xlab('Position') +
  theme_linedraw() + theme(panel.grid.minor.y = element_blank())



other.tfbs.all.avgFootprint <- other.tfbs.all %>% 
  filter(highScoring) %>% 
  group_by(gene_id,X3) %>% 
  mutate(avgFootprint = mean(X11)) %>% 
  mutate(maxFootprint = max(X11)) %>% 
  select(gene_id,X3,avgFootprint,maxFootprint) %>% 
  mutate(round = 'Secondary') %>% 
  left_join(tf_fams,by=c('X3'='name')) %>% 
  group_by(family) %>% 
  mutate(order = length(unique(gene_id))) %>% 
  arrange(-order) %>% 
  mutate(family = factor(family,levels=unique(.$family))) 

other.all.1st2nd.exp <- other.tfbs.all.avgFootprint %>% 
  mutate(tfbs = str_to_upper(X3)) %>% 
  mutate(tfbs2 = tfbs) %>% 
  mutate(tfbs2 = ifelse(str_detect(tfbs2,'\\:'),str_split_fixed(tfbs2,'\\:\\:',2),tfbs2)) %>% 
  left_join(tf.exp,by=c('tfbs2'='id')) %>% 
  select(family,tfbs,avg1DPE) %>% 
  unique()




# ggplot(other.tfbs.all.avgFootprint,aes(x=gene_id,y=reorder(X3,avgFootprint),size=avgFootprint,fill=avgFootprint,shape=round)) +
#   geom_point() +
#   scale_size_continuous(range = c(2,8)) +
#   scale_shape_manual(values=c(23)) +
#   # scale_fill_discrete_sequential(palette = "Blues", nmax = 6,order = c(2,6)) +
#   scale_fill_continuous_sequential('Greens') +
#   facet_wrap(~cluster,ncol=1,strip.position = 'right',scales='free_y') +
#   ylab('TFBS') +
#   xlab('Venom Gene') +
#   theme_linedraw() + theme(panel.grid.minor.x = element_blank(),legend.position = 'right',axis.text.x = element_text(angle=45,vjust=1,hjust=1))


o.p1 <- ggplot(other.tfbs.all.avgFootprint,aes(x=gene_id,y=X3,size=maxFootprint,fill=maxFootprint,shape=round)) +
  geom_point() +
  scale_size_continuous(range = c(2,8)) +
  scale_shape_manual(values=c(23)) +
  # scale_fill_discrete_sequential(palette = "Blues", nmax = 6,order = c(2,6)) +
  scale_fill_continuous_sequential('Greens') +
  facet_grid(rows=vars(family),scales='free_y',space = 'free_y')+
  ylab('TFBS') +
  xlab('Venom Gene') +
  theme_linedraw(base_size = 16) + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1),panel.grid.minor.x = element_blank(),legend.position = 'bottom',
                           strip.text.y = element_text(angle = 0),strip.background.y = element_rect(fill='grey20'))
o.p1


# o.p2 <- ggplot(other.all.1st2nd.exp %>% mutate(avg1DPE = ifelse(is.na(avg1DPE),0,avg1DPE)),aes(x=(avg1DPE+1),y=tfbs)) +
#   geom_bar(stat='identity',orientation = 'y',fill='firebrick4') +
#   geom_point(data=subset(other.all.1st2nd.exp,is.na(avg1DPE)),aes(x=3,y=tfbs),pch=4,color='red') +
#   geom_point(data=subset(other.all.1st2nd.exp,tfbs %in% tf.de$id),aes(x=avg1DPE*3,y=tfbs),pch=8,color='black') +
#   facet_grid(rows=vars(cluster),scales='free_y',space = 'free_y')+
#   scale_x_log10(expand=c(0,0)) +
#   xlab('Avg. 1DPE Expression\n(Log10 Axis)') +
#   theme_classic() + theme(axis.title.y = element_blank(), axis.text.y=element_blank(), strip.background = element_blank(),strip.text = element_blank(),legend.position = 'bottom')
# 
# o.p1 + o.p2 + plot_layout(widths = c(3,1))
# 
# op3 <- o.p1 + o.p2 + plot_layout(widths = c(3,1))

p1 | o.p1
