
library(tidyverse)
library(patchwork)
library(ggforce)
library(gggenes)

# TE Composition ----------------------------------------------------------

# te.fullDetail <- read_tsv('~/Downloads/BED_extended/CroVir_GenomeLib_TE.bed',col_names = F) %>% 
#   select(X4) %>% 
#   unique() %>% 
#   separate(X4,sep='#',into = c('id','family'))

te.venomRegs <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/RepeatAnnotation.focalVenomRegs.txt', col_names = F) %>% 
  select(chr=1,te=7,family=8,start=5,end=6) %>% 
  mutate(venom_cluster = case_when(
    chr == 'scaffold-mi1' ~ 'SVMP',
    chr == 'scaffold-mi2' ~ 'SVSP',
    chr == 'scaffold-mi7' ~ 'PLA2',
  ))


te.composition <- te.venomRegs %>% 
  group_by(venom_cluster,family) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(venom_cluster) %>% 
  top_n(n = 10,wt=n) %>%  # Taking top 10 most abundant element types per group %>% 
  mutate(group = case_when(
    str_detect(family,'LTR') ~ 'LTR',
    str_detect(family,'INE') ~ 'Non-LTR',
    str_detect(family,'DNA') ~ 'DNA',
    str_detect(family,'rich') ~ 'Other'
  )) %>% 
  arrange(group,family)

te.composition.heatdata <- te.composition %>% 
  group_by(venom_cluster) %>% 
  mutate(venom_cluster = factor(venom_cluster,levels=c('SVMP','SVSP','PLA2'))) %>% 
  mutate(percent = n / sum(n)) %>% 
  select(-n,-group) %>% 
  pivot_wider(names_from = venom_cluster,values_from = percent,values_fill = NA) %>% 
  as.data.frame() %>% 
  select(family,SVMP,SVSP,PLA2) %>% 
  column_to_rownames('family')

te.composition.annot <- te.composition %>% 
  ungroup() %>% 
  select(family,group) %>% 
  unique() %>% 
  column_to_rownames('family')


pheatmap::pheatmap(te.composition.heatdata,
                   cluster_cols = F,
                   cluster_rows = F,
                   na_col = 'white',
                   cellwidth = 25,cellheight = 25,
                   color = viridis::viridis(50),
                   # annotation_row = te.composition.annot
                   )



# TE Abundance ------------------------------------------------------------

te.cov.10kb <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/repeatCoverage.10kb.bedgraph',col_names = c('chr','start','end','coverage'))
te.cov.100kb <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/repeatCoverage.100kb.bedgraph',col_names = c('chr','start','end','coverage'))

# Relevant regions (venom clusters)
#SVMP
# scaffold-mi1 13881006 14424729

#SVSP
# scaffold-mi2 8548728 8981362

#PLA2
# scaffold-mi7 3007531 3043778


te.cov.svmp <- te.cov.10kb %>% 
  filter(chr=='scaffold-mi1') %>% 
  mutate(venom_overlap = ifelse(
    (start > 13881006 & start < 14424729) |
    (end > 13881006 & end < 14424729),T,F)
  )

p.svmp <- ggplot(te.cov.svmp,aes(y=coverage,x=venom_overlap,fill=venom_overlap)) +
  geom_boxplot(show.legend = F) +
  labs(x='Overlap with venom region',y='% bases annotated as repeat',title='SVMP') +
  theme_linedraw()

te.cov.svsp <- te.cov.10kb %>% 
  filter(chr=='scaffold-mi2') %>% 
  mutate(venom_overlap = ifelse(
    (start > 8548728 & start < 8981362) |
      (end > 8548728 & end   < 8981362),T,F)
  )

p.svsp <- ggplot(te.cov.svsp,aes(y=coverage,x=venom_overlap,fill=venom_overlap)) +
  geom_boxplot(show.legend = F) +
  labs(x='Overlap with venom region',y='% bases annotated as repeat',title='SVSP') +
  theme_linedraw()

te.cov.pla2 <- te.cov.10kb %>% 
  filter(chr=='scaffold-mi7') %>% 
  mutate(venom_overlap = ifelse(
    (start > 3007531 & start < 3043778) |
      (end > 3007531 & end   < 3043778),T,F)
  )

p.pla2 <- ggplot(te.cov.pla2,aes(y=coverage,x=venom_overlap,fill=venom_overlap)) +
  geom_boxplot(show.legend = F) +
  labs(x='Overlap with venom region',y='% bases annotated as repeat',title='PLA2') +
  theme_linedraw()

p.svmp + p.svsp + p.pla2


# Plot and summarize Giggle results ---------------------------------------

colnames <- c('TE',	'genome_count',	'peak_count',	'odds_ratio',	'fishers_two_tail',	'fishers_left_tail',	'fishers_right_tail',	'combo_score','ignore')

svmp.giggle <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/giggle_rerun/SVMP.PromotersEnhancers.bed.giggleStats.sort',col_names = colnames) %>% 
  select(-ignore) %>% 
  mutate(TE = str_remove_all(TE,'/'))

pla2.giggle <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/giggle_rerun/PLA2.PromotersEnhancers.bed.giggleStats.sort',col_names = colnames) %>% 
  select(-ignore)%>% 
  mutate(TE = str_remove_all(TE,'/'))

svmp.giggle.sig <- svmp.giggle %>% filter(fishers_right_tail<0.05)
pla2.giggle.sig <- pla2.giggle %>% filter(fishers_right_tail<0.05)

svmp.teOverlaps <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/SVMP.TEoverlap.txt',col_names = F) %>% 
  filter(X5 != '.') %>% 
  select(chr=1,reg_start=2,reg_end=3,reg_id=4,te_start=6,te_end=7,te_id=8,te_family=9,overlap_bases=11) %>% 
  filter(te_id %in% svmp.giggle.sig$TE)

# Note: neither overlap with SVMP vPERs that have bound TFBS

pla2.teOverlaps <- read_tsv('_ms/__GR_Reviews_R1/TE_revisionAnalyses/PLA2.TEoverlap.txt',col_names = F) %>% 
  filter(X5 != '.') %>% 
  select(chr=1,reg_start=2,reg_end=3,reg_id=4,te_start=6,te_end=7,te_id=8,te_family=9,overlap_bases=11) %>% 
  filter(te_id %in% pla2.giggle.sig$TE)  %>% 
  filter(str_detect(reg_id,'PLA2_gIIE',negate = T))

# Agkistrodon_rnd-1_family-415 DNA/TcMar-Tigger overlaps with PLA2 A and B promoters
# Overlaps with PER36:
#   - Cmitchellii_rnd-5_family-57 LINE/CR1
#   - ERV1-4_SSc-I LTR/ERV1
# Overlaps with PER38:
#   - Thamnophis454_rnd-1_family-109 DNA/TcMar (only by 9 bp)
# A Catrox_rnd-3_family-1141 overlaps with per39, but that per has no bound tfbs



# Supp MultiTrack TE plots ------------------------------------------------

# svmp  -------------------------------------------------------------------

pri_venom_genes <- read_tsv('data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)

exp <- read_tsv('analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  mutate(gene = str_replace_all(txid,'_',' ')) %>% 
  mutate(log10avg1DPE = log10(Median1DPE + 1))


all_info <- read_tsv('data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
  filter(str_detect(X9,'trnascan',negate = T)) %>% 
  mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>% 
  mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>% 
  filter(tx_id %in% pri_venom_genes$X6) %>% 
  left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>% 
  select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>% 
  mutate(strand = ifelse(strand == '+','forward','reverse')) %>% 
  mutate(direction = ifelse(strand == 'forward',1,-1)) %>% 
  mutate(gene = str_replace_all(gene,'_',' ')) %>% 
  left_join(exp,by=c('tx_id'='txid')) %>% 
  mutate(prom_start = ifelse(strand=='forward',start,end))

svmp.info <- all_info %>% 
  filter(str_detect(tx_id,'SVMP')) %>% 
  mutate(molecule = str_remove(molecule,'scaffold-')) %>% 
  mutate(molecule = factor(molecule,levels=c('Promoter','mi1','vPERs','Enhancer','Other')))

svmp.reg.start <- min(svmp.info$start)-100000
svmp.reg.end <- max(svmp.info$end)+100000

svmp.reg.length <- paste(c(round((svmp.reg.end-svmp.reg.start)/1000,digits = 2),'kb'),collapse = ' ')


# Read in vPERs and super-enhancers

svmp.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(svmp.info,by=c('gene'='gene.x')) %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start) %>% 
  mutate(molecule = str_remove(molecule,'scaffold-')) %>% 
  mutate(molecule = factor(molecule,levels=c('Promoter','mi1','vPERs','Enhancer','Other')))


# Read in TE bedfile


svmp.dummy <-  tribble(
  ~molecule, ~start,  ~end,
  "Promoter", svmp.reg.start, svmp.reg.end,
  "mi1", svmp.reg.start, svmp.reg.end,
  "vPERs", svmp.reg.start, svmp.reg.end,
  "Enhancer", svmp.reg.start, svmp.reg.end,
) %>% 
  mutate(molecule = factor(molecule,levels=rev(c('Promoter','mi1','vPERs','Enhancer'))))



# Plotting

ggplot(svmp.dummy,aes(x = start, xend = end, y = molecule, yend=molecule)) +
  geom_segment(lwd=1,color='grey70') +
  geom_diagonal(data=svmp.vpers,aes(x=prom_start,xend=start,y='mi1',yend='vPERs',alpha = stat(index)),strength = -0.2,show.legend = F) +
  geom_gene_arrow(data=svmp.info,aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10avg1DPE),arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  geom_point(inherit.aes = F, data=svmp.vpers, aes(x=(start+end)/2, y='vPERs'),size=2) +
  # geom_segment(inherit.aes = F, data=svmp.tes.all %>% filter(type=='Promoter'), aes(x=(X2+X3)/2,xend=(X2+X3)/2,y='Promoter',yend='mi2',color=type),alpha=0.5) +
  geom_segment(inherit.aes = F, data=svmp.teOverlaps, aes(x=(te_start+te_end)/2,xend=(te_start+te_end)/2,y='Enhancer',yend='vPERs',color=te_id),alpha=0.5) +
  geom_point(inherit.aes = F, data=svmp.teOverlaps, aes(x=(te_start+te_end)/2, y='Enhancer',color=te_id),size=5,pch=18) +
  # scale_color_manual(values = c('Promoter'='#2C8790','Enhancer'='#82B47A','Other'='grey40')) +
  scale_x_continuous(labels = scales::comma,limits=c(svmp.reg.start,svmp.reg.end),expand=c(0,0)) +
  ylab('') +
  xlab('') +
  ggtitle('SVMP Tandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  theme_classic(base_size = 14) + theme(panel.grid = element_blank(),axis.line = element_blank(),
                                        plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())

#



# pla2  -------------------------------------------------------------------


pla2.info <- all_info %>% 
  filter(str_detect(tx_id,'PLA2')) %>% 
  mutate(molecule = str_remove(molecule,'scaffold-')) %>% 
  mutate(molecule = factor(molecule,levels=c('Promoter','mi7','vPERs','Enhancer','Other')))

pla2.reg.start <- min(pla2.info$start)-5000
pla2.reg.end <- max(pla2.info$end)+5000

pla2.reg.length <- paste(c(round((pla2.reg.end-pla2.reg.start)/1000,digits = 2),'kb'),collapse = ' ')


# Read in vPERs and super-enhancers

pla2.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/pla2_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(pla2.info,by=c('gene'='gene.x')) %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start) %>% 
  mutate(molecule = str_remove(molecule,'scaffold-')) %>% 
  mutate(molecule = factor(molecule,levels=c('Promoter','mi7','vPERs','Enhancer','Other')))


# Read in TE bedfile


pla2.dummy <-  tribble(
  ~molecule, ~start,  ~end,
  "Promoter", pla2.reg.start, pla2.reg.end,
  "mi7", pla2.reg.start, pla2.reg.end,
  "vPERs", pla2.reg.start, pla2.reg.end,
  "Enhancer", pla2.reg.start, pla2.reg.end,
) %>% 
  mutate(molecule = factor(molecule,levels=rev(c('Promoter','mi7','vPERs','Enhancer'))))



# Plotting

ggplot(pla2.dummy,aes(x = start, xend = end, y = molecule, yend=molecule)) +
  geom_segment(lwd=1,color='grey70') +
  geom_diagonal(data=pla2.vpers,aes(x=prom_start,xend=start,y='mi7',yend='vPERs',alpha = stat(index)),strength = -0.2,show.legend = F) +
  geom_gene_arrow(data=pla2.info,aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10avg1DPE),arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  geom_point(inherit.aes = F, data=pla2.vpers, aes(x=(start+end)/2, y='vPERs'),size=2) +
  geom_segment(inherit.aes = F, data=pla2.teOverlaps %>% filter(str_detect(reg_id,'PROMOTER')), aes(x=(te_start+te_end)/2,xend=(te_start+te_end)/2,y='Promoter',yend='mi7',color=te_id),alpha=0.5) +
  geom_segment(inherit.aes = F, data=pla2.teOverlaps %>% filter(str_detect(reg_id,'PER')), aes(x=(te_start+te_end)/2,xend=(te_start+te_end)/2,y='Enhancer',yend='vPERs',color=te_id),alpha=0.5) +
  geom_point(inherit.aes = F, data=pla2.teOverlaps %>% filter(str_detect(reg_id,'PER')), aes(x=(te_start+te_end)/2, y='Enhancer',color=te_id),size=5,pch=18) +
  geom_point(inherit.aes = F, data=pla2.teOverlaps %>% filter(str_detect(reg_id,'PROMOTER')), aes(x=(te_start+te_end)/2, y='Promoter',color=te_id),size=5,pch=18) +
  # scale_color_manual(values = c('Promoter'='#2C8790','Enhancer'='#82B47A','Other'='grey40')) +
  scale_x_continuous(labels = scales::comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  ylab('') +
  xlab('') +
  ggtitle('PLA2 Tandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  theme_classic(base_size = 14) + theme(panel.grid = element_blank(),axis.line = element_blank(),
                                        plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())
#




  




