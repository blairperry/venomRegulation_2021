
library(tidyverse)
library(patchwork)
library(ggridges)
library(ggforce)
library(gggenes)



# Plot of number elements per chromosome ----------------------------------

all.tes <- read_tsv('analysis/11_TEs/bed/CroMit_hAT-Tip100_genomic.bed.txt',col_names = F)

chrom_sizes <- read_tsv('data/misc/Cvv_ChromSizes_OldID.txt',col_names = F) %>% 
  filter(str_detect(X1,'un',negate = T))

all.tes.summary <- all.tes %>% 
  group_by(X1) %>% 
  tally() %>% 
  left_join(chrom_sizes) %>% 
  mutate(norm = n/X2) %>% 
  mutate(chrom = str_remove(X1,'scaffold-')) %>% 
  mutate(chrom_type = ifelse(str_detect(chrom,'mi'),'Micro','Macro'))
 

ggplot(all.tes.summary,aes(y=reorder(chrom,X2),x=norm)) +
  # geom_bar(stat='identity',orientation = 'y') +
  geom_segment(inherit.aes = F,aes(x=0,xend=norm,y=reorder(chrom,X2),yend=reorder(chrom,X2)),lwd=1,color='grey')+
  geom_point(size=4,aes(color=chrom_type)) +
  labs(y='Chromosome',x='Number of Annotated Elements / Chrom Length',color='Chromosome Type') +
  scale_x_continuous(expand = c(0,0),limits=c(0,max(all.tes.summary$norm)*1.1)) +
  ggtitle('Number of Cv1-hAt-Tip100 Elements (Normalized by Chrom Length)') +
  theme_linedraw(base_size = 12) + theme(panel.grid = element_blank(), legend.position = c(0.85,0.85),plot.title.position = 'plot',plot.title=element_text(face='bold')) 



# Parse TE/SVSP region overlaps -------------------------------------------

te.enhOverlap <- read_tsv('analysis/11_TEs/bed/overlap/overlap_HatTipEnhancers.bed',col_names = F) %>% mutate(TE_ID = paste(paste(X1,X2,sep = ':'),X3,sep = '-')) %>% mutate(group = 'Enhancer') %>% 
  select(chr=X1,start=X2,end=X3,feature.start=X11,feature.end=X12,TE_ID,group)
te.promOverlap <- read_tsv('analysis/11_TEs/bed/overlap/overlap_HatTipPromoters.bed',col_names = F) %>% mutate(TE_ID = paste(paste(X1,X2,sep = ':'),X3,sep = '-')) %>% mutate(group = 'Promoter') %>% 
  select(chr=X1,start=X2,end=X3,feature.start=X11,feature.end=X12,TE_ID,group)
te.otherOverlap <- read_tsv('analysis/11_TEs/bed/overlap/overlap_HatTipWholeRegion.bed',col_names = F) %>% mutate(TE_ID = paste(paste(X1,X2,sep = ':'),X3,sep = '-')) %>% mutate(group = 'Other') %>% 
  select(chr=X1,start=X2,end=X3,feature.start=X11,feature.end=X12,TE_ID,group) %>% 
  mutate(seen = ifelse(TE_ID %in% te.enhOverlap$TE_ID,'yes',
                       ifelse(TE_ID %in% te.promOverlap$TE_ID,'yes','no'))) %>% 
  filter(seen == 'no') %>% 
  select(-seen)

te.combined <- te.enhOverlap %>% 
  bind_rows(te.promOverlap,te.otherOverlap) %>% 
  arrange(start)

# write_tsv(te.combined,'analysis/11_TEs/bed/overlap/SVSP_FunctionalRegionTEOverlaps_08.31.21.bed')



# Distribution of pairwise divergences to consensus -----------------------

mut_rate <- 2.4 * 10^9

focal_peaks <- read_tsv('analysis/11_TEs/Pairwise_PI_results/SVSP_FunctionalRegionTEOverlaps_wPI_08.31.21.txt',col_names = F) %>% 
  select(id=X6,type=X7)

all.pi <- read_tsv('analysis/11_TEs/Pairwise_PI_results//AllCopies_PIwithID.txt',col_names = c('id','pi')) %>% 
  mutate(chrom = str_split_fixed(str_remove(id,'scaffold-'),':',2)[,1]) %>% 
  left_join(focal_peaks) %>% 
  mutate(type = ifelse(is.na(type),'Genomic',type)) %>% 
  mutate(type = factor(type,levels=c('Promoter','Enhancer','Other','Genomic'))) %>% 
  mutate(time = pi / 2 * mut_rate)



# Pi on X
p.piDensity <- ggplot(all.pi,aes(x=pi)) +
  geom_density(fill='grey',alpha=0.5) +
  # geom_segment(inherit.aes = F, data=all.pi %>% filter(type != 'Genomic'), aes(x=pi,xend=pi,y=0,yend=23,color=type)) +
  # geom_point(inherit.aes = F, data=all.pi %>% filter(type != 'Genomic'), aes(x=pi,y=23,color=type)) +
  labs(x='Pairwise Pi',y='Density',title='Pairwise divergence from genome-wide consensus sequence') +
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  scale_x_continuous(limits=c(0,0.15)) +
  theme_linedraw(base_size = 12) + theme(plot.title = element_text(face='bold'),plot.title.position = 'plot',axis.title.x = element_blank())

p.focalPi <- ggplot(all.pi %>% filter(type != 'Genomic'),aes(x=pi,y=reorder(type,desc(type)),color=type)) +
  geom_point(show.legend = F,size=3,alpha=0.75) +
  scale_x_continuous(limits=c(0,0.15)) +
  labs(x='Pairwise Pi',y='Cv1-hAT-Tip100s\nin SVSP region') +
  scale_color_manual(values=c('Promoter'='#2C8790','Enhancer'='#82B47A','Other'='grey70','Genomic'='grey30')) +
  theme_linedraw(base_size = 12)

p.piDensity / p.focalPi + plot_layout(heights = c(5,1))


ggplot(all.pi,aes(x=pi,y=reorder(type,desc(type)),fill=reorder(type,desc(type)))) +
  geom_boxplot(outlier.alpha = 0,alpha=0.5,show.legend = F) +
  geom_jitter(height = 0.15,alpha=0.75,pch=1,show.legend = F) +
  labs(x='Pairwise Pi',y='',title='Pairwise divergence from genome-wide consensus sequence') +
  scale_fill_manual(values=c('Promoter'='#2C8790','Enhancer'='#82B47A','Other'='grey90','Genomic'='grey30')) +
  # scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  # scale_x_continuous(limits=c(0,0.15)) +
  theme_linedraw(base_size = 12) + theme(plot.title = element_text(face='bold'),plot.title.position = 'plot')

ggplot(all.pi,aes(x=time/1000000,y=reorder(type,desc(type)),fill=reorder(type,desc(type)))) +
  geom_boxplot(outlier.alpha = 0,alpha=0.5,show.legend = F) +
  geom_jitter(height = 0.15,alpha=0.75,pch=1,show.legend = F) +
  labs(x='Estimated age based on pairwise pi (MYA)',y='',title='Estimated age based on pairwise divergence from\ngenome-wide consensus sequence') +
  scale_fill_manual(values=c('Promoter'='#2C8790','Enhancer'='#82B47A','Other'='grey90','Genomic'='grey30')) +
  # scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  scale_x_continuous(labels = scales::comma) +
  theme_linedraw(base_size = 12) + theme(plot.title = element_text(face='bold'),plot.title.position = 'plot')



# Multitrack plotting of SVSP region w/ TEs -------------------------------


# SVSP  -------------------------------------------------------------------

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

svsp.info <- all_info %>% 
  filter(str_detect(tx_id,'SVSP')) %>% 
  mutate(molecule = str_remove(molecule,'scaffold-')) %>% 
  mutate(molecule = factor(molecule,levels=c('Promoter','mi2','vPERs','Enhancer','Other')))

svsp.reg.start <- min(svsp.info$start)-100000
svsp.reg.end <- max(svsp.info$end)+100000

svsp.reg.length <- paste(c(round((svsp.reg.end-svsp.reg.start)/1000,digits = 2),'kb'),collapse = ' ')


# Read in vPERs and super-enhancers

svsp.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(svsp.info,by=c('gene'='gene.x')) %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start) %>% 
  mutate(molecule = str_remove(molecule,'scaffold-')) %>% 
  mutate(molecule = factor(molecule,levels=c('Promoter','mi2','vPERs','Enhancer','Other')))


# Read in TE bedfile

svsp.tes.enhan <- read_tsv('analysis/11_TEs/bed/overlap/SVSP_FunctionalRegionTEOverlaps_08.31.21.bed',col_names = T) %>% 
  filter(group == 'Enhancer')

svsp.tes.promo <- read_tsv('analysis/11_TEs/bed/overlap/SVSP_FunctionalRegionTEOverlaps_08.31.21.bed',col_names = T) %>% 
  filter(group == 'Promoter')

svsp.tes.all <- read_tsv('analysis/11_TEs/bed/CroMit_hAT-Tip100_genomic.bed.txt',col_names = F) %>% 
  filter(X1 == 'scaffold-mi2') %>% 
  mutate(TE_ID = paste(paste(X1,X2,sep = ':'),X3,sep = '-')) %>% 
  mutate(type = ifelse(TE_ID %in% svsp.tes.enhan$TE_ID,'Enhancer','Other')) %>% 
  mutate(type = ifelse(TE_ID %in% svsp.tes.promo$TE_ID,'Promoter',type)) %>% 
  mutate(type = factor(type,levels=c('Promoter','mi2','vPERs','Enhancer','Other')))

svsp.focalTEs <- svsp.tes.all %>% 
  filter(X3 >= svsp.reg.start & X2 <= svsp.reg.end)

svsp.dummy <-  tribble(
  ~molecule, ~start,  ~end,
  "Promoter", svsp.reg.start, svsp.reg.end,
  "mi2", svsp.reg.start, svsp.reg.end,
  "vPERs", svsp.reg.start, svsp.reg.end,
  "Enhancer", svsp.reg.start, svsp.reg.end,
  "Other", svsp.reg.start, svsp.reg.end,
) %>% 
  mutate(molecule = factor(molecule,levels=rev(c('Promoter','mi2','vPERs','Enhancer','Other'))))



# Plotting

ggplot(svsp.dummy,aes(x = start, xend = end, y = molecule, yend=molecule)) +
  geom_segment(lwd=1,color='grey70') +
  geom_diagonal(data=svsp.vpers,aes(x=prom_start,xend=start,y='mi2',yend='vPERs',alpha = stat(index)),strength = -0.2,show.legend = F) +
  geom_gene_arrow(data=svsp.info,aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10avg1DPE),arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  geom_point(inherit.aes = F, data=svsp.vpers, aes(x=(start+end)/2, y='vPERs'),size=2) +
  geom_segment(inherit.aes = F, data=svsp.tes.all %>% filter(type=='Promoter'), aes(x=(X2+X3)/2,xend=(X2+X3)/2,y='Promoter',yend='mi2',color=type),alpha=0.5) +
  geom_segment(inherit.aes = F, data=svsp.tes.all %>% filter(type=='Enhancer'), aes(x=(X2+X3)/2,xend=(X2+X3)/2,y='Enhancer',yend='vPERs',color=type),alpha=0.5) +
  geom_point(inherit.aes = F, data=svsp.tes.all, aes(x=(X2+X3)/2, y=type,color=type),size=5,pch=18) +
  scale_color_manual(values = c('Promoter'='#2C8790','Enhancer'='#82B47A','Other'='grey40')) +
  scale_x_continuous(labels = scales::comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  ylab('') +
  xlab('') +
  ggtitle('b) SVSP\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  theme_classic(base_size = 14) + theme(panel.grid = element_blank(),axis.line = element_blank(),
                                        plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())


