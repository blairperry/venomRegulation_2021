library(gggenes)
library(viridis)
library(patchwork)
library(tidyverse)
library(scales)
library(ggforce)


svmp.blastRes <- read_tsv('analysis/8_NVP_Analyses/enhancers/coreVPERs_vs_genome.blastResults.08.25.21.txt',col_names = F) %>% 
  filter(str_detect(X1,'SVMP')) %>%  mutate(group = 'SVMP')
svsp.blastRes <- read_tsv('analysis/8_NVP_Analyses/enhancers/coreVPERs_vs_genome.blastResults.08.25.21.txt',col_names = F) %>% 
  filter(str_detect(X1,'SVSP')) %>%  mutate(group = 'SVSP')
pla2.blastRes <- read_tsv('analysis/8_NVP_Analyses/enhancers/coreVPERs_vs_genome.blastResults.08.25.21.txt',col_names = F) %>% 
  filter(str_detect(X1,'PLA2')) %>%  mutate(group = 'PLA2')

svsp.blastRes.tally <- svsp.blastRes %>% 
  bind_rows(svmp.blastRes,pla2.blastRes) %>% 
  filter(str_detect(X2,'un',negate = T)) %>% 
  group_by(group,X2) %>% 
  mutate(X2 = str_remove(X2,'scaffold-')) %>% 
  tally() 


ggplot(svsp.blastRes.tally,aes(y=reorder(X2,desc(X2)),x=n,label=n)) +
  geom_bar(stat='identity',orientation = 'y',alpha=0.5) +
  geom_text(aes(x=250)) +
  facet_wrap(~group) +
  labs(y='Chromosome',x='# Blast Hits') +
  theme_linedraw() + theme(panel.grid = element_blank())


##########

pri_venom_genes <- read_tsv('data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
exp <- read_tsv('analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
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
  mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>% 
  left_join(exp,by=c('tx_id'='txid')) %>% 
  mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>% 
  mutate(prom_start = ifelse(strand=='forward',start,end))



# SVMP  -------------------------------------------------------------------

SVMP.info <- all_info %>% 
  filter(str_detect(gene,'SVMP|ADAM')) 

SVMP.reg.start <- min(SVMP.info$start)-100000
SVMP.reg.end <- max(SVMP.info$end)+100000

SVMP.reg.length <- paste(c(round((SVMP.reg.end-SVMP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k27ac chipseq 
svmp.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi1.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

svmp.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

# Read in vPERs and super-enhancers
svmp.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(SVMP.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

svmp.ses <- read_tsv('analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 


# Read in ATACseq
svmp.atac.vg3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi1.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

svmp.atac.vg3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 


# Plotting


p.svmp.genes <-  ggplot(SVMP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVMP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'ADAM')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = -0.5,size=3,color='grey60') +
  
  geom_segment(aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=svmp.vpers,aes(x=prom_start,xend=start,y='mi1',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=svmp.vpers,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svmp.vpers, aes(x=(start+end)/2, y=type),size=2) +

  geom_segment(inherit.aes = F,data=svmp.vpers,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y= ' BLAST Hits',yend=' BLAST Hits'),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svmp.blastRes,aes(x=(X9+X10)/2,y=' BLAST Hits'),pch=5,size=3,color='red') +
  
  ylab('') +
  xlab('') +
  ggtitle('a) SVMP\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())

p.svmp.genes


p.svmp.h3k27ac.chip <- ggplot(svmp.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svmp.h3k27ac.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.h3k27ac.chip$density)*1.2)) +
  ggtitle('Chip-seq (H3K27ac)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svmp.atac.vg3 <- ggplot(svmp.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.atac.vg3.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svmp.atac.vg3$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.atac.vg3$density)*1.2)) +
  ggtitle('ATAC-seq (Post-Extraction Venom Gland)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svmp.scale <- ggplot(SVMP.info,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=SVMP.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=SVMP.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(SVMP.reg.start,SVMP.reg.end)),label=SVMP.reg.length)) +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_void()

svmp.fullPlots <- p.svmp.genes  / p.svmp.h3k27ac.chip / p.svmp.atac.vg3 + p.svmp.scale + plot_layout(ncol = 1,heights = c(2,2,2,1))
svmp.fullPlots



# SVSP  -------------------------------------------------------------------

svsp.info <- all_info %>% 
  filter(str_detect(gene,'SVSP')) 

svsp.reg.start <- min(svsp.info$start)-200000
svsp.reg.end <- max(svsp.info$end)+200000

svsp.reg.length <- paste(c(round((svsp.reg.end-svsp.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k27ac chipseq 
svsp.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi2.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

svsp.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

# Read in vPERs and super-enhancers

svsp.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(svsp.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

svsp.ses <- read_tsv('analysis/2_super_enhancers/bed//SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 


# Read in ATACseq
svsp.atac.vg3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi2.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

svsp.atac.vg3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 



# Plotting


p.svsp.genes <-  ggplot(svsp.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVSP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_segment(aes(x=svsp.reg.start,xend=svsp.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=svsp.vpers,aes(x=prom_start,xend=start,y='mi2',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=svsp.vpers,aes(x=svsp.reg.start,xend=svsp.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svsp.vpers, aes(x=(start+end)/2, y=type),size=2) +
  
  
  geom_segment(inherit.aes = F,data=svsp.vpers,aes(x=svsp.reg.start,xend=svsp.reg.end,y= ' BLAST Hits',yend=' BLAST Hits'),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svsp.blastRes,aes(x=(X9+X10)/2,y=' BLAST Hits'),pch=5,size=3,color='red') +
  
  ylab('') +
  xlab('') +
  ggtitle('b) SVSP\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())
p.svsp.genes


p.svsp.h3k27ac.chip <- ggplot(svsp.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.h3k27ac.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  geom_segment(inherit.aes = F,
               data=svsp.ses,
               aes(x=start,xend=end,y=max(svsp.h3k27ac.chip$density)*1.25,yend=max(svsp.h3k27ac.chip$density)*1.25),
               lwd=2,
               color='goldenrod2') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.h3k27ac.chip$density)*1.3)) +
  ggtitle('Chip-seq (H3K27ac) - Super-Enhancers Shown in Yellow') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.svsp.h3k27ac.chip

p.svsp.atac.vg3 <- ggplot(svsp.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.atac.vg3.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.atac.vg3$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.atac.vg3$density)*1.2)) +
  ggtitle('ATAC-seq (Post-Extraction Venom Gland)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svsp.scale <- ggplot(svsp.info,aes(x=svsp.reg.start,xend=svsp.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=svsp.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=svsp.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(svsp.reg.start,svsp.reg.end)),label=svsp.reg.length)) +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  theme_void()

svsp.fullPlots <- p.svsp.genes / p.svsp.h3k27ac.chip / p.svsp.atac.vg3 + p.svsp.scale + plot_layout(ncol = 1,heights = c(2,2,2,1))
svsp.fullPlots





# PLA2  ------------------------------------------------------------------- <- 

pla2.info <- all_info %>% 
  filter(str_detect(gene,'PLA2')) 

pla2.reg.start <- 3001808-10000
pla2.reg.end <- 3050307+10000
pla2.reg.length <- paste(c(round((pla2.reg.end-pla2.reg.start)/1000,digits = 2),'kb'),collapse = ' ')



# Read in  H3k27ac chipseq 
pla2.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi7.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

pla2.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

# Read in vPERs and super-enhancers

pla2.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/PLA2_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(pla2.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

pla2.ses <- read_tsv('analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 


# Read in ATACseq
pla2.atac.vg3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi7.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

pla2.atac.vg3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 



# Plotting


p.pla2.genes <-  ggplot(pla2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_segment(aes(x=pla2.reg.start,xend=pla2.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=pla2.vpers,aes(x=prom_start,xend=start,y='mi7',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=pla2.vpers,aes(x=pla2.reg.start,xend=pla2.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=pla2.vpers, aes(x=(start+end)/2, y=type),size=2) +
  
  geom_segment(inherit.aes = F,data=pla2.vpers,aes(x=pla2.reg.start,xend=pla2.reg.end,y= ' BLAST Hits',yend=' BLAST Hits'),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=pla2.blastRes,aes(x=(X9+X10)/2,y=' BLAST Hits'),pch=5,size=3,color='red') +
  
  ylab('') +
  xlab('') +
  ggtitle('c) PLA2\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())
p.pla2.genes



p.pla2.h3k27ac.chip <- ggplot(pla2.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = pla2.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.h3k27ac.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  geom_segment(inherit.aes = F,
               data=pla2.ses,
               aes(x=start,xend=end,y=max(pla2.h3k27ac.chip$density)*1.25,yend=max(pla2.h3k27ac.chip$density)*1.25),
               lwd=2,
               color='goldenrod2') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.h3k27ac.chip$density)*1.3)) +
  ggtitle('Chip-seq (H3K27ac) - Super-Enhancers Shown in Yellow') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.pla2.h3k27ac.chip

p.pla2.atac.vg3 <- ggplot(pla2.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = pla2.atac.vg3.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.atac.vg3$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.atac.vg3$density)*1.2)) +
  ggtitle('ATAC-seq (Post-Extraction Venom Gland)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.pla2.scale <- ggplot(pla2.info,aes(x=pla2.reg.start,xend=pla2.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=pla2.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=pla2.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(pla2.reg.start,pla2.reg.end)),label=pla2.reg.length)) +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  theme_void()

pla2.fullPlots <- p.pla2.genes  / p.pla2.h3k27ac.chip / p.pla2.atac.vg3 + p.pla2.scale + plot_layout(ncol = 1,heights = c(2,2,2,1))
pla2.fullPlots


svmp.fullPlots 
svsp.fullPlots 
pla2.fullPlots
