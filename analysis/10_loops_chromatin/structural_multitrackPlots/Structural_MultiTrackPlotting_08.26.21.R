
#install.packages('gggenes')
library(gggenes)
library(viridis)
library(patchwork)
library(tidyverse)
library(scales)
library(ggforce)
library(colorspace)


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

SVMP.reg.start <- min(SVMP.info$start)-30000
SVMP.reg.end <- max(SVMP.info$end)+30000

SVMP.reg.length <- paste(c(round((SVMP.reg.end-SVMP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k4me3 chipseq
svmp.h3k4me3.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi1.bed',
                      col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) %>% 
  mutate(type = 'H3K4me3')

svmp.h3k4me3.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) %>% 
  mutate(type = 'H3K4me3')

svmp.h3k4me3.regionMax <- svmp.h3k4me3.chip %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) %>% 
  summarise(max=max(density))

# Read in  H3k27ac chipseq 
svmp.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi1.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) %>% 
  mutate(type = 'H3K27ac')

svmp.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) %>% 
  mutate(type = 'H3K27ac')

svmp.h3k27ac.regionMax <- svmp.h3k27ac.chip %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) %>% 
  summarise(max=max(density))

# Read in  CTCF chipseq and related files
svmp.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/CTCF_CroVir_i87_uniq_signal_mi1.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) 

svmp.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) 

svmp.ctcf.bound <- read_tsv('analysis/10_loops_chromatin/contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) 

svmp.ctcf.regionMax <- svmp.ctcf.chip %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) %>% 
  summarise(max=max(density))


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
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) %>% 
  mutate(type='Super-Enhancers')


# Read in ATACseq
svmp.atac.mean <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi1.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1)  %>% 
  mutate(type = 'ATAC-seq (Mean)')

svmp.atac.mean.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                          col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start*0.9 & start <= SVMP.reg.end*1.1) %>% 
  mutate(type = 'ATAC-seq (Mean)')

svmp.atac.mean.regionMax <- svmp.atac.mean %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) %>% 
  summarise(max=max(density))


# Read in loops and TADs

svmp.CTCFloops <- read_tsv('analysis/10_loops_chromatin/contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>% 
  select(-7,-8) %>% 
  filter(chrA=='scaffold-mi1') 

svmp.loops <- read_tsv('analysis/10_loops_chromatin/contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>% 
  select(-7,-8) %>% 
  filter(chrA=='scaffold-mi1') 

svmp.TADs <- read_tsv('analysis/10_loops_chromatin/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) %>% 
  rowid_to_column() %>% 
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


p.svmp.genes <-  ggplot(SVMP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVMP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'ADAM')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = -0.25,size=3) +
  
  geom_segment(aes(x=SVMP.reg.start*0.9,xend=SVMP.reg.end*1.1,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=svmp.vpers,aes(x=prom_start,xend=start,y='mi1',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=svmp.vpers,aes(x=SVMP.reg.start*0.9,xend=SVMP.reg.end*1.1,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svmp.vpers, aes(x=(start+end)/2, y=type),size=2) +
  ylab('') +
  xlab('') +
  ggtitle('Tandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(SVMP.reg.start,SVMP.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())


# TADS plot
p.svmp.tads <- ggplot(svmp.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(SVMP.reg.start,SVMP.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.svmp.loops <- ggplot(svmp.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  coord_cartesian(xlim=c(SVMP.reg.start,SVMP.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.svmp.h3k4me3.chip <- ggplot(svmp.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.h3k4me3.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svmp.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  coord_cartesian(expand = T,xlim=c(SVMP.reg.start,SVMP.reg.end),ylim=c(0,svmp.h3k4me3.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K4me3)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svmp.h3k27ac.chip <- ggplot(svmp.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=svmp.h3k27ac.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  coord_cartesian(expand = T,xlim=c(SVMP.reg.start,SVMP.reg.end),ylim=c(0,svmp.h3k27ac.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K27ac)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svmp.ctcf.chip <- ggplot(svmp.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svmp.ctcf.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=svmp.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = svmp.ctcf.chip.peaks,aes(x=start-200,xend=end+200,y=svmp.ctcf.regionMax$max*1.2,yend=svmp.ctcf.regionMax$max*1.2),lwd=5,alpha=1) +
  geom_point(inherit.aes = F,data=svmp.ctcf.bound,aes(x=(start+end)/2),y=svmp.ctcf.regionMax$max*1.1,pch=25,size=4,fill='purple4',alpha=0.75) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(SVMP.reg.start,SVMP.reg.end),ylim=c(0,svmp.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.svmp.ctcf.chip

p.svmp.atac.mean <- ggplot(svmp.atac.mean,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.atac.mean.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=svmp.atac.mean.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  coord_cartesian(expand = T,xlim=c(SVMP.reg.start,SVMP.reg.end),ylim=c(0,svmp.atac.mean.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('ATAC-seq (Post-Extraction Venom Gland)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svmp.scale <- ggplot(SVMP.info,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=SVMP.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=SVMP.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(SVMP.reg.start,SVMP.reg.end)),label=SVMP.reg.length)) +
  coord_cartesian(expand = T,xlim=c(SVMP.reg.start,SVMP.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

# svmp.fullPlots <- p.svmp.tads / p.svmp.ctcf.chip / p.svmp.loops / p.svmp.genes  /  p.svmp.h3k27ac.chip / p.svmp.h3k4me3.chip/ p.svmp.atac.mean + p.svmp.scale + plot_layout(ncol = 1,heights = c(1,2,2,2,2,2,2,1))
# svmp.fullPlots

svmp.all.peaks <- svmp.h3k27ac.chip.peaks %>% 
  bind_rows(svmp.h3k4me3.chip.peaks,svmp.atac.mean.peaks)


p.svmp.all.peaks <- ggplot(svmp.all.peaks,aes(x=start,xend=end, y=type,yend=type,color=type)) +
  geom_point(show.legend = F) +
  scale_color_manual(values = c('H3K4me3'='#2C8790','H3K27ac'='#82B47A','ATAC-seq (Mean)'='#824029','Super-Enhancers'='goldenrod2')) +
  coord_cartesian(expand = T,xlim=c(SVMP.reg.start,SVMP.reg.end)) +
  scale_x_continuous(labels=comma) +
  xlab('Position') +
  theme_linedraw(base_size = 14)+ theme(axis.title.y = element_blank())


p.svmp.condensed <- p.svmp.tads / p.svmp.ctcf.chip / p.svmp.loops / p.svmp.genes  / p.svmp.all.peaks + plot_layout(heights = c(1,2,2,2,2))
p.svmp.condensed


# SVSP  -------------------------------------------------------------------

SVSP.info <- all_info %>% 
  filter(str_detect(gene,'SVSP')) 

SVSP.reg.start <- min(SVSP.info$start)-50000
SVSP.reg.end <- max(SVSP.info$end)+50000

SVSP.reg.length <- paste(c(round((SVSP.reg.end-SVSP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k4me3 chipseq
svsp.h3k4me3.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi2.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) %>% 
  mutate(type = 'H3K4me3')

svsp.h3k4me3.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) %>% 
  mutate(type = 'H3K4me3')

svsp.h3k4me3.regionMax <- svsp.h3k4me3.chip %>% 
  filter(end >= SVSP.reg.start & start <= SVSP.reg.end) %>% 
  summarise(max=max(density))

# Read in  H3k27ac chipseq 
svsp.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi2.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) %>% 
  mutate(type = 'H3K27ac')

svsp.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) %>% 
  mutate(type = 'H3K27ac')

svsp.h3k27ac.regionMax <- svsp.h3k27ac.chip %>% 
  filter(end >= SVSP.reg.start & start <= SVSP.reg.end) %>% 
  summarise(max=max(density))

# Read in  CTCF chipseq and related files
svsp.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/CTCF_CroVir_i87_uniq_signal_mi2.bed',
                           col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) 

svsp.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                 col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) 

svsp.ctcf.bound <- read_tsv('analysis/10_loops_chromatin/contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) 

svsp.ctcf.regionMax <- svsp.ctcf.chip %>% 
  filter(end >= SVSP.reg.start & start <= SVSP.reg.end) %>% 
  summarise(max=max(density))


# Read in vPERs and super-enhancers
svsp.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(SVSP.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

svsp.ses <- read_tsv('analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) %>% 
  mutate(type = 'Super-Enhancers')


# Read in ATACseq
svsp.atac.mean <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi2.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1)  %>% 
  mutate(type = 'ATAC-seq (Mean)')

svsp.atac.mean.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start*0.9 & start <= SVSP.reg.end*1.1) %>% 
  mutate(type = 'ATAC-seq (Mean)')

svsp.atac.mean.regionMax <- svsp.atac.mean %>% 
  filter(end >= SVSP.reg.start & start <= SVSP.reg.end) %>% 
  summarise(max=max(density))


# Read in loops and TADs

svsp.CTCFloops <- read_tsv('analysis/10_loops_chromatin/contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>% 
  select(-7,-8) %>% 
  filter(chrA=='scaffold-mi2') 

svsp.loops <- read_tsv('analysis/10_loops_chromatin/contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>% 
  select(-7,-8) %>% 
  filter(chrA=='scaffold-mi2') 

svsp.TADs <- read_tsv('analysis/10_loops_chromatin/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= SVSP.reg.start & start <= SVSP.reg.end) %>% 
  rowid_to_column() %>% 
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


p.svsp.genes <-  ggplot(SVSP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVSP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_segment(aes(x=SVSP.reg.start*0.9,xend=SVSP.reg.end*1.1,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=svsp.vpers,aes(x=prom_start,xend=start,y='mi2',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=svsp.vpers,aes(x=SVSP.reg.start*0.9,xend=SVSP.reg.end*1.1,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svsp.vpers, aes(x=(start+end)/2, y=type),size=2) +
  ylab('') +
  xlab('') +
  ggtitle('Tandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(SVSP.reg.start,SVSP.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())


# TADS plot
p.svsp.tads <- ggplot(svsp.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(SVSP.reg.start,SVSP.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.svsp.loops <- ggplot(svsp.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  geom_curve(data=svsp.CTCFloops,color='red',lwd=1) +
  coord_cartesian(xlim=c(SVSP.reg.start,SVSP.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.svsp.h3k4me3.chip <- ggplot(svsp.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.h3k4me3.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  coord_cartesian(expand = T,xlim=c(SVSP.reg.start,SVSP.reg.end),ylim=c(0,svsp.h3k4me3.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K4me3)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svsp.h3k27ac.chip <- ggplot(svsp.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=svsp.h3k27ac.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  coord_cartesian(expand = T,xlim=c(SVSP.reg.start,SVSP.reg.end),ylim=c(0,svsp.h3k27ac.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K27ac)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svsp.ctcf.chip <- ggplot(svsp.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svsp.ctcf.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=svsp.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = svsp.ctcf.chip.peaks,aes(x=start-200,xend=end+200,y=svsp.ctcf.regionMax$max*1.2,yend=svsp.ctcf.regionMax$max*1.2),lwd=5,alpha=1) +
  geom_point(inherit.aes = F,data=svsp.ctcf.bound,aes(x=(start+end)/2),y=svsp.ctcf.regionMax$max*1.1,pch=25,size=4,fill='purple4',alpha=0.75) +
  # geom_point(inherit.aes = F,data=svsp.ctcf.bound,aes(x=(start+end)/2),y=svsp.ctcf.regionMax$max*1.1,pch=18,size=2.5) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(SVSP.reg.start,SVSP.reg.end),ylim=c(0,svsp.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svsp.atac.mean <- ggplot(svsp.atac.mean,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.atac.mean.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=svsp.atac.mean.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  coord_cartesian(expand = T,xlim=c(SVSP.reg.start,SVSP.reg.end),ylim=c(0,svsp.atac.mean.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('ATAC-seq (Post-Extraction Venom Gland)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svsp.scale <- ggplot(SVSP.info,aes(x=SVSP.reg.start,xend=SVSP.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=SVSP.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=SVSP.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(SVSP.reg.start,SVSP.reg.end)),label=SVSP.reg.length)) +
  coord_cartesian(expand = T,xlim=c(SVSP.reg.start,SVSP.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

# svsp.fullPlots <- p.svsp.tads / p.svsp.ctcf.chip / p.svsp.loops / p.svsp.genes  /  p.svsp.h3k27ac.chip / p.svsp.h3k4me3.chip/ p.svsp.atac.mean / p.svsp.atac.vgu + p.svsp.scale + plot_layout(ncol = 1,heights = c(1,2,2,2,2,2,2,2,1))
# svsp.fullPlots


svsp.all.peaks <- svsp.h3k27ac.chip.peaks %>% 
  bind_rows(svsp.h3k4me3.chip.peaks,svsp.atac.mean.peaks,svsp.ses) %>% 
  mutate(type = factor(type,levels = c('H3K4me3','H3K27ac','Super-Enhancers','ATAC-seq (Mean)')))

p.svsp.all.peaks <- ggplot(svsp.all.peaks,aes(x=start,xend=end, y=type,yend=type,color=type)) +
  geom_point(data=svsp.all.peaks %>% filter(type!='Super-Enhancers'),show.legend = F,size=1,alpha=0.8) +
  geom_segment(inherit.aes = T, data=svsp.all.peaks %>% filter(type=='Super-Enhancers'),lwd=2,show.legend = F) +
  scale_color_manual(values = c('H3K4me3'='#2C8790','H3K27ac'='#82B47A','ATAC-seq (Mean)'='#824029','Super-Enhancers'='goldenrod2')) +
  coord_cartesian(expand = T,xlim=c(SVSP.reg.start,SVSP.reg.end)) +
  scale_x_continuous(labels=comma) +
  xlab('Position') +
  theme_linedraw(base_size = 14)+ theme(axis.title.y = element_blank())
p.svsp.all.peaks

p.svsp.condensed <- p.svsp.tads / p.svsp.ctcf.chip / p.svsp.loops / p.svsp.genes  / p.svsp.all.peaks + plot_layout(heights = c(1,2,2,2,2))

p.svsp.condensed


# PLA2 --------------------------------------------------------------------


PLA2.info <- all_info %>% 
  filter(str_detect(gene,'PLA2')) 

PLA2.reg.start <- min(PLA2.info$start)-80000
PLA2.reg.end <- max(PLA2.info$end)+80000

PLA2.reg.length <- paste(c(round((PLA2.reg.end-PLA2.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k4me3 chipseq
pla2.h3k4me3.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi7.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type = 'H3K4me3')

pla2.h3k4me3.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type = 'H3K4me3')

pla2.h3k4me3.regionMax <- pla2.h3k4me3.chip %>% 
  filter(end >= PLA2.reg.start & start <= PLA2.reg.end) %>% 
  summarise(max=max(density))

# Read in  H3k27ac chipseq 
pla2.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi7.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type = 'H3K27ac')

pla2.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type = 'H3K27ac')

pla2.h3k27ac.regionMax <- pla2.h3k27ac.chip %>% 
  filter(end >= PLA2.reg.start & start <= PLA2.reg.end) %>% 
  summarise(max=max(density))

# Read in  CTCF chipseq and related files
pla2.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/CTCF_CroVir_i87_uniq_signal_mi7.bed',
                           col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) 

pla2.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                 col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) 

pla2.ctcf.bound <- read_tsv('analysis/10_loops_chromatin/contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) 

pla2.ctcf.regionMax <- pla2.ctcf.chip %>% 
  filter(end >= PLA2.reg.start & start <= PLA2.reg.end) %>% 
  summarise(max=max(density))


# Read in vPERs and super-enhancers
pla2.vpers <- read_tsv('analysis/6_ABC_Enhancers/ABC_output/_reformat/PLA2_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(PLA2.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

pla2.ses <- read_tsv('analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type='Super-Enhancers')


# Read in ATACseq
pla2.atac.mean <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ATACseq/VenomATACseq_ActiveMotif/BW_ATAC/bw2bedgraph/mi7_VG3_ATAC.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1)  %>% 
  mutate(type = 'ATAC-seq (Mean)')

pla2.atac.mean.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type = 'ATAC-seq (Mean)')

pla2.atac.mean.regionMax <- pla2.atac.mean %>% 
  filter(end >= PLA2.reg.start & start <= PLA2.reg.end) %>% 
  summarise(max=max(density))


# Read in loops and TADs

pla2.CTCFloops <- read_tsv('analysis/10_loops_chromatin/contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>% 
  select(-7,-8) %>% 
  filter(chrA=='scaffold-mi7') 

pla2.loops <- read_tsv('analysis/10_loops_chromatin/contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>% 
  select(-7,-8) %>% 
  filter(chrA=='scaffold-mi7') 

pla2.TADs <- read_tsv('analysis/10_loops_chromatin/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start & start <= PLA2.reg.end) %>% 
  rowid_to_column() %>% 
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


p.pla2.genes <-  ggplot(PLA2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_segment(aes(x=PLA2.reg.start*0.9,xend=PLA2.reg.end*1.1,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=pla2.vpers,aes(x=prom_start,xend=start,y='mi7',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=pla2.vpers,aes(x=PLA2.reg.start*0.9,xend=PLA2.reg.end*1.1,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=pla2.vpers, aes(x=(start+end)/2, y=type),size=2) +
  ylab('') +
  xlab('') +
  ggtitle('Tandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(PLA2.reg.start,PLA2.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())


# TADS plot
p.pla2.tads <- ggplot(pla2.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(PLA2.reg.start,PLA2.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.pla2.loops <- ggplot(pla2.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  geom_curve(data=pla2.CTCFloops,color='red',lwd=1) +
  coord_cartesian(xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.pla2.h3k4me3.chip <- ggplot(pla2.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = pla2.h3k4me3.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.h3k4me3.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K4me3)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.pla2.h3k27ac.chip <- ggplot(pla2.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = pla2.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=pla2.h3k27ac.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.h3k27ac.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K27ac)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.pla2.ctcf.chip <- ggplot(pla2.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_rect(inherit.aes = F,data = pla2.ctcf.chip.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=pla2.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  # geom_point(inherit.aes = F,data=pla2.ctcf.bound,aes(x=(start+end)/2),y=pla2.ctcf.regionMax$max*1.1,pch=18,size=2.5) +
  geom_segment(data = pla2.ctcf.chip.peaks,aes(x=start,xend=end,y=pla2.ctcf.regionMax$max*1.2,yend=pla2.ctcf.regionMax$max*1.2),lwd=5,alpha=1) +
  geom_point(inherit.aes = F,data=pla2.ctcf.bound,aes(x=(start+end)/2),y=pla2.ctcf.regionMax$max*1.1,pch=25,size=4,fill='purple4',alpha=0.75) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.pla2.atac.mean <- ggplot(pla2.atac.mean,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = pla2.atac.mean.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=pla2.atac.mean.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.atac.mean.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('ATAC-seq (Post-Extraction Venom Gland)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.pla2.scale <- ggplot(PLA2.info,aes(x=PLA2.reg.start,xend=PLA2.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=PLA2.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=PLA2.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(PLA2.reg.start,PLA2.reg.end)),label=PLA2.reg.length)) +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

# pla2.fullPlots <- p.pla2.tads / p.pla2.ctcf.chip / p.pla2.loops / p.pla2.genes  /  p.pla2.h3k27ac.chip / p.pla2.h3k4me3.chip/ p.pla2.atac.mean / p.pla2.atac.vgu + p.pla2.scale + plot_layout(ncol = 1,heights = c(1,2,2,2,2,2,2,2,1))
# pla2.fullPlots

pla2.all.peaks <- pla2.h3k27ac.chip.peaks %>% 
  bind_rows(pla2.h3k4me3.chip.peaks,pla2.atac.mean.peaks,pla2.ses) %>% 
  mutate(type = factor(type,levels = c('H3K4me3','H3K27ac','Super-Enhancers','ATAC-seq (Mean)')))

p.pla2.all.peaks <- ggplot(pla2.all.peaks,aes(x=start,xend=end, y=type,yend=type,color=type)) +
  geom_point(data=pla2.all.peaks %>% filter(type!='Super-Enhancers'),show.legend = F,size=1,alpha=0.8) +
  geom_segment(inherit.aes = T, data=pla2.all.peaks %>% filter(type=='Super-Enhancers'),lwd=2,show.legend = F) +
  scale_color_manual(values = c('H3K4me3'='#2C8790','H3K27ac'='#82B47A','ATAC-seq (Mean)'='#824029','Super-Enhancers'='goldenrod2')) +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end)) +
  scale_x_continuous(labels=comma) +
  xlab('Position') +
  theme_linedraw(base_size = 14)+ theme(axis.title.y = element_blank())
p.pla2.all.peaks


p.pla2.condensed <- p.pla2.tads / p.pla2.ctcf.chip / p.pla2.loops / p.pla2.genes  / p.pla2.all.peaks + plot_layout(heights = c(1,2,2,2,2))


### Plot all!

p.svmp.condensed | p.svsp.condensed | p.pla2.condensed





# Hi-C Plotting -----------------------------------------------------------

#BiocManager::install("Sushi")
library(Sushi)

# SVMP
#Load and format HiC matrix
svmp.hic.long <- read.table('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/z_old_and_misc/multiTrack_plotting/data/misc_data/mi1_intra_10kb_contact_matrix_01.16.19.txt',sep='\t',stringsAsFactors = F)
svmp.hic.long$V3 <- log10(svmp.hic.long$V3+.1)
colnames(svmp.hic.long) <- c('b1','b2','count')

svmp.hic.wide <- spread(svmp.hic.long,key = b2,value = count, fill=0)
row.names(svmp.hic.wide) <- svmp.hic.wide$b1
svmp.hic.wide <- svmp.hic.wide[,-1]
svmp.hic.wide[lower.tri(svmp.hic.wide)] = t(svmp.hic.wide)[lower.tri(svmp.hic.wide)]

# SVSP
#Load and format HiC matrix
svsp.hic.long <- read.table('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/z_old_and_misc/multiTrack_plotting/data/misc_data/mi2_intra_10kb_contact_matrix_01.16.19.txt',sep='\t',stringsAsFactors = F)
svsp.hic.long$V3 <- log10(svsp.hic.long$V3+.1)
colnames(svsp.hic.long) <- c('b1','b2','count')

svsp.hic.wide <- spread(svsp.hic.long,key = b2,value = count, fill=0)
row.names(svsp.hic.wide) <- svsp.hic.wide$b1
svsp.hic.wide <- svsp.hic.wide[,-1]
svsp.hic.wide[lower.tri(svsp.hic.wide)] = t(svsp.hic.wide)[lower.tri(svsp.hic.wide)]

# PLA2
#Load and format HiC matrix
pla2.hic.long <- read.table('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/z_old_and_misc/multiTrack_plotting/data/misc_data/mi7_intra_10kb_contact_matrix_01.16.19.txt',sep='\t',stringsAsFactors = F)
pla2.hic.long$V3 <- log10(pla2.hic.long$V3+.1)
colnames(pla2.hic.long) <- c('b1','b2','count')

pla2.hic.wide <- spread(pla2.hic.long,key = b2,value = count, fill=0)
row.names(pla2.hic.wide) <- pla2.hic.wide$b1
pla2.hic.wide <- pla2.hic.wide[,-1]
pla2.hic.wide[lower.tri(pla2.hic.wide)] = t(pla2.hic.wide)[lower.tri(pla2.hic.wide)]



# Plot

chrom_svmp = 'scaffold-mi1'
chrom_svsp = 'scaffold-mi2'
chrom_pla2 = 'scaffold-mi7'

svmp.hic.start = (SVMP.reg.start+SVMP.reg.end)/2 - 1000000
svmp.hic.end = (SVMP.reg.start+SVMP.reg.end)/2 + 1000000

svsp.hic.start = (SVSP.reg.start+SVSP.reg.end)/2 - 1000000
svsp.hic.end = (SVSP.reg.start+SVSP.reg.end)/2 + 1000000

pla2.hic.start = (PLA2.reg.start+PLA2.reg.end)/2 - 1000000
pla2.hic.end = (PLA2.reg.start+PLA2.reg.end)/2 + 1000000


#SVMP
plotHic(svmp.hic.wide,chrom_svmp,svmp.hic.start,svmp.hic.end,max_y=36,zrange = c(0,1),palette = viridis)
labelgenome(chrom_svmp,svmp.hic.start,svmp.hic.end,n=10,scale="Mb",chromcex = .75,scalecex = .75,edgeblankfraction = 0.1)
zoomsregion(c(p.svmp.tads$coordinates$limits$x[1],p.svmp.tads$coordinates$limits$x[2]),wideextend=0.14,extend=c(0.01,0.13),
            offsets=c(0,0),lty=3,lwd=2)


#SVSP
plotHic(svsp.hic.wide,chrom_svsp,svsp.hic.start,svsp.hic.end,max_y=36,zrange = c(0,1),palette = viridis)
labelgenome(chrom_svsp,svsp.hic.start,svsp.hic.end,n=10,scale="Mb",chromcex = .75,scalecex = .75,edgeblankfraction = 0.1)
zoomsregion(c(p.svsp.tads$coordinates$limits$x[1],p.svsp.tads$coordinates$limits$x[2]),wideextend=0.14,extend=c(0.01,0.13),
            offsets=c(0,0),lty=3,lwd=2)

#PLA2
plotHic(pla2.hic.wide,chrom_pla2,pla2.hic.start,pla2.hic.end,max_y=36,zrange = c(0,1),palette = viridis)
labelgenome(chrom_pla2,pla2.hic.start,pla2.hic.end,n=10,scale="Mb",chromcex = .75,scalecex = .75,edgeblankfraction = 0.1)
zoomsregion(c(p.pla2.tads$coordinates$limits$x[1],p.pla2.tads$coordinates$limits$x[2]),wideextend=0.14,extend=c(0.01,0.13),
            offsets=c(0,0),lty=3,lwd=2)
