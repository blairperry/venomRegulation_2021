
#install.packages('gggenes')
library(gggenes)
library(viridis)
library(patchwork)
library(tidyverse)
library(scales)
library(ggforce)



pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes,by=c('txid'='X6')) %>% 
  mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))


all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
  filter(str_detect(X9,'trnascan',negate = T)) %>% 
  mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>% 
  mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>% 
  filter(tx_id %in% pri_venom_genes$X6) %>% 
  left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>% 
  select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>% 
  mutate(strand = ifelse(strand == '+','forward','reverse')) %>% 
  mutate(direction = ifelse(strand == 'forward',1,-1)) %>% 
  mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>% 
  left_join(exp) %>% 
  mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>% 
  mutate(prom_start = ifelse(strand=='forward',start,end))



# SVMP  -------------------------------------------------------------------

SVMP.info <- all_info %>% 
  filter(str_detect(gene,'SVMP|ADAM')) 

SVMP.reg.start <- min(SVMP.info$start)-20000
SVMP.reg.end <- max(SVMP.info$end)+20000

SVMP.reg.length <- paste(c(round((SVMP.reg.end-SVMP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k4me3 chipseq
svmp.h3k4me3.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi1.bed',
                      col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

svmp.h3k4me3.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 


# Read in  H3k27ac chipseq 
svmp.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi1.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

svmp.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

# Read in vPERs and super-enhancers

svmp.vpers <- read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(SVMP.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

svmp.ses <- read_tsv('./analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 


# Read in ATACseq
svmp.atac.rep1 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG1_mi1_scaleFactorNorm_atacseqDensity.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end)  %>% 
  mutate(replicate = 'rep1')

svmp.atac.rep1.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG1_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) %>% 
  mutate(replicate = 'rep1')

svmp.atac.rep2 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG11_mi1_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end)  %>% 
  mutate(replicate = 'rep2')

svmp.atac.rep2.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG11_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end)  %>% 
  mutate(replicate = 'rep2')

svmp.atac.rep3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG4_mi1_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end)  %>% 
  mutate(replicate = 'rep3')

svmp.atac.rep3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG4_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end)  %>% 
  mutate(replicate = 'rep3')


svmp.atac.allReps <- svmp.atac.rep1 %>% 
  bind_rows(svmp.atac.rep2,svmp.atac.rep3)

svmp.atac.allReps.peaks <- svmp.atac.rep1.peaks %>% 
  bind_rows(svmp.atac.rep2.peaks,svmp.atac.rep3.peaks)



svmp.atac.mean <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi1.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end)  

svmp.atac.mean.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3) %>% 
  filter(chr=='scaffold-mi1') %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 



# Plotting
# NOTE: widening peak regions by 200bp both directions to make more visible in plots

p.svmp.genes <-  ggplot(SVMP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(Median1DPE+1))) +
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
  
  ylab('') +
  xlab('') +
  ggtitle('a) SVMP\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())


p.svmp.h3k4me3.chip <- ggplot(svmp.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.h3k4me3.chip.peaks,aes(x=start-200,xend=end+200,y=max(svmp.h3k4me3.chip$density)*1.2,yend=max(svmp.h3k4me3.chip$density)*1.2),lwd=5,alpha=1) +
  geom_segment(color='#2C8790') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.h3k4me3.chip$density)*1.2)) +
  ggtitle('Chip-seq (H3K4me3)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.svmp.h3k27ac.chip <- ggplot(svmp.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.h3k27ac.chip.peaks,aes(x=start-200,xend=end+200,y=max(svmp.h3k27ac.chip$density)*1.2,yend=max(svmp.h3k27ac.chip$density)*1.2),lwd=5,alpha=1) + 
  geom_segment(color='#82B47A') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.h3k27ac.chip$density)*1.2)) +
  ggtitle('Chip-seq (H3K27ac)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svmp.atac.mean <- ggplot(svmp.atac.mean,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.atac.mean.peaks,aes(x=start-200,xend=end+200,y=max(svmp.atac.mean$density)*1.2,yend=max(svmp.atac.mean$density)*1.2),lwd=5,alpha=1) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.atac.mean$density)*1.2)) +
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
# 
# svmp.fullPlots <- p.svmp.genes  / p.svmp.h3k27ac.chip / p.svmp.h3k4me3.chip/ p.svmp.atac.rep1 + p.svmp.scale + plot_layout(ncol = 1,heights = c(2,2,2,2,1))
# svmp.fullPlots


# Testing plotting multiple ATAC reps

p.svmp.atac.allReps <- ggplot(svmp.atac.allReps,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svmp.atac.allReps.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svmp.atac.allReps$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#C56A4B') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  # scale_y_continuous(limits=c(0,max(svmp.atac.allReps$density)*1.2)) +
  facet_grid(rows=vars(replicate),scales = 'free_y') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),
                                        strip.text.y = element_blank(),strip.background = element_blank())


svmp.fullPlots.v2 <- p.svmp.genes  / p.svmp.h3k27ac.chip / p.svmp.h3k4me3.chip/ p.svmp.atac.mean / p.svmp.atac.allReps + p.svmp.scale + plot_layout(ncol = 1,heights = c(2,2,2,2,3,1))
svmp.fullPlots.v2

# SVSP  -------------------------------------------------------------------

svsp.info <- all_info %>% 
  filter(str_detect(gene,'SVSP')) 

svsp.reg.start <- min(svsp.info$start)-100000
svsp.reg.end <- max(svsp.info$end)+100000

svsp.reg.length <- paste(c(round((svsp.reg.end-svsp.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  H3k4me3 chipseq
svsp.h3k4me3.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi2.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

svsp.h3k4me3.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 


# Read in  H3k27ac chipseq 
svsp.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi2.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

svsp.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

# Read in vPERs and super-enhancers

svsp.vpers <- read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(svsp.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

svsp.ses <- read_tsv('./analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 


# Read in ATACseq
svsp.atac.rep1 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG1_mi2_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end)  %>% 
  mutate(replicate = 'rep1')

svsp.atac.rep1.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG1_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) %>% 
  mutate(replicate = 'rep1')

svsp.atac.rep2 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG11_mi2_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end)  %>% 
  mutate(replicate = 'rep2')

svsp.atac.rep2.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG11_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end)  %>% 
  mutate(replicate = 'rep2')

svsp.atac.rep3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG4_mi2_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end)  %>% 
  mutate(replicate = 'rep3')

svsp.atac.rep3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG4_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end)  %>% 
  mutate(replicate = 'rep3')


svsp.atac.allReps <- svsp.atac.rep1 %>% 
  bind_rows(svsp.atac.rep2,svsp.atac.rep3)

svsp.atac.allReps.peaks <- svsp.atac.rep1.peaks %>% 
  bind_rows(svsp.atac.rep2.peaks,svsp.atac.rep3.peaks)



svsp.atac.mean <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi2.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end)  

svsp.atac.mean.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 



# Plotting


p.svsp.genes <-  ggplot(svsp.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVSP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_segment(aes(x=svsp.reg.start,xend=svsp.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=svsp.vpers,aes(x=prom_start,xend=start,y='mi2',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=svsp.vpers,aes(x=svsp.reg.start,xend=svsp.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svsp.vpers, aes(x=(start+end)/2, y=type),size=2) +
  
  ylab('') +
  xlab('') +
  ggtitle('b) SVSP\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())
# p.svsp.genes



p.svsp.h3k4me3.chip <- ggplot(svsp.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svsp.h3k4me3.chip.peaks,aes(x=(start+end)/2,
  #                                                 xend=(start+end)/2,
  #                                                 y=0,
  #                                                 yend=max(svsp.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = svsp.h3k4me3.chip.peaks,aes(x=start-200,xend=end+200,y=max(svsp.h3k4me3.chip$density)*1.2,yend=max(svsp.h3k4me3.chip$density)*1.2),lwd=5,alpha=1) +
  geom_segment(color='#2C8790') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.h3k4me3.chip$density)*1.2)) +
  ggtitle('Chip-seq (H3K4me3)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.svsp.h3k4me3.chip


p.svsp.h3k27ac.chip <- ggplot(svsp.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svsp.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.h3k27ac.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = svsp.h3k27ac.chip.peaks,aes(x=start-200,xend=end+200,y=max(svsp.h3k27ac.chip$density)*1.2,yend=max(svsp.h3k27ac.chip$density)*1.2),lwd=5,alpha=1) +
  geom_segment(color='#82B47A') +
  geom_segment(inherit.aes = F,
               data=svsp.ses,
               aes(x=start,xend=end,y=max(svsp.h3k27ac.chip$density)*1.45,yend=max(svsp.h3k27ac.chip$density)*1.45),
               lwd=3,
               color='goldenrod2') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.h3k27ac.chip$density)*1.5)) +
  ggtitle('Chip-seq (H3K27ac) - Super-Enhancers Shown in Yellow') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.svsp.h3k27ac.chip

p.svsp.atac.mean <- ggplot(svsp.atac.mean,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svsp.atac.mean.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.atac.mean$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = svsp.atac.mean.peaks,aes(x=start-200,xend=end+200,y=max(svsp.atac.mean$density)*1.2,yend=max(svsp.atac.mean$density)*1.2),lwd=5,alpha=1) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.atac.mean$density)*1.2)) +
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

# svsp.fullPlots <- p.svsp.genes / p.svsp.h3k27ac.chip / p.svsp.h3k4me3.chip/ p.svsp.atac.rep1 + p.svsp.scale + plot_layout(ncol = 1,heights = c(2,2,2,2,1))
# svsp.fullPlots


# Testing plotting multiple ATAC reps

p.svsp.atac.allReps <- ggplot(svsp.atac.allReps,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = svsp.atac.allReps.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.atac.allReps$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#C56A4B') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  # scale_y_continuous(limits=c(0,max(svsp.atac.allReps$density)*1.2)) +
  facet_grid(rows=vars(replicate),scales = 'free_y') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),
                                        strip.text.y = element_blank(),strip.background = element_blank())


svsp.fullPlots.v2 <- p.svsp.genes  / p.svsp.h3k27ac.chip / p.svsp.h3k4me3.chip/ p.svsp.atac.mean / p.svsp.atac.allReps + p.svsp.scale + plot_layout(ncol = 1,heights = c(2,2,2,2,3,1))
svsp.fullPlots.v2


# PLA2  ------------------------------------------------------------------- <- 

pla2.info <- all_info %>% 
  filter(str_detect(gene,'PLA2')) 

pla2.reg.start <- 3001808-5000
pla2.reg.end <- 3050307+5000
pla2.reg.length <- paste(c(round((pla2.reg.end-pla2.reg.start)/1000,digits = 2),'kb'),collapse = ' ')


# Read in  H3k4me3 chipseq
pla2.h3k4me3.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi7.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

pla2.h3k4me3.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 


# Read in  H3k27ac chipseq 
pla2.h3k27ac.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K27Ac_CroVir_i88_uniq_signal_mi7.bed',
                              col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

pla2.h3k27ac.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K27Ac_VenomGland_intervals_simple.bed',
                                    col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

# Read in vPERs and super-enhancers

pla2.vpers <- read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/PLA2_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(pla2.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

pla2.ses <- read_tsv('./analysis/2_super_enhancers/bed/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 


# Read in ATACseq
pla2.atac.rep1 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG1_mi7_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end)  %>% 
  mutate(replicate = 'rep1')

pla2.atac.rep1.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG1_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) %>% 
  mutate(replicate = 'rep1')

pla2.atac.rep2 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG11_mi7_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end)  %>% 
  mutate(replicate = 'rep2')

pla2.atac.rep2.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG11_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end)  %>% 
  mutate(replicate = 'rep2')

pla2.atac.rep3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/RVG4_mi7_scaleFactorNorm_atacseqDensity.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end)  %>% 
  mutate(replicate = 'rep3')

pla2.atac.rep3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/sorted/RVG4_PostExt_ATAC.dedup.filtered.unique.norm_peaks.sorted.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3,sample=4) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end)  %>% 
  mutate(replicate = 'rep3')


pla2.atac.allReps <- pla2.atac.rep1 %>% 
  bind_rows(pla2.atac.rep2,pla2.atac.rep3)

pla2.atac.allReps.peaks <- pla2.atac.rep1.peaks %>% 
  bind_rows(pla2.atac.rep2.peaks,pla2.atac.rep3.peaks)



pla2.atac.mean <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi7.bedgraph',
                           col_names = c('chr','start','end','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end)  

pla2.atac.mean.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',col_names = F) %>% 
  select(chr=1,start=2,end=3) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 





# Plotting


p.pla2.genes <-  ggplot(pla2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_segment(aes(x=pla2.reg.start,xend=pla2.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=pla2.vpers,aes(x=prom_start,xend=start,y='mi7',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=pla2.vpers,aes(x=pla2.reg.start,xend=pla2.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=pla2.vpers, aes(x=(start+end)/2, y=type),size=2) +
  
  ylab('') +
  xlab('') +
  ggtitle('c) PLA2\n\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())
p.pla2.genes

p.pla2.h3k4me3.chip <- ggplot(pla2.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_rect(inherit.aes = F, data = pla2.h3k4me3.chip.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=max(pla2.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = pla2.h3k4me3.chip.peaks,aes(x=start,xend=end,y=max(pla2.h3k4me3.chip$density)*1.2,yend=max(pla2.h3k4me3.chip$density)*1.2),lwd=5,alpha=1) +
  # geom_segment(data = pla2.h3k4me3.chip.peaks,aes(x=(start+end)/2,
  #                                                 xend=(start+end)/2,
  #                                                 y=0,
  #                                                 yend=max(pla2.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.h3k4me3.chip$density)*1.2)) +
  ggtitle('Chip-seq (H3K4me3)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.pla2.h3k4me3.chip


p.pla2.h3k27ac.chip <- ggplot(pla2.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_rect(inherit.aes = F, data = pla2.h3k27ac.chip.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=max(pla2.h3k27ac.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = pla2.h3k27ac.chip.peaks,aes(x=start,xend=end,y=max(pla2.h3k27ac.chip$density)*1.2,yend=max(pla2.h3k27ac.chip$density)*1.2),lwd=5,alpha=1) +
  # geom_segment(data = pla2.h3k27ac.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.h3k27ac.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#82B47A') +
  geom_segment(inherit.aes = F,
               data=pla2.ses,
               aes(x=start,xend=end,y=max(pla2.h3k27ac.chip$density)*1.45,yend=max(pla2.h3k27ac.chip$density)*1.45),
               lwd=3,
               color='goldenrod2') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.h3k27ac.chip$density)*1.5)) +
  ggtitle('Chip-seq (H3K27ac) - Super-Enhancers Shown in Yellow') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.pla2.h3k27ac.chip


p.pla2.atac.mean <- ggplot(pla2.atac.mean,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_rect(inherit.aes = F, data = pla2.atac.mean.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=max(pla2.atac.mean$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(data = pla2.atac.mean.peaks,aes(x=start,xend=end,y=max(pla2.atac.mean$density)*1.2,yend=max(pla2.atac.mean$density)*1.2),lwd=5,alpha=1) +
  # geom_segment(data = pla2.atac.mean.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.atac.mean$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.atac.mean$density)*1.2)) +
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

# pla2.fullPlots <- p.pla2.genes  / p.pla2.h3k27ac.chip / p.pla2.h3k4me3.chip/ p.pla2.atac.rep1 + p.pla2.scale + plot_layout(ncol = 1,heights = c(2,2,2,2,1))
# pla2.fullPlots


p.pla2.atac.allReps <- ggplot(pla2.atac.allReps,aes(x=start,xend=start,y=0,yend=density)) +
  # geom_segment(data = pla2.atac.allReps.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.atac.allReps$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#C56A4B') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  # scale_y_continuous(limits=c(0,max(pla2.atac.allReps$density)*1.2)) +
  facet_grid(rows=vars(replicate),scales = 'free_y') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),
                                        strip.text.y = element_blank(),strip.background = element_blank())


pla2.fullPlots.v2 <- p.pla2.genes  / p.pla2.h3k27ac.chip / p.pla2.h3k4me3.chip/ p.pla2.atac.mean / p.pla2.atac.allReps + p.pla2.scale + plot_layout(ncol = 1,heights = c(2,2,2,2,3,1))
pla2.fullPlots.v2


### Plot all 3

svmp.fullPlots.v2 | svsp.fullPlots.v2 | pla2.fullPlots.v2
