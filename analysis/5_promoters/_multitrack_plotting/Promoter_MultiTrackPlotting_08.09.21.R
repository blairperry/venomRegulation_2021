
#install.packages('gggenes')
library(gggenes)
library(viridis)
library(patchwork)
library(tidyverse)
library(scales)


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

# Read in ChIPseq
svmp.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi1.bed',
                      col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= SVMP.reg.start & start <= SVMP.reg.end) 

svmp.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                            col_names = c('chr','start','end','id')) %>% 
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

p.svmp.genes <- ggplot(SVMP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVMP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'ADAM')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = -0.5,size=3,color='grey40') +
  geom_segment(aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('SVMP\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())
p.svmp.genes

p.svmp.chip <- ggplot(svmp.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svmp.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.chip$density)*1.2)) +
  ggtitle('ChIP-seq - H3K4me3') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svmp.atac.vg3 <- ggplot(svmp.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svmp.atac.vg3.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svmp.atac.vg3$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svmp.atac.vg3$density)*1.2)) +
  ggtitle('ATAC-seq - Post-Extraction VG') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
p.svmp.atac.vg3

p.svmp.scale <- ggplot(SVMP.info,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_label(aes(x=mean(c(SVMP.reg.start,SVMP.reg.end)),label=SVMP.reg.length)) +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_void()

svmp.fullPlots <- p.svmp.genes / p.svmp.chip / p.svmp.atac.vg3  + p.svmp.scale + plot_layout(ncol = 1,heights = c(1,1,1,0.5))

svmp.fullPlots


# SVSP  -------------------------------------------------------------------

svsp.info <- all_info %>% 
  filter(str_detect(gene,'SVSP')) 

svsp.reg.start <- min(svsp.info$start)-20000
svsp.reg.end <- max(svsp.info$end)+20000
svsp.reg.length <- paste(c(round((svsp.reg.end-svsp.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in ChIPseq
svsp.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi2.bed',
                      col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

svsp.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

# Read in ATACseq
svsp.atac.vg3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi2.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 

svsp.atac.vg3.peaks <- read_tsv('./analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi2') %>% 
  filter(end >= svsp.reg.start & start <= svsp.reg.end) 


# Plotting

p.svsp.genes <- ggplot(svsp.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVSP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  geom_segment(aes(x=svsp.reg.start,xend=svsp.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('SVSP\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())

p.svsp.chip <- ggplot(svsp.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  # geom_point(inherit.aes = F,
  #            data=svsp.chip.peaks,
  #            aes(x=(start+end)/2,y=max(svsp.chip$density)*1.1),
  #            color='#2C8790',
  #            pch=25,
  #            size=2) +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.chip$density)*1.2)) +
  ggtitle('ChIP-seq - H3K4me3') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svsp.atac.vg3 <- ggplot(svsp.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = svsp.atac.vg3.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(svsp.atac.vg3$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  # geom_point(inherit.aes = F,
  #            data=svsp.atac.vg3.peaks,
  #            aes(x=(start+end)/2,y=max(svsp.atac.vg3$density)*1.1),
  #            color='#824029',
  #            pch=25,
  #            size=2) +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(svsp.atac.vg3$density)*1.2)) +
  ggtitle('ATAC-seq - Post-Extraction VG') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.svsp.scale <- ggplot(svsp.info,aes(x=svsp.reg.start,xend=svsp.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_label(aes(x=mean(c(svsp.reg.end,svsp.reg.start)),label=svsp.reg.length)) +
  scale_x_continuous(labels = comma,limits=c(svsp.reg.start,svsp.reg.end),expand=c(0,0)) +
  theme_void()

svsp.fullPlots <- p.svsp.genes / p.svsp.chip / p.svsp.atac.vg3 / p.svsp.scale + plot_layout(ncol = 1,heights =c(1,1,1,0.5))

svsp.fullPlots




# PLA2  -------------------------------------------------------------------

pla2.info <- all_info %>% 
  filter(str_detect(gene,'PLA2')) 

pla2.reg.start <- min(pla2.info$start)-5000
pla2.reg.end <- max(pla2.info$end)+5000
pla2.reg.length <- paste(c(round((pla2.reg.end-pla2.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in ChIPseq
pla2.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/venomScaffoldSubsets/H3K4me3_CroVir_i89_uniq_signal_mi7.bed',
                      col_names = c('chr','start','end','id','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

pla2.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/H3K4me3_VenomGland_intervals_simple.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

# Read in ATACseq
pla2.atac.vg3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi7.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 

pla2.atac.vg3.peaks <- read_tsv('./analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.simple.bed',
                                col_names = c('chr','start','end')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= pla2.reg.start & start <= pla2.reg.end) 


# Plotting

p.pla2.genes <- ggplot(pla2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  geom_segment(aes(x=pla2.reg.start,xend=pla2.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('PLA2\nTandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x = element_blank())

p.pla2.chip <- ggplot(pla2.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_rect(inherit.aes = F, data = pla2.chip.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=max(pla2.chip$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#2C8790') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.chip$density)*1.2)) +
  ggtitle('ChIP-seq - H3K4me3') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))


p.pla2.atac.vg3 <- ggplot(pla2.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_rect(inherit.aes = F, data = pla2.atac.vg3.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=max(pla2.atac.vg3$density)*1.2),lwd=0.5,alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,max(pla2.atac.vg3$density)*1.2)) +
  ggtitle('ATAC-seq - Post-Extraction VG') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.pla2.scale <- ggplot(pla2.info,aes(x=pla2.reg.start,xend=pla2.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_label(aes(x=mean(c(pla2.reg.end,pla2.reg.start)),label=pla2.reg.length)) +
  scale_x_continuous(labels = comma,limits=c(pla2.reg.start,pla2.reg.end),expand=c(0,0)) +
  theme_void()

pla2.fullPlots <- p.pla2.genes / p.pla2.chip / p.pla2.atac.vg3 / p.pla2.scale + plot_layout(ncol = 1,heights =c(1,1,1,0.5))



# Plotting all 3 families -------------------------------------------------

svmp.fullPlots | svsp.fullPlots | pla2.fullPlots
