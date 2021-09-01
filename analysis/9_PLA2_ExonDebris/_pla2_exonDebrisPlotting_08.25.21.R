library(scales)
library(tidyverse)
library(gggenes)
library(ggforce)
library(colorspace)
library(patchwork)


blast.res <- read_tsv('analysis/9_PLA2_ExonDebris/PLA2_coreVPERBlast_ExonDebrisOverlap_08.25.21.txt',col_names = c('chr.a','start.a','end.a','name.a','chr.b','start.b','end.b','name.b'))

ex.debris <- read_tsv('analysis/9_PLA2_ExonDebris/PLA2_ExonDebrisLocations.bed',col_names = c('chr','start','end','name')) %>% 
  mutate(desc= str_remove(name,'Name=ExonDebris_'))

pla2giie.exons <- read_tsv('analysis/9_PLA2_ExonDebris/Pla2_annotation_plus_exons.gff',skip = 1,col_names = c('chr','program','type','start','end','score','strand','score2','desc')) %>% 
  filter(str_detect(desc,'PLA2gIIE') & type == 'exon') %>% 
  mutate(desc= str_remove(desc,'Name='))

pri_venom_genes <- read_tsv('data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)

exp <- read_tsv('analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  # mutate(txid = ifelse(str_detect(txid,'ADAM28',negate = T),str_replace_all(txid,'_',' '),txid)) %>% 
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




PLA2.info <- all_info %>% 
  filter(str_detect(gene,'PLA2')) 

PLA2.reg.start <- min(PLA2.info$start)-3000
PLA2.reg.end <- max(PLA2.info$end)+3000

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

pla2.ses <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_H3K27Ac_CSeq_SE_analysis/SuperEnhancer_VenomGland_intervals_simple.bed',col_names = 
                       c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type='Super-Enhancers')


# Read in ATACseq
pla2.atac.vg3 <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/6_deepTools/scaleFactorNorm/bedgraph/focal_scaffolds/allSampleMean.scaleFactorNorm_mi7.bedgraph',
                          col_names = c('chr','start','end','density')) %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1)  %>% 
  mutate(type = 'ATAC-seq (VG3)')

pla2.atac.vg3.peaks <- read_tsv('analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.tsv',
                                col_types = list(chr=col_character(), start=col_double(), end=col_double(), samples=col_character())) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= PLA2.reg.start*0.9 & start <= PLA2.reg.end*1.1) %>% 
  mutate(type = 'ATAC-seq (VG3)')

pla2.atac.vg3.regionMax <- pla2.atac.vg3 %>% 
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


## Plotting

p.pla2.genes <- ggplot(PLA2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  
  geom_rect(inherit.aes = F,data=pla2giie.exons,aes(xmin=start,xmax=end,ymin='  vPER Blast Hits',ymax='mi7'),fill='darkblue',alpha=0.5) +
  
  geom_segment(aes(x=PLA2.reg.start*0.9,xend=PLA2.reg.end*1.1,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  
  geom_diagonal(inherit.aes = F,data=pla2.vpers,aes(x=prom_start,xend=(start+end)/2,y='mi7',yend=type,alpha = stat(index)),strength = 0.1,show.legend = F) +
  
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  
  geom_segment(inherit.aes = F,data=pla2.vpers,aes(x=PLA2.reg.start*0.9,xend=PLA2.reg.end*1.1,y=type,yend=type),lwd=1,color='grey70') +
  geom_segment(inherit.aes = F, data=pla2.vpers, aes(x=start,xend=end, y=type,yend=type),size=2) +
  
  geom_point(inherit.aes = F, data=ex.debris %>% filter(str_detect(desc,'g2E')), aes(x=(start+end)/2, y=' Exon Debris',color=desc),size=4) +
  scale_color_manual(values = c('g2E_e2'='#8EC164','g2E_e3'='#197E51')) +
  
  geom_point(inherit.aes = F, data=blast.res, aes(x=(start.a+end.a)/2,y='  vPER Blast Hits'),color='firebrick4',lwd=4) +
  
  ylab('') +
  xlab('') +
  ggtitle('Tandem Array with Gene Expression at 1DPE') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(PLA2.reg.start,PLA2.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())

p.pla2.genes

# TADS plot
# p.pla2.tads <- ggplot(pla2.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
#   geom_segment(size=4,show.legend = F) +
#   scale_color_continuous_sequential('OrRd',begin = 0.65) +
#   coord_cartesian(xlim=c(PLA2.reg.start,PLA2.reg.end),ylim = c(-10,10),expand = T) +
#   scale_x_continuous(labels = comma)+
#   ggtitle('TADs') +
#   theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
#                                         plot.title = element_text(color='black',face='bold',size = 14),
#                                         axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())
# 
# p.pla2.loops <- ggplot(pla2.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
#   # geom_segment(size=8,show.legend = F) +
#   geom_curve() +
#   geom_curve(data=pla2.CTCFloops,color='red',lwd=1) +
#   coord_cartesian(xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(-20,0),expand = T) +
#   scale_x_continuous(labels = comma)+
#   ggtitle('Loops') +
#   theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
#                                         plot.title = element_text(color='black',face='bold',size = 14),
#                                         axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())
# 
# 
# p.pla2.h3k4me3.chip <- ggplot(pla2.h3k4me3.chip,aes(x=start,xend=start,y=0,yend=density)) +
#   geom_segment(data = pla2.h3k4me3.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=max(pla2.h3k4me3.chip$density)*1.2),lwd=0.5,alpha=0.2) +
#   geom_segment(color='#2C8790') +
#   coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.h3k4me3.regionMax$max*1.2)) +
#   scale_x_continuous(labels = comma) +
#   ggtitle('Chip-seq (H3K4me3)') +
#   ylab('Read\nDensity') +
#   theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.pla2.h3k27ac.chip <- ggplot(pla2.h3k27ac.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_rect(data = pla2.h3k27ac.chip.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=pla2.h3k27ac.regionMax$max*1.2),inherit.aes = F,fill='#82B47A',alpha=0.2) +
  geom_segment(color='#82B47A') +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.h3k27ac.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (H3K27ac)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))
 


p.pla2.atac.vg3 <- ggplot(pla2.atac.vg3,aes(x=start,xend=start,y=0,yend=density)) +
  geom_rect(data = pla2.atac.vg3.peaks,aes(xmin=start,xmax=end,ymin=0,ymax=pla2.atac.vg3.regionMax$max*1.2),inherit.aes = F,fill='#824029',alpha=0.2) +
  geom_segment(color='#824029') +
  scale_x_continuous(labels = comma) +
  coord_cartesian(expand = T,xlim=c(PLA2.reg.start,PLA2.reg.end),ylim=c(0,pla2.atac.vg3.regionMax$max*1.2)) +
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


pla2.fullPlots <- p.pla2.genes  /  p.pla2.h3k27ac.chip /  p.pla2.atac.vg3 / p.pla2.scale + plot_layout(ncol = 1,heights = c(3,2,2,1))
pla2.fullPlots

