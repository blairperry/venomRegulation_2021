
#install.packages('gggenes')
library(gggenes)
library(viridis)
library(patchwork)
library(tidyverse)
library(scales)
library(ggforce)
library(colorspace)

###### NOTE: abandoned this before it was complete

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/transcription_factors/_tfbs_analysis/_Enhancer_CIIIDER_TFBS_EnrichmentAnalyses_12.11.20/_multiTrackPlots")


pri_venom_genes <- read_tsv('../../../../../data/venom_annotations/___PriorityVenomGenes_Sept2020/PriorityVenomGeneList_withNVPs_02.25.21.txt',col_names = F)
exp <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/gene_expression/VG_1DPEAvgExpression_wNVPs.txt') %>% 
  mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>% 
  mutate(log10avg1DPE = log10(avg1DPE + 0.1))


all_info <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
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
  mutate(prom_start = ifelse(strand=='forward',start,end)) %>% 
  filter(str_detect(gene,'PLA2|NVP|SVMP|SVSP',negate = T))



# CRISPs ------------------------------------------------------------------

crisp.info <- all_info %>% 
  filter(str_detect(gene,'CRISP')) 

crisp.reg.start <- min(crisp.info$start)-120000
crisp.reg.end <- max(crisp.info$end)+120000

crisp.reg.length <- paste(c(round((crisp.reg.end-crisp.reg.start)/1000,digits = 2),'kb'),collapse = ' ')


# Read in  CTCF chipseq and related files
crisp.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/1_047U_00JKUTA_VenomGland-495_CTCF_CroVir_i87_uniq_signal.bed',
                           col_names = c('chr','start','end','id','density')) %>% 
  filter(chr=='scaffold-ma1') %>% 
  filter(end >= crisp.reg.start & start <= crisp.reg.end) 

crisp.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                 col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-ma1') %>% 
  filter(end >= crisp.reg.start & start <= crisp.reg.end) 

crisp.ctcf.bound <- read_tsv('../../../../contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-ma1') %>% 
  filter(end >= crisp.reg.start & start <= crisp.reg.end) 

crisp.ctcf.regionMax <- crisp.ctcf.chip %>% 
  filter(end >= crisp.reg.start & start <= crisp.reg.end) %>% 
  summarise(max=max(density))

# Read in loops and TADs

crisp.CTCFloops <- read_tsv('../../../../contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-ma1')

crisp.loops <- read_tsv('../../../../contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-ma1')

crisp.TADs <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>%
  filter(chr=='scaffold-ma1') %>%
  filter(end >= crisp.reg.start & start <= crisp.reg.end) %>%
  rowid_to_column() %>%
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


# Plotting

p.crisp.genes <-  ggplot(crisp.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'CRISP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  geom_segment(aes(x=crisp.reg.start,xend=crisp.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('CRISPs') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(crisp.reg.start,crisp.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())

p.crisp.genes


# TADS plot
p.crisp.tads <- ggplot(crisp.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(crisp.reg.start,crisp.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.crisp.loops <- ggplot(crisp.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  geom_curve(data=crisp.CTCFloops,color='red',lwd=1) +
  coord_cartesian(xlim=c(crisp.reg.start,crisp.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.crisp.ctcf.chip <- ggplot(crisp.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = crisp.ctcf.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=crisp.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_point(inherit.aes = F,data=crisp.ctcf.bound,aes(x=(start+end)/2),y=crisp.ctcf.regionMax$max*1.1,pch=18,size=2.5) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(crisp.reg.start,crisp.reg.end),ylim=c(0,crisp.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.crisp.scale <- ggplot(crisp.info,aes(x=crisp.reg.start,xend=crisp.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=crisp.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=crisp.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(crisp.reg.start,crisp.reg.end)),label=crisp.reg.length)) +
  coord_cartesian(expand = T,xlim=c(crisp.reg.start,crisp.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

crisp.fullPlots <- p.crisp.genes  /  p.crisp.tads / p.crisp.ctcf.chip / p.crisp.loops / p.crisp.scale + plot_layout(ncol = 1,heights = c(1,1,1,4,1))
crisp.fullPlots


# CTL ---------------------------------------------------------------------


CTL.info <- all_info %>% 
  filter(str_detect(gene,'CTL')) 
CTL.reg.start <- min(CTL.info$start)-120000
CTL.reg.end <- max(CTL.info$end)+120000
CTL.reg.length <- paste(c(round((CTL.reg.end-CTL.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  CTCF chipseq and related files
CTL.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/1_047U_00JKUTA_VenomGland-495_CTCF_CroVir_i87_uniq_signal.bed',
                           col_names = c('chr','start','end','id','density')) %>% 
  filter(chr=='scaffold-mi5') %>% 
  filter(end >= CTL.reg.start & start <= CTL.reg.end) 

CTL.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                 col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi5') %>% 
  filter(end >= CTL.reg.start & start <= CTL.reg.end) 

CTL.ctcf.bound <- read_tsv('../../../../contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi5') %>% 
  filter(end >= CTL.reg.start & start <= CTL.reg.end) 

CTL.ctcf.regionMax <- CTL.ctcf.chip %>% 
  filter(end >= CTL.reg.start & start <= CTL.reg.end) %>% 
  summarise(max=max(density))

# Read in loops and TADs
CTL.CTCFloops <- read_tsv('../../../../contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-mi5')

CTL.loops <- read_tsv('../../../../contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-mi5')

CTL.TADs <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>%
  filter(chr=='scaffold-mi5') %>%
  filter(end >= CTL.reg.start & start <= CTL.reg.end) %>%
  rowid_to_column() %>%
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


# Plotting

p.CTL.genes <-  ggplot(CTL.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'CTL')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  geom_segment(aes(x=CTL.reg.start,xend=CTL.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('C-Type Lectin') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(CTL.reg.start,CTL.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())

p.CTL.genes


# TADS plot
p.CTL.tads <- ggplot(CTL.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(CTL.reg.start,CTL.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.CTL.loops <- ggplot(CTL.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  geom_curve(data=CTL.CTCFloops,color='red',lwd=1) +
  coord_cartesian(xlim=c(CTL.reg.start,CTL.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.CTL.ctcf.chip <- ggplot(CTL.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = CTL.ctcf.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=CTL.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_point(inherit.aes = F,data=CTL.ctcf.bound,aes(x=(start+end)/2),y=CTL.ctcf.regionMax$max*1.1,pch=18,size=2.5) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(CTL.reg.start,CTL.reg.end),ylim=c(0,CTL.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.CTL.scale <- ggplot(CTL.info,aes(x=CTL.reg.start,xend=CTL.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=CTL.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=CTL.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(CTL.reg.start,CTL.reg.end)),label=CTL.reg.length)) +
  coord_cartesian(expand = T,xlim=c(CTL.reg.start,CTL.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

CTL.fullPlots <- p.CTL.genes  /  p.CTL.tads / p.CTL.ctcf.chip / p.CTL.loops / p.CTL.scale + plot_layout(ncol = 1,heights = c(1,1,1,4,1))
CTL.fullPlots

# Exonuclease 1 -----------------------------------------------------------

exo1.info <- all_info %>% 
  filter(str_detect(gene,'Exonuclease 1')) 
exo1.reg.start <- min(exo1.info$start)-120000
exo1.reg.end <- max(exo1.info$end)+120000
exo1.reg.length <- paste(c(round((exo1.reg.end-exo1.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  CTCF chipseq and related files
exo1.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/1_047U_00JKUTA_VenomGland-495_CTCF_CroVir_i87_uniq_signal.bed',
                           col_names = c('chr','start','end','id','density')) %>% 
  filter(chr=='scaffold-ma6') %>% 
  filter(end >= exo1.reg.start & start <= exo1.reg.end) 

exo1.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                 col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-ma6') %>% 
  filter(end >= exo1.reg.start & start <= exo1.reg.end) 

exo1.ctcf.bound <- read_tsv('../../../../contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-ma6') %>% 
  filter(end >= exo1.reg.start & start <= exo1.reg.end) 

exo1.ctcf.regionMax <- exo1.ctcf.chip %>% 
  filter(end >= exo1.reg.start & start <= exo1.reg.end) %>% 
  summarise(max=max(density))

# Read in loops and TADs
exo1.CTCFloops <- read_tsv('../../../../contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-ma6')

exo1.loops <- read_tsv('../../../../contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-ma6')

exo1.TADs <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>%
  filter(chr=='scaffold-ma6') %>%
  filter(end >= exo1.reg.start & start <= exo1.reg.end) %>%
  rowid_to_column() %>%
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


# Plotting

p.exo1.genes <-  ggplot(exo1.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'Exonuclease 1')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  geom_segment(aes(x=exo1.reg.start,xend=exo1.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('Exonuclease 1') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(exo1.reg.start,exo1.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())

p.exo1.genes


# TADS plot
p.exo1.tads <- ggplot(exo1.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(exo1.reg.start,exo1.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.exo1.loops <- ggplot(exo1.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  geom_curve(data=exo1.CTCFloops,color='red',lwd=1) +
  coord_cartesian(xlim=c(exo1.reg.start,exo1.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.exo1.ctcf.chip <- ggplot(exo1.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = exo1.ctcf.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=exo1.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_point(inherit.aes = F,data=exo1.ctcf.bound,aes(x=(start+end)/2),y=exo1.ctcf.regionMax$max*1.1,pch=18,size=2.5) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(exo1.reg.start,exo1.reg.end),ylim=c(0,exo1.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.exo1.scale <- ggplot(exo1.info,aes(x=exo1.reg.start,xend=exo1.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=exo1.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=exo1.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(exo1.reg.start,exo1.reg.end)),label=exo1.reg.length)) +
  coord_cartesian(expand = T,xlim=c(exo1.reg.start,exo1.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

exo1.fullPlots <- p.exo1.genes  /  p.exo1.tads / p.exo1.ctcf.chip / p.exo1.loops / p.exo1.scale + plot_layout(ncol = 1,heights = c(1,1,1,4,1))
exo1.fullPlots


# Exonuclease 3 -----------------------------------------------------------

exo3.info <- all_info %>% 
  filter(str_detect(gene,'Exonuclease 3')) 
exo3.reg.start <- min(exo3.info$start)-120000
exo3.reg.end <- max(exo3.info$end)+120000
exo3.reg.length <- paste(c(round((exo3.reg.end-exo3.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in  CTCF chipseq and related files
exo3.ctcf.chip <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/bw2bed/bed/1_047U_00JKUTA_VenomGland-495_CTCF_CroVir_i87_uniq_signal.bed',
                           col_names = c('chr','start','end','id','density')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= exo3.reg.start & start <= exo3.reg.end) 

exo3.ctcf.chip.peaks <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/ChIP-seq/AM_UTA_CTCF_H3K27Ac_H3K4me3_CSeq_34205/intervals_txt/CTCF_VenomGland_intervals_simple.bed',
                                 col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= exo3.reg.start & start <= exo3.reg.end) 

exo3.ctcf.bound <- read_tsv('../../../../contact_loops/CTCF_BoundMotifs_02.03.21.bed',
                            col_names = c('chr','start','end','id')) %>% 
  filter(chr=='scaffold-mi7') %>% 
  filter(end >= exo3.reg.start & start <= exo3.reg.end) 

exo3.ctcf.regionMax <- exo3.ctcf.chip %>% 
  filter(end >= exo3.reg.start & start <= exo3.reg.end) %>% 
  summarise(max=max(density))

# Read in loops and TADs
exo3.CTCFloops <- read_tsv('../../../../contact_loops/CTCFLoops_merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                           col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-mi7')

exo3.loops <- read_tsv('../../../../contact_loops/merged_loops_loweracse_simple_10kbGreyArea.bedpe',
                       col_names = c('chrA','startA','endA','chrB','startB','endB','na1','na2')) %>%
  select(-7,-8) %>%
  filter(chrA=='scaffold-mi7')

exo3.TADs <- read_tsv('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/tad_identification/TADs_01.18.19/10kb/BED/allChrom_10kb_HierTADs_hierBOTH_01.22.19.bed',
                      col_names = c('chr','start','end','id','level')) %>%
  filter(chr=='scaffold-mi7') %>%
  filter(end >= exo3.reg.start & start <= exo3.reg.end) %>%
  rowid_to_column() %>%
  mutate(plot_y = ifelse(rowid%%2==0,1,-1))


# Plotting

p.exo3.genes <-  ggplot(exo3.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10avg1DPE)) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'Exonuclease 3')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  geom_segment(aes(x=exo3.reg.start,xend=exo3.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  ggtitle('Exonuclease 3') +
  scale_fill_viridis_c(option = 'B') +
  coord_cartesian(xlim=c(exo3.reg.start,exo3.reg.end),expand = T) +
  scale_x_continuous(labels = comma)+
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14),axis.title.x=element_blank())

p.exo3.genes


# TADS plot
p.exo3.tads <- ggplot(exo3.TADs,aes(x=start,xend=end,y=plot_y,yend=plot_y,color=plot_y)) +
  geom_segment(size=4,show.legend = F) +
  scale_color_continuous_sequential('OrRd',begin = 0.65) +
  coord_cartesian(xlim=c(exo3.reg.start,exo3.reg.end),ylim = c(-10,10),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('TADs') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

p.exo3.loops <- ggplot(exo3.loops,aes(x=(startA + endA)/2,xend=(startB + endB)/2,y=0,yend=0)) +
  # geom_segment(size=8,show.legend = F) +
  geom_curve() +
  geom_curve(data=exo3.CTCFloops,color='red',lwd=1) +
  coord_cartesian(xlim=c(exo3.reg.start,exo3.reg.end),ylim=c(-20,0),expand = T) +
  scale_x_continuous(labels = comma)+
  ggtitle('Loops') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',
                                        plot.title = element_text(color='black',face='bold',size = 14),
                                        axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


p.exo3.ctcf.chip <- ggplot(exo3.ctcf.chip,aes(x=start,xend=start,y=0,yend=density)) +
  geom_segment(data = exo3.ctcf.chip.peaks,aes(x=(start+end)/2,xend=(start+end)/2,y=0,yend=exo3.ctcf.regionMax$max*1.2),lwd=0.5,alpha=0.2) +
  geom_point(inherit.aes = F,data=exo3.ctcf.bound,aes(x=(start+end)/2),y=exo3.ctcf.regionMax$max*1.1,pch=18,size=2.5) +
  geom_segment(color='purple4') +
  coord_cartesian(expand = T,xlim=c(exo3.reg.start,exo3.reg.end),ylim=c(0,exo3.ctcf.regionMax$max*1.2)) +
  scale_x_continuous(labels = comma) +
  ggtitle('Chip-seq (CTCF)') +
  ylab('Read\nDensity') +
  theme_classic(base_size = 14) + theme(axis.title.x = element_blank(),plot.title.position = 'plot',plot.title = element_text(color='black',face='bold',size = 14))

p.exo3.scale <- ggplot(exo3.info,aes(x=exo3.reg.start,xend=exo3.reg.end,y=0,yend=0)) +
  geom_segment() +
  geom_point(aes(x=exo3.reg.start,y=0),pch='|',size=5)+
  geom_point(aes(x=exo3.reg.end,y=0),pch='|',size=5)+
  geom_label(aes(x=mean(c(exo3.reg.start,exo3.reg.end)),label=exo3.reg.length)) +
  coord_cartesian(expand = T,xlim=c(exo3.reg.start,exo3.reg.end)) +
  scale_x_continuous(labels = comma) +
  theme_void()

exo3.fullPlots <- p.exo3.genes  /  p.exo3.tads / p.exo3.ctcf.chip / p.exo3.loops / p.exo3.scale + plot_layout(ncol = 1,heights = c(1,1,1,4,1))
exo3.fullPlots
