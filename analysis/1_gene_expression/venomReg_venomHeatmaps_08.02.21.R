
library(pheatmap)
library(viridis)
library(stringr)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(ggrepel)
library(tidyverse)
library(gggenes)
library(patchwork)
library(readxl)


# Read in normalized counts
all.normcounts <- read_csv('./analysis/1_gene_expression/norm_counts/CvvVenomReg_RepRNAseq_wNonVen_VSTNormCounts_08.02.21.csv')

# Read in venom annotation
ven.annot <- read_xlsx('./data/venom_annotation/fullVenomAnnotation_02.20.20.xlsx',col_names = c('chr','start','end','description','messy_name','txid')) %>% 
  filter(description != 'Myotoxin/crotamine')   # excluding this annotation because it is wrong
  
# Read in pairwise results
ven.vs.nonven <- read_csv('./analysis/1_gene_expression/pairwise_results/cvv_Venom.vs.NonVenom_VG_PairwiseResult_08.02.21.csv') %>% 
  filter(IHW_pvalue < 0.05) %>%
  filter(log2FoldChange > 0) %>%
  mutate(X1=`...1`) %>% 
  mutate(txid = str_split_fixed(X1,'_',2)[,2]) 
  
  
# Read in conversion table for simpler venom gene names

ven.nameConvert <- read_tsv('./data/venom_annotation/VenomIDConvert.txt',col_names = c('txid','symbol'))



# New expression plot

VGOnly_ven.vs.nonven <- ven.vs.nonven %>% 
  filter(txid %in% ven.annot$txid) %>% 
  left_join(ven.annot) %>% 
  left_join(ven.nameConvert) %>% 
  mutate(family = str_split_fixed(symbol,'[_]',2)[,1]) %>% 
  mutate(family = ifelse(str_detect(family,'SVSP|SVMP|PLA2'),family,'Other')) %>% 
  mutate(symbol = str_replace(symbol,'_',' ')) %>% 
  mutate(family = factor(family,levels = c('SVMP','SVSP','PLA2','Other')))

ggplot(VGOnly_ven.vs.nonven,aes(x=reorder(symbol,-log2FoldChange),y=log2FoldChange,fill=family)) +
  geom_bar(stat='identity',show.legend = F) +
  xlab('Venom Gene') + ylab('Log2FoldChange\n(1DPE Venom Gland vs. Non-Venom Tissues)')+
  facet_grid(cols = vars(family),scales = 'free_x',space = "free") +
  theme_linedraw() +
  scale_y_continuous(expand = c(0,0),limits=c(0,20)) +
  scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),panel.grid = element_blank())






# # Filter norm counts for venom genes
# vg.normcounts <- all.normcounts %>% 
#   mutate(txid = str_split_fixed(X1,'[_]',2)[,2]) %>% 
#   filter(txid %in% ven.annot$...6) %>% 
#   left_join(ven.nameConvert) %>% 
#   select(gene=X1,symbol,4,8,9,2,3,5,6,7,10)
# 
# vg.1DPEavg <- vg.normcounts %>% 
#   mutate(avg1DPE = rowMeans(vg.normcounts[,3:5])) %>% 
#   select(gene,avg1DPE)
# 


# write_tsv(vg.1DPEavg,'VG_1DPEAvgExpression_08.02.21.tsv')

# 
# # Split by venom gene family/group
# 
# pla2.normcounts <- vg.normcounts %>% filter(str_detect(gene,'PLA2'))
# svmp.normcounts <- vg.normcounts %>% filter(str_detect(gene,'SVMP'))
# svsp.normcounts <- vg.normcounts %>% filter(str_detect(gene,'SVSP'))
# other.normcounts <- vg.normcounts %>% filter(str_detect(gene,'PLA2|SVMP|SVSP',negate = T)) %>% 
#   filter(gene %in% ven.vs.nonven$X1)
# 
# 
# # Format for pheatmap and reorder by exp magnitude
# pla2.heatdata <- as.data.frame(pla2.normcounts[,c(-1,-2)])
# row.names(pla2.heatdata) <- pla2.normcounts$symbol
# pla2.heatdata <- pla2.heatdata[order(pla2.heatdata$ODPE_1,decreasing = T),]
# 
# pla2.annotdata <- as_tibble(row.names(pla2.heatdata)) %>% 
#   left_join(ven.nameConvert,by=c('value'='symbol')) %>% 
#   left_join(ven.vs.nonven,by=c('txid'))
# 
# 
# svmp.heatdata <- as.data.frame(svmp.normcounts[,c(-1,-2)])
# row.names(svmp.heatdata) <- svmp.normcounts$symbol
# svmp.heatdata <- svmp.heatdata[order(svmp.heatdata$ODPE_1,decreasing = T),]
# 
# svmp.annotdata <- as_tibble(row.names(svmp.heatdata)) %>% 
#   left_join(ven.nameConvert,by=c('value'='symbol')) %>% 
#   left_join(ven.vs.nonven,by=c('txid'))
# 
# svsp.heatdata <- as.data.frame(svsp.normcounts[,c(-1,-2)])
# row.names(svsp.heatdata) <- svsp.normcounts$symbol
# svsp.heatdata <- svsp.heatdata[order(svsp.heatdata$ODPE_1,decreasing = T),]
# 
# svsp.annotdata <- as_tibble(row.names(svsp.heatdata)) %>% 
#   left_join(ven.nameConvert,by=c('value'='symbol')) %>% 
#   left_join(ven.vs.nonven,by=c('txid'))
# 
# other.heatdata <- as.data.frame(other.normcounts[,c(-1,-2)])
# row.names(other.heatdata) <- other.normcounts$symbol
# other.heatdata <- other.heatdata[order(other.heatdata$ODPE_1,decreasing = T),]
# 
# # Plot heatmaps
# 
# breaksList = seq(min(vg.normcounts[,3:11]), 19, by = 0.01)
# 
# pheatmap(pla2.heatdata,
#          scale='none',
#          # clustering_distance_rows = 'manhattan',
#          cluster_cols = F,cluster_rows = F,
#          col=inferno(length(breaksList)),
#          breaks = breaksList,
#          show_rownames = T,
#          border_color = NA,
#          fontsize_row = 9,
#          cellwidth = 20,cellheight = 10,
#          treeheight_row = 0,
#          gaps_col = c(3))
# 
# pheatmap(svmp.heatdata,
#          scale='none',
#          # clustering_distance_rows = 'manhattan',
#          cluster_cols = F,cluster_rows = F,
#          col=inferno(length(breaksList)),
#          breaks = breaksList,
#          show_rownames = T,
#          border_color = NA,
#          fontsize_row = 9,
#          cellwidth = 20,cellheight = 10,
#          treeheight_row = 0,
#          gaps_col = 3)
# 
# pheatmap(svsp.heatdata,
#          scale='none',
#          # clustering_distance_rows = 'manhattan',
#          cluster_cols = F,cluster_rows = F,
#          col=inferno(length(breaksList)),
#          breaks = breaksList,
#          show_rownames = T,
#          border_color = NA,
#          fontsize_row = 9,
#          cellwidth = 20,cellheight = 10,
#          treeheight_row = 0,
#          gaps_col = 3)
# 
# pheatmap(other.heatdata,
#          scale='none',
#          cluster_cols = F,cluster_rows = F,
#          col=inferno(length(breaksList)),
#          breaks = breaksList,
#          show_rownames = T,
#          border_color = NA,
#          fontsize_row = 9,
#          cellwidth = 20,cellheight = 10,
#          treeheight_row = 0,
#          gaps_col = 3)


###
### NOTE: "Priority" list is comprised of VGs from three main families, their NVPs, and "Other" VGs that are upregulated in venom gland tissues
###       Here, I am adding the second CRISP in so that we can vizualize that cluster in Figure 1. 



### Gene cluster plotting w/ TPM Normed counts

CRISP2 <- c('scaffold-ma1',	'169434958',	'169437996',	'Cysteine-rich secretory protein',	'maker-scaffold-ma1-augustus-gene-564.2',	'crovir-transcript-8572',	'CRISP_2')

pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F) %>% 
  add_row(X1='scaffold-ma1',X2=169434958,X3=169437996,X4='Cysteine-rich secretory protein',X5='maker-scaffold-ma1-augustus-gene-564.2',X6='crovir-transcript-8572',X7='CRISP_2')

tpm.exp <- read_tsv('analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  select(txid,Median1DPE) %>% 
  filter(txid %in% pri_venom_genes$X6) %>% 
  left_join(ven.nameConvert)

max(log10(tpm.exp$Median1DPE+1))


all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
  filter(str_detect(X9,'trnascan',negate = T)) %>% 
  mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>% 
  mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>% 
  filter(tx_id %in% pri_venom_genes$X6) %>% 
  left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>% 
  select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>% 
  mutate(strand = ifelse(strand == '+','forward','reverse')) %>% 
  mutate(direction = ifelse(strand == 'forward',1,-1)) %>% 
  left_join(tpm.exp,by=c('tx_id'='txid')) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'ADAM28',negate = T),str_replace_all(symbol,'_',' '),symbol)) %>% 
  mutate(symbol = ifelse(str_detect(symbol,'ADAM28|gIIE'),paste('NVP: ',symbol,sep = ''),symbol))




SVMP <- all_info %>% 
  filter(str_detect(symbol,'SVMP|ADAM')) %>% 
  ggplot(aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(symbol,'SVMP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = molecule, label = symbol), inherit.aes = F, nudge_y = 1) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(symbol,'ADAM')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = molecule, label = symbol), inherit.aes = F, nudge_y = 0.25,nudge_x = 100000) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),show.legend = T) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = scales::comma,limits=c(13890000,14495000)) +
  theme_genes()
SVMP

# region length (kb)
(14495000 - 13890000) / 1000


SVSP <- all_info %>% 
  filter(str_detect(symbol,'SVSP')) %>% 
  ggplot(aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(symbol,'SVSP')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = molecule, label = symbol), inherit.aes = F, nudge_y = 1) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = scales::comma,limits=c(8560000,9000000)) +
  theme_genes()
SVSP

# region length (kb)
(9000000 - 8560000) / 1000


PLA2 <- all_info %>% 
  filter(str_detect(symbol,'PLA2')) %>% 
  ggplot(aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(symbol,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = molecule, label = symbol), inherit.aes = F, nudge_y = 1) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = scales::comma,limits=c(3019000,3045000)) +
  theme_genes()
PLA2

# region length (kb)
(3045000 - 3019000) / 1000

CRISP <- all_info %>% 
  filter(str_detect(symbol,'CRISP')) %>% 
  ggplot(aes(xmin = start, xmax = end, y = molecule, forward = direction, fill = log10(Median1DPE+1))) +
  ggrepel::geom_text_repel(data = all_info %>%
                          filter(str_detect(symbol,'CRISP')) %>%
                          mutate(start = (start + end)/2),
                          aes(x = start, y = molecule, label = symbol), inherit.aes = F, nudge_y = 1) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = scales::comma,limits=c(169423000,169439000)) +
  theme_genes()
CRISP

SVMP / SVSP / PLA2 / CRISP 

