
library(tidyverse)
library(pheatmap)
library(viridis)

##
## Analyses to identify candidate TFs involved in the regulation of venom genes. 
##

'%!in%' <- function(x,y)!('%in%'(x,y))


##

# Read in longID to transcript ID table

longID.to.txID <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.annotation_table.final.txt',col_names = F) %>% 
  select(longid = 1, txid = 2) %>% 
  mutate(longid = str_remove_all(longid, '-mRNA-1'))


# Read in normalized counts

all.normCounts <- read_csv('./analysis/1_gene_expression/norm_counts/CvvVenomReg_RepRNAseq_wNonVen_VSTNormCounts_08.02.21.csv') %>% 
  mutate(txid = str_split_fixed(X1,'_',2)[,2]) %>% 
  mutate(id = str_split_fixed(X1,'_',2)[,1]) %>% 
  mutate(id = ifelse(id=='Venom',str_split_fixed(X1,'-',3)[,3],id))


# Read in pairwise result tables from DEseq2 and filter for significantly upregulated genes

Ven.vs.NonVen.pw <- read_csv('./analysis/1_gene_expression/pairwise_results/cvv_Venom.vs.NonVenom_VG_PairwiseResult_08.02.21.csv') %>% 
  filter(IHW_pvalue < 0.05) %>% 
  filter(IHW_pvalue < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, 'Up','Down')) %>%
  mutate(txid = str_split_fixed(X1,'_',2)[,2]) %>%
  mutate(id = str_split_fixed(X1,'_',2)[,1]) %>%
  mutate(id = ifelse(id=='Venom',str_split_fixed(X1,'-',3)[,3],id))

Ven.vs.NonVen.up <- Ven.vs.NonVen.pw %>% filter(direction == 'Up' & log2FoldChange > 1)
Ven.vs.NonVen.down <- Ven.vs.NonVen.pw %>% filter(direction == 'Down' & log2FoldChange < -1)



# Read in TF lists

tf.dnaBinding <- read_tsv('./data/annotation/TF_lists/Cvv_TF.DNABinding.txt',col_names = F) %>% 
  mutate(type = 'DNA Binding') %>% 
  select(longid=1,id=3, type)

tf.protBinding <- read_tsv('./data/annotation/TF_lists/Cvv_TF.ProtBinding.txt',col_names = F) %>% 
  mutate(type = 'Protein Binding') %>% 
  select(longid=1,id=2, type)

tf.coreg <- read_tsv('./data/annotation/TF_lists/Cvv_TF.Coregulators.txt',col_names = F) %>% 
  mutate(type = 'Coregulator') %>% 
  select(longid=1,id=2, type)

# Combine all TF tables and add T/F columns for the three types
tf.all <- tf.dnaBinding %>% 
  bind_rows(tf.protBinding) %>% 
  bind_rows(tf.coreg) %>% 
  left_join(longID.to.txID) %>% 
  select(longid,txid,id) %>% 
  unique() %>% 
  mutate(DnaBinding = ifelse(longid %in% tf.dnaBinding$longid,T,F)) %>% 
  mutate(ProtBinding = ifelse(longid %in% tf.protBinding$longid,T,F)) %>% 
  mutate(Coregulator = ifelse(longid %in% tf.coreg$longid,T,F))

# write_tsv(tf.all,'allTFs_withTXIDs_03.21.21.txt')
#write_tsv(as.data.frame(tf.all$id),'allTFs_forBackground_09.24.20.txt')


# Read in SE-associated genes

se.genes <- read_tsv('./analysis/2_super_enhancers/SEtoGene_Association/Cvv_SEassocPairs_Genes.bed',col_names = F) %>%
  mutate(longid = str_split_fixed(X4,'=',3)[,2]) %>%
  mutate(longid = str_remove_all(longid,';Name')) %>%
  left_join(longID.to.txID)


# Adding columns denoting if/where TFs are upregulated, if they are SE-associated, and TF type

tf.normCounts <- all.normCounts %>% 
  filter(txid %in% tf.all$txid) %>% 
  mutate(Upreg.Ven.vs.NonVen = ifelse(txid %in% Ven.vs.NonVen.up$txid,T,F)) %>% 
  mutate(SEassociated = ifelse(txid %in% se.genes$txid,T,F)) %>%
  left_join(tf.all,by='txid') %>% 
  select(-longid) %>% 
  unique()

tf.reformat <- tf.normCounts %>% 
  select(X1,id.y) %>% 
  mutate(cvvID = str_split_fixed(X1,'_',2)[,2]) %>% 
  select(cvvID,tfID = id.y)

# write_tsv(tf.reformat,'../tfbs_network_analysis/input/cvvTFIDs_01.28.21.tsv')
#write_csv(tf.normCounts,'UpregTFs_CountsAndInfo_09.23.20.csv')



# Upregulated TFs
# Filter details: upregulated in 1DPE compared to Unextracted

upreg.TFs <- tf.normCounts %>% filter(Upreg.Ven.vs.NonVen)

upreg.TFs.pw <- Ven.vs.NonVen.up %>% filter(X1 %in% upreg.TFs$X1)
# 111 upreg TFs


write_csv(upreg.TFs,'./analysis/3_transcription_factors/UpregulatedTFs_08.02.21.csv')

upreg.TFs.heatdata <- as.data.frame(upreg.TFs[,c(4,8,9,2,3,5,6,7,10)])
row.names(upreg.TFs.heatdata) <- upreg.TFs$id.y

upreg.TFs.annot <- as.data.frame(upreg.TFs[,c(14,16:18)]*1)
row.names(upreg.TFs.annot) <- upreg.TFs$id.y

pheatmap(upreg.TFs.heatdata,
         cluster_cols = F,
         scale='none',
         color=viridis(50),
         border_color = NA,
         treeheight_row = 10,
         cellwidth = 7,
         cellheight = 7,
         main = 'Upregulated TFs',
         fontsize = 6,
         annotation_row = upreg.TFs.annot,
         gaps_col = 3,
         # filename='./figures/fig_pieces/Upreg_TFheatmap_02.12.21.pdf'
)





# SE-associated TFs
# Filter Details: SE-associated 

seAssoc.TFs <- tf.normCounts %>% filter(SEassociated)

#write_csv(seAssoc.TFs,'./SEassociatedTFs_09.25.20.csv')

seAssoc.TFs.heatdata <- as.data.frame(seAssoc.TFs[,c(16,17,18,19,2,3,4,5,12,13,14,15,6,7,8,9,10,11)])
row.names(seAssoc.TFs.heatdata) <- seAssoc.TFs$id.y

seAssoc.TFs.annot <- as.data.frame(seAssoc.TFs[,28:26]*1)
row.names(seAssoc.TFs.annot) <- seAssoc.TFs$id.y

pheatmap(seAssoc.TFs.heatdata,
         cluster_cols = F,
         scale='row',
         color=viridis(50),
         border_color = NA,
         treeheight_row = 10,
         cellwidth = 10,
         cellheight = 10,
         main = 'SE-Associated TFs',
         annotation_row = seAssoc.TFs.annot,
         gaps_col = 12,
         # filename='./figures/fig_pieces/SEassoc_TFheatmap_02.12.21.pdf'
)


# Upregulated AND SE-associated TFs
# Filter Details: Overlap of above datasets

test <- (union(upreg.TFs$X1,seAssoc.TFs$X1))

both.upreg.SE.TFs <- tf.normCounts %>% filter(Upreg.Unext.vs.ODPE & SEassociated)

both.upreg.SE.TFs.heatdata <- as.data.frame(both.upreg.SE.TFs[,c(16,17,18,19,2,3,4,5,12,13,14,15,6,7,8,9,10,11)])
row.names(both.upreg.SE.TFs.heatdata) <- both.upreg.SE.TFs$id.y

both.upreg.SE.TFs.annot <- as.data.frame(both.upreg.SE.TFs[,28:26]*1)
row.names(both.upreg.SE.TFs.annot) <- both.upreg.SE.TFs$id.y

pheatmap(both.upreg.SE.TFs.heatdata,
         cluster_cols = F,
         scale='row',
         color=viridis(50),
         border_color = NA,
         treeheight_row = 10,
         cellwidth = 10,
         cellheight = 10,
         main = 'Upregulated AND SE-Associated TFs',
         annotation_row = both.upreg.SE.TFs.annot,
         gaps_col = 12,
         filename='./figures/fig_pieces/bothUpregAndSE_TFheatmap_09.23.20.pdf'
)

