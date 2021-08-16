
library(tidyverse)

raw.counts <- read_tsv('./data/rnaseq/raw_counts/cvv_VenomRegulationRNAseq_wNonVenom_rawCounts_08.02.21.txt',comment = '#') %>% 
  select(7,6,8:35) %>% 
  filter(!is.na(Crovir_Transcript_ID)) %>% 
  select(1,2,5,13,14,3,4,9,11,12,19)

raw.matrix <- as.matrix(raw.counts[,-1])
row.names(raw.matrix) <- raw.counts$Crovir_Transcript_ID
colnames(raw.matrix) <- str_remove_all(colnames(raw.matrix),'./STAR_mapped/') %>% str_remove_all(.,'Aligned.sortedByCoord.out.bam')

geneLengths <- as.vector(subset(raw.matrix, select = c(Length)))

rpk <- apply( subset(raw.matrix, select = c(-Length)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

colSums(tpm)

tpm.avg.median <- as_tibble(tpm) %>% 
  mutate(txid = row.names(tpm)) %>% 
  select(txid,1,2,3) %>% 
  mutate(Avg1DPE = rowMeans(.[,2:4])) %>% 
  rowwise() %>% 
  mutate(Median1DPE = median(ODPE_1,RVG_11,RVG_4))

write_tsv(tpm.avg.median,'./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv')
