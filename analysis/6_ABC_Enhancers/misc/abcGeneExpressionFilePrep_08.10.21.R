
library(tidyverse)

# Need to get avg counts at 1DPE and also fill in genes with no norm counts (i.e., filtered out for low expression in DESeq analysis) with expression = 0 so that ABC considers them "not expressed"

# Read in normalized counts
all.normcounts <- read_csv('./analysis/1_gene_expression/norm_counts/CvvVenomReg_RepRNAseq_wNonVen_VSTNormCounts_08.02.21.csv')

# Read in full gene annotation
all.annot <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
  select(X9) %>% 
  filter(str_detect(X9,'trna',negate = T),str_detect(X9,'scaffold-un',negate = T)) %>% 
  mutate(txid = str_split_fixed(X9,'[;]',4)[,3] %>% str_remove_all('Crovir_Transcript_ID=')) %>% 
  select(txid)

all.1DPEavg.abc_temp <- all.normcounts %>% 
  mutate(txid = str_split_fixed(X1,'[_]',2)[,2]) %>% 
  select(gene=X1,4,8,9,2,3,5,6,7,10) %>% 
  mutate(gene = str_split_fixed(gene,'[_]',2)[,2]) %>% 
  mutate(avg1DPE = rowMeans(.[,2:4])) %>% 
  select(gene,avg1DPE)

all.1DPEavg.abc <- all.annot %>% 
  left_join(all.1DPEavg.abc_temp,by=c('txid'='gene')) %>% 
  mutate(avg1DPE = ifelse(is.na(avg1DPE),0,avg1DPE))

write_tsv(all.1DPEavg.abc,'./analysis/6_ABC_analysis/misc/allGenes_1DPEAvgExpression_08.10.21.txt',col_names = F)
