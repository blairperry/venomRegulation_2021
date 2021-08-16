
library(tidyverse)

# Need to get avg median TPM counts at 1DPE and also fill in genes with no norm counts (i.e., filtered out for low expression in DESeq analysis) with expression = 0 so that ABC considers them "not expressed"

# Read in normalized counts
all.normcounts <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  select(txid,Median1DPE)

# Read in full gene annotation
all.annot <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>%
  select(X9) %>%
  filter(str_detect(X9,'trna',negate = T),str_detect(X9,'scaffold-un',negate = T)) %>%
  mutate(txid = str_split_fixed(X9,'[;]',4)[,3] %>% str_remove_all('Crovir_Transcript_ID=')) %>%
  select(txid)

all.1DPEavg.abc <- all.annot %>% 
  left_join(all.normcounts) %>% 
  mutate(Median1DPE = ifelse(is.na(Median1DPE),0,Median1DPE))

# write_tsv(all.1DPEavg.abc,'./analysis/6_ABC_Enhancers/misc/allGenes_1DPE_medianTPM_forABC_08.10.21.txt',col_names = F)


all.annot.full <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>%
  filter(str_detect(X9,'trna',negate = T),str_detect(X9,'scaffold-un',negate = T)) %>%
  mutate(txid = str_split_fixed(X9,'[;]',4)[,3] %>% str_remove_all('Crovir_Transcript_ID='))%>% 
  select(txid,X1, X4, X5)

all.1DPEavg.bedgraph <- all.normcounts %>% 
  left_join(all.annot.full) %>% 
  select(3,4,5,2,1) %>% 
  filter(!is.na(X4))

# write_tsv(all.1DPEavg.bedgraph,'./analysis/6_ABC_Enhancers/misc/allGenes_1DPE_MedianTPM_08.12.21.bedgraph',col_names = F)
