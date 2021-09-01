
library(tidyverse)
library(scico)


# Read in all candidate TFs -----------------------------------------------

cand_tfs <- read_tsv('analysis/3_transcription_factors/allCandidateTFs_08.03.21.tsv') 


# Read in candidate subset (functional categories) ------------------------

funcChar_tfs <- read_tsv('analysis/3_transcription_factors/_functional_characterization/CandTFs_InterestingSubset_FULL_08.27.21.tsv') %>% 
  pivot_wider(names_from = feature,values_from = value)


# Read in pairwise expression results and normalized counts  ------------------------------------

ven.vs.nonven <- read_csv('analysis/1_gene_expression/pairwise_results/cvv_Venom.vs.NonVenom_VG_PairwiseResult_08.02.21.csv') %>% 
  mutate(txid = str_split_fixed(X1,'_',2)[,2]) %>% 
  select(txid,
         `Non-Venom vs. Venom - Log2FoldChange`=3,
         `Non-Venom vs. Venom - p-value`=6,
         `Non-Venom vs. Venom - IHW p-value`=8)



# Build table -------------------------------------------------------------

cand_tfs.annot.wExp <- cand_tfs %>% 
  left_join(funcChar_tfs,by=c('id.y'='id')) %>% 
  left_join(ven.vs.nonven,by=c('txid'))


# Export table for formatting in Excel ------------------------------------

write_csv(cand_tfs.annot.wExp,'analysis/3_transcription_factors/CandTFInfo_forSupp_08.31.21.csv')
