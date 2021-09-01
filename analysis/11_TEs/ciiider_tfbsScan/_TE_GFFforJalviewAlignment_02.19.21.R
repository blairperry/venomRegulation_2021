
library(tidyverse)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/_tes_giulia/_CroMit_hAT-Tip_100_NEW_02.13.21/alignments")

ciiider.res <- read_csv('ciiider_TFBSscan/Promoter Panel.csv') %>% 
  janitor::clean_names() %>% 
  mutate(gene_name = str_replace_all(gene_name,'[()]','_'))

te.enh.cand_tfs <- read_tsv('../../../transcription_factors/_tfbs_analysis/_Enhancer_CIIIDER_TFBS_EnrichmentAnalyses_12.11.20/All3_PrimSecHighTFs_Enh_01.28.21.tsv') %>% 
  filter(group == 'SVSP' & round == 'Primary') %>% 
  select(tfbs,M_id) %>% 
  unique()

te.prom.cand_tfs <- read_tsv('../../../_new_promoterATACpeaks/_venom_promoterPeaks/_tfbs_enrichment_binding/All3_PrimSecHighTFs_01.27.21.tsv') %>% 
  filter(group == 'SVSP' & round == 'Primary') %>% 
  select(tfbs,M_id) %>% 
  unique()

te.all.cand_tfs <- te.enh.cand_tfs %>% 
  bind_rows(te.prom.cand_tfs) %>% 
  unique()


ciiider.res.gff <- ciiider.res %>% 
  filter(transcription_factor_id %in% te.all.cand_tfs$M_id) %>% 
  mutate(source = 'ciiider',
         type = 'tfbs',
         phase = '.',
         attrib = paste('ID=',transcription_factor_name,sep = ''),
         strand = ifelse(strand == 1,'+','-')) %>% 
  select(gene_name,source,type,start_position,end_position,matrix_match_score,strand,phase,attrib)
  
write_tsv(ciiider.res.gff,'ciiider_TFBSscan/ciiiderResult_PromEnhancerPrimaryTFBS.gff',col_names = F)
