
library(tidyverse)
library(ggvenn)

svmp.enrich <- read_csv('analysis/5_promoters/_tfbs_analysis/SVMP.vs.NonVen/Enrichment: Text4550150718958082102_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')

svsp.enrich <- read_csv('analysis/5_promoters/_tfbs_analysis/SVSP.vs.NonVen/Enrichment: Text8267836149309786442_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')

pla2.enrich <- read_csv('analysis/5_promoters/_tfbs_analysis/PLA2.vs.NonVen/Enrichment: Text3896374089949899309_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')

forVenn <- list('SVMP'=svmp.enrich$transcription_factor_name,
                  'SVSP'=svsp.enrich$transcription_factor_name,
                  'PLA2'=pla2.enrich$transcription_factor_name)

ggvenn(forVenn)
